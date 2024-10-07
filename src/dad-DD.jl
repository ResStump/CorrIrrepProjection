# %%########################################################################################
# dad-DD.jl
#
# Project correlators from DD and diquark-antidiquark operators to their irreducible
# representations.
#
# Usage:
#   dad-DD.jl -i <parms file>
#
# where <parms file> is a toml file containing the required parameters
#
############################################################################################

import MKL
import LinearAlgebra as LA
import MPI
import HDF5
import DelimitedFiles as DF
import FilePathsBase: /, Path
import BenchmarkTools.@btime
import CorrIrrepProjection as CIP

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
myrank = MPI.Comm_rank(comm)
N_ranks = MPI.Comm_size(comm)

if myrank != 0
    redirect_stdout(devnull)
end


# %%###############
# Global Parameters
###################

# Set global parameters
CIP.read_parameters()


# %%#######
# Operators
###########

include("operators/DD_operators.jl")
include("operators/dad_operators.jl")

# Operator array
O_arr(i) = [DDstarₛᵢ_A₁⁺0_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_A₁⁺1_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_A₁⁺p0_T₁⁺_I0_local(i),
            dadᵢ_A₁⁺p0_T₁⁺_I0_local(i)]

N_op = length(O_arr(1))

# Operator labels
operator_labels = ["DD*s nonlocal A1+(0)", "DD*s nonlocal A1+(1)", "DD*s local A1+",
                   "dad local A1+"]


# Skipped correlator matrix entries
skipped_entries = CIP.parms_toml["Corr Matrix Entries"]["skip_entries"]

# Correlator matrix entries from external file (if they are already computed)
ext_corr_types = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["types"]
ext_corr_file = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["file"]
ext_corr_group = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["group"]
ext_corr_entries = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["entries"]



# %%###############################
# Paths and Functions to Read Files
###################################

# Path to result directory
result_dir = Path(CIP.parms_toml["Directories and Files"]["result_dir"])


function read_raw_corr(label, n_cnfg, t₀)
    base_name = CIP.parms_toml["File base names"][label]
    dir = Path(CIP.parms_toml["Directories and Files"]["raw_corr_$(label)_dir"])
    file_path = dir/"$(base_name)_n$(n_cnfg)_tsrc$(t₀).hdf5"

    file = HDF5.h5open(string(file_path), "r")
    correlators = read(file["Correlators"])
    spin_structure = read(file["Spin Structure"])
    close(file)

    return Dict("Correlators" => correlators, "Spin Structure" => spin_structure)
end

function get_raw_corr_dict(n_cnfg, t₀)
    raw_corr_dict = Dict()

    types_arr = ["DD_local", "DD_nonlocal", "DD_mixed", "dad_local",
                 "dad-DD_local", "dad-DD_local-nonlocal"]
    
    # Remove types that are skipped
    types_arr = filter(type -> type ∉ ext_corr_types, types_arr)

    for type in types_arr
        raw_corr_dict[type] = read_raw_corr(type, n_cnfg, t₀)
    end

    return raw_corr_dict
end


#%% ########################
# Allocate Correlator Arrays
############################

corr_matrix_size = CIP.parms.Nₜ÷2+1, N_op, N_op, 2, CIP.parms.N_src, CIP.parms.N_cnfg
corr_matrix = Array{ComplexF64}(undef, corr_matrix_size)

# Dimension labels (reversed order in julia)
dimension_labels = ["config", "tsrc", "fwd/bwd", "op_src", "op_snk", "t"]


# %%#########
# Calculation
#############

function compute_corr_matrix_entries(raw_corr_dict)
    # Allocate corr matrix with zeros
    Cₜ = zeros(ComplexF64, CIP.parms.Nₜ, N_op, N_op)

    # Loop over spin index
    for i in 1:3
        # Loop over all operators
        for (i_O_sink, O_sink) in enumerate(O_arr(i))
            for (i_O_src, O_src) in enumerate(O_arr(i))
                # Check if this entrie should be computed
                if !([i_O_sink, i_O_src] in skipped_entries)
                    Cₜ[:, i_O_sink, i_O_src] .+= 
                        1/3*CIP.project_tetraquark_corr(O_sink, O_src, raw_corr_dict)
                end
            end
        end
    end
    
    return Cₜ
end

function main()
    # Loop over all configurations
    for (i_cnfg, n_cnfg) in enumerate(CIP.parms.cnfg_indices)
        # Skip the cnfgs this rank doesn't have to compute
        if !CIP.is_my_cnfg(i_cnfg)
            continue
        end

        #println("Configuration $n_cnfg")
        @time "Finished configuration $n_cnfg" begin
            # Loop over all sources
            for (i_src, t₀) in enumerate(CIP.parms.tsrc_arr[i_cnfg, :])
                # println("  Source: $i_src of $(CIP.parms.N_src)")

                # Compute correlator matrix
                raw_corr_dict = get_raw_corr_dict(n_cnfg, t₀)
                Cₜ = compute_corr_matrix_entries(raw_corr_dict) 

                # Bring correlator in forward/backward shape and store it
                #########################################################
                Nₜ = CIP.parms.Nₜ

                # Forward correlator
                Cₜ_fwd = Cₜ[1:Nₜ÷2+1, :, :]
                corr_matrix[:, :, :, 1, i_src, i_cnfg] = Cₜ_fwd

                # Backward correlator (reversed time order and operator indices permuted)
                bwd_time_indices = [1, (Nₜ:-1:Nₜ÷2+1)...]
                Cₜ_bwd = Cₜ[bwd_time_indices, :, :]
                corr_matrix[:, :, :, 2, i_src, i_cnfg] = permutedims(Cₜ_bwd, [1, 3, 2])
            end
        end
        println("\n")
    end

    # Broadcast correlators to all ranks
    @time "Broadcast correlators" begin
        CIP.broadcast_correlators!(corr_matrix, cnfg_dim=6)
    end

    # Finalize and write correlator
    if myrank == 0
        # Take entries from external file
        if length(skipped_entries) != 0
            if length(skipped_entries) != length(ext_corr_entries)
                throw(ArgumentError("There must be as many entries that are skipped as " *
                                    "entries that are taken from the external file."))
            end

            # Read external correlator
            ext_corr = HDF5.h5read(ext_corr_file, ext_corr_group)

            for ((i_O_sink, i_O_src), (i_O_sink_ext, i_O_src_ext)) in 
                zip(skipped_entries, ext_corr_entries)

                # Store external correlator entries in correlator
                # (assume all sorces are used)
                corr_matrix[:, i_O_sink, i_O_src, :, :, :] = 
                    ext_corr[:, i_O_sink_ext, i_O_src_ext, :, :, CIP.parms.cnfg_indices]
            end
        end

        # Write correlator
        correlator_file = result_dir/CIP.parms_toml["File base names"]["corr_result"]
        group = "c2_DD/I0/P000/T1+"
        CIP.write_corr_matrix(correlator_file, corr_matrix, group, operator_labels,
                              dimension_labels, "w")
    end

end

main()


# %%
