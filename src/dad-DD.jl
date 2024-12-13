# %%########################################################################################
# dad-DD.jl
#
# Project correlators from diquark-antidiquark and DD operators to their irreducible
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
O_i_arr = [[DDstarₛᵢ_A₁⁺0_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_A₁⁺1_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_E⁺1_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_A₁⁺2_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_J1⁺2_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_J3⁺2_T₁⁺_I0_nonlocal(i),
            DstarDstarₛᵢ_A₁⁺0_S1_T₁⁺_I0_nonlocal(i),
            DstarDstarₛᵢ_A₁⁺1_S1_T₁⁺_I0_nonlocal(i),
            DstarDstarₛᵢ_E⁺1_S1_T₁⁺_I0_nonlocal(i),
            DDstarₛᵢ_A₁⁺_T₁⁺_I0_local(i),
            dadᵢ_A₁⁺_T₁⁺_I0_local(i),
            DstarDstarₛᵢ_A₁⁺_S1_T₁⁺_I0_local(i)]
            for i in 1:3]

N_op = length(O_i_arr[1])

# Operator labels
operator_labels = ["DD*s nonlocal A1+(0)", "DD*s nonlocal A1+(1)", "DD*s nonlocal E+(1)",
                   "DD*s nonlocal A1+(2)", "DD*s nonlocal J1+(2)", "DD*s nonlocal J3+(2)",
                   "D*D*s nonlocal A1+(0),1", "D*D*s nonlocal A1+(1),1",
                   "D*D*s nonlocal E+(1),1",
                   "DD*s local A1+", "dad local A1+", "D*D*s local A1+"]


# Skipped correlator matrix entries
skip_operators = CIP.parms_toml["Corr Matrix Entries"]["skip_operators"]
skipped_entries = [[i, j] for (i, j) in Iterators.product(skip_operators, skip_operators)]

# Correlator matrix entries from external file (if they are already computed)
ext_corr_types = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["types"]
ext_corr_file = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["file"]
ext_corr_group = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["group"]
ext_corr_operators = CIP.parms_toml["Corr Matrix Entries"]["external corr"]["operators"]
ext_corr_entries = [[i, j]
                    for (i, j) in Iterators.product(ext_corr_operators, ext_corr_operators)]



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
    if length(skipped_entries) != 0
        types_arr = filter(type -> type ∉ ext_corr_types, types_arr)
    end

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
    # Allocate corr matrix with zeros (forward and backward shape)
    Cₜ_fb = zeros(ComplexF64, CIP.parms.Nₜ÷2+1, N_op, N_op, 2)

    # Loop over spin index
    for i in 1:3
        # Loop over all operators
        for (i_O_sink, O_sink) in enumerate(O_i_arr[i])
            for (i_O_src, O_src) in enumerate(O_i_arr[i])
                # Check if this entrie should be computed
                if !([i_O_sink, i_O_src] in skipped_entries)
                    Cₜ_fb[:, i_O_sink, i_O_src, :] .+= 
                        1/3*CIP.project_tetraquark_corr(O_sink, O_src, raw_corr_dict)
                end
            end
        end
    end
    
    return Cₜ_fb
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
                Cₜ_fb = compute_corr_matrix_entries(raw_corr_dict)

                # Store correlator matrix entries (transpose backward correlator)
                corr_matrix[:, :, :, 1, i_src, i_cnfg] = Cₜ_fb[:, :, :, 1]
                corr_matrix[:, :, :, 2, i_src, i_cnfg] = 
                    permutedims(Cₜ_fb[:, :, :, 2], [1, 3, 2])
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
        group = "c2_dad-DD/I0/P000/T1+"
        CIP.write_corr_matrix(correlator_file, corr_matrix, group, operator_labels,
                              dimension_labels, "w")
    end

end

main()


# %%
