# %%########################################################################################
# DD_projections.jl
#
# Project DD correlators to irreps.
#
# Usage:
#   DD_projections.jl -i <parms file>
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
# Global Carameters
###################

# Set global parameters
CIP.read_parameters()


#%% #######################
# Allocate Correlator Array
###########################

correlator_size = CIP.parms.Nₜ, CIP.parms.N_src, CIP.parms.N_cnfg
correlator_I0 = Array{ComplexF64}(undef, correlator_size)
correlator_I1 = Array{ComplexF64}(undef, correlator_size)


# %%#########
# Calculation
#############

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
                raw_correlator_file = CIP.parms.raw_correlator_dir/
                    "$(CIP.parms_toml["File base name"]["name"])_n$(n_cnfg)_tsrc$(t₀).hdf5"
                
                correlator_I0[:, i_src, i_cnfg] = 
                    1/3*(CIP.DD_local_T₁⁺(raw_correlator_file, 1, 0) +
                         CIP.DD_local_T₁⁺(raw_correlator_file, 2, 0) +
                         CIP.DD_local_T₁⁺(raw_correlator_file, 3, 0))
                correlator_I1[:, i_src, i_cnfg] =
                    1/3*(CIP.DD_local_T₁⁺(raw_correlator_file, 1, 1) +
                         CIP.DD_local_T₁⁺(raw_correlator_file, 2, 1) +
                         CIP.DD_local_T₁⁺(raw_correlator_file, 3, 1))

                end

        end
        println("\n")
    end

    # Broadcast correlators to all ranks
    @time "Broadcast correlators" begin
        CIP.broadcast_correlators!(correlator_I0)
        CIP.broadcast_correlators!(correlator_I1)
    end

    if myrank == 0
        proj_correlator_file = CIP.parms.result_dir/
            "$(CIP.parms_toml["File base name"]["name"]).hdf5"

        CIP.write_correlator(proj_correlator_file, correlator_I0,
                             "c2_DD_local/I0/P000/T1+", "w")
        CIP.write_correlator(proj_correlator_file, correlator_I1,
                             "c2_DD_local/I1/P000/T1+", "r+")
    end

end

main()


# %%