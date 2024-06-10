module CorrIrrepProjection

import MKL
import LinearAlgebra as LA
import MPI
import TOML
import HDF5
import DelimitedFiles as DF
import FilePathsBase: /, Path

include("functions/IO.jl")
include("functions/mpi_utils.jl")
include("functions/projections_DD.jl")

end # module CorrIrrepProjection
