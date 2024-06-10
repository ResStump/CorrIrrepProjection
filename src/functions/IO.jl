@doc raw"""
    Parms

Important parameters for the perambulator contractions.
"""
struct Parms
    # String containt in parameter toml file passed to program
    parms_toml_string::String

    # Paths
    raw_correlator_dir
    result_dir

    # Configuration numbers and source times
    cnfg_indices::Vector{Int}
    tsrc_arr::Array{Int, 2}

    # Lattice size in time
    Nₜ::Int

    # Number of configurations and sources
    N_cnfg::Int
    N_src::Int
end

parms = nothing
parms_toml = nothing

@doc raw"""
    read_parameters()

Read the parameters stored in the parameter file passed to the program with the flag -i and
return the dictonary `parms_toml` and the Parms instance `parms`.
"""
function read_parameters()
    # Search for parameter file in arguments passed to program
    parms_file_index = findfirst(arg -> arg == "-i", ARGS)
    if isnothing(parms_file_index)
        throw(ArgumentError("argument -i not provided to the program."))
    elseif parms_file_index == size(ARGS)[1]
        throw(ArgumentError("argument after -i not provided to the program."))
    end

    parms_file = ARGS[parms_file_index+1]

    # Read parameters from parameter file and store them in the `parms_toml`
    parms_toml_string = read(parms_file, String)
    global parms_toml = TOML.parse(parms_toml_string)

    # Add information about program to `parms_toml` (as string)
    parms_toml["Program Information"] =
        "Julia version = $VERSION\n"*
        "$(@__MODULE__) version = $(pkgversion(@__MODULE__))\n"*
        "Program file = $PROGRAM_FILE\n"

    # Read source times
    tsrc_list = DF.readdlm(
        parms_toml["Directories and Files"]["tsrc_list"], ' ', Int
    )
    
    # Set paths
    raw_correlator_dir = Path(parms_toml["Directories and Files"]["raw_correlator_dir"])
    result_dir = Path(parms_toml["Directories and Files"]["result_dir"])

    # Lattice size in time
    Nₜ = parms_toml["Geometry"]["N_t"]

    # Set `cnfg_indices`
    first_cnfg = parms_toml["Configurations"]["first"]
    step_cnfg = parms_toml["Configurations"]["step"]
    last_cnfg = parms_toml["Configurations"]["last"]
    N_cnfg = (last_cnfg - first_cnfg) ÷ step_cnfg + 1
    cnfg_indices = Array(first_cnfg:step_cnfg:last_cnfg)

    # Find first source times for specified configurations
    tsrc_first = Array{Int, 1}(undef, N_cnfg)
    for (i_cnfg, n_cnfg) in enumerate(cnfg_indices)
        idx = findfirst(n -> n == n_cnfg, tsrc_list[:, 1])
        tsrc_first[i_cnfg] = tsrc_list[idx, 2]
    end

    # Store all source times in `tsrc_arr`
    N_src = parms_toml["Sources"]["N_src"]
    src_separation = parms_toml["Sources"]["src_separation"]
    tsrc_arr = hcat([tsrc_first .+ i_src*src_separation
                     for i_src in 0:N_src-1]...)
    tsrc_arr = mod.(tsrc_arr, Nₜ) # periodically shift values >= Nₜ

    # Store all parameters
    global parms = Parms(parms_toml_string, raw_correlator_dir, result_dir, cnfg_indices,
                         tsrc_arr, Nₜ, N_cnfg, N_src)
    
    return
end

@doc raw"""
    write_correlator(correlator_file, correlator, dataset)

Write `correlator` and its dimension labels to the dataset `dataset` in the HDF5 file
`correlator_file`. Additionally, also write the parameter file `parms_toml_string` and
program information to it.
"""
function write_correlator(correlator_file, correlator, dataset)
    hdf5_file = HDF5.h5open(string(correlator_file), "w")

    # Write correlator with dimension labels
    hdf5_file[dataset] = correlator
    HDF5.attributes(hdf5_file[dataset])["DIMENSION_LABELS"] = ["t", "source", "cnfg"]

    # Write parameter file
    hdf5_file["parms.toml"] = parms.parms_toml_string

    # Write program information
    hdf5_file["Program Information"] = parms_toml["Program Information"]
    
    close(hdf5_file)

    return
end
