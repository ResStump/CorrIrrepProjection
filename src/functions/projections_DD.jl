#= @doc raw"""
    doc ...
""" =#
function DD_local_T₁⁺(raw_corr_file, spin_idx::Int, isospin::Int)
    # Define Operators
    ##################
    if spin_idx == 1
        γᵢ = :γ₁
    elseif spin_idx == 2
        γᵢ = :γ₂
    elseif spin_idx == 3
        γᵢ = :γ₃
    else
        throw(ArgumentError("spin_idx has to be 1, 2, or 3."))
    end

    # Momentum
    p = [0, 0, 0]
    p_str = "p$(join(p, ","))"

    # Indices of the monomials of gamma matrices
    Γ_indices = Dict(:γ₅=>1, :γ₁=>2, :γ₂=>3, :γ₃=>4, :mi=>5)

    # Sink operator
    O_sink = [
        Dict("factor"=>1/√2, :Γ₁=>:γ₅, :Γ₂=>γᵢ),
        Dict("factor"=>1/√2, :Γ₁=>γᵢ, :Γ₂=>:γ₅)
    ]
    
    # Isospin
    if isospin == 0
        O_sink[2]["factor"] *= -1
    elseif  isospin != 1
        throw(ArgumentError("isospin has to be 0 or 1."))
    end

    # Source operator
    O_src = copy(O_sink)


    # Read correlators and project to irrep
    #######################################

    # Create Correlator
    Cₜ = zeros(ComplexF64, parms.Nₜ)

    hdf5_file = HDF5.h5open(string(raw_corr_file))

    for (O₁, O₂) in Iterators.product(O_sink, O_src)
        # Indices for gamma matrices
        Γ_idx = [Γ_indices[O₁[:Γ₁]], Γ_indices[O₁[:Γ₂]],
                 Γ_indices[O₂[:Γ₁]], Γ_indices[O₂[:Γ₂]]]

        Cₜ_ = hdf5_file["Correlator/$p_str"][:, Γ_idx...]
        Cₜ += O₁["factor"]*O₂["factor"]*Cₜ_
    end

    close(hdf5_file)

    return Cₜ
end