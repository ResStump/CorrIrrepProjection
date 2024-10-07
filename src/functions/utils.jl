@doc raw"""
    scalar_mul(O_arr, λ) -> O_arr_new

Multiply each operator in `O_arr` with the scalar `λ` and return the new array `O_arr_new`.
"""
function scalar_mul(O_arr, λ)
    O_arr_new = copy(O_arr)
    map(O -> O["N"] *= λ, O_arr_new)
    return O_arr_new
end

# Gamma matrix strings
γ = Dict(1=>"gamma_1", 2=>"gamma_2", 3=>"gamma_3", 4=>"gamma_4", 5=>"gamma_5")

@doc raw"""
    Γ_str_to_idx(Γ_str, Γ_labels) -> idx

Return the index of the gamma matrix string `Γ_str` in the gamma matrix labels `Γ_labels`.
"""
function Γ_str_to_idx(Γ_str, Γ_labels)
    idx = findfirst(isequal(Γ_str), Γ_labels)
    
    if isnothing(idx)
        throw(ArgumentError("Gamma matrix string not found in the labels."))
    end

    return idx
end