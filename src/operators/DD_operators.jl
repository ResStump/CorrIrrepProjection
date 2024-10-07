γ = CIP.γ

# %%################
# Nonlocal Operators
####################

# Generic DD operator
#####################

DDstarₐ_I0_nonlocal(p₁, p₂) = begin
    @assert p₁ != p₂ "operator is zero for p₁ == p₂"
    O = [
        Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
             "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2),
        Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
             "Γ₁"=>γ[5], "Γ₂"=>γ[5],"N"=>-1/√2)
        ]
    return O
end

DDstarₛ_I1_nonlocal(p₁, p₂) = begin
    if p₁ == p₂
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2)
        ]
    else
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2)
        ]
    end
    return O
end


DDstarₛᵢ_I0_nonlocal(i, p₁, p₂) = begin
    if p₁ == p₂
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/√2)
        ]
    else
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/2),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/2),
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>-1/2),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>1/2)
        ]
    end
    return O
end

# T₁⁺ operators
###############

# Angular momentum: A₁⁺, Ptot²=0; spin: T₁⁺; I=0
DDstarₛᵢ_A₁⁺0_T₁⁺_I0_nonlocal(i) = DDstarₛᵢ_I0_nonlocal(i, [0, 0, 0], [0, 0, 0])

# Angular momentum: A₁⁺, Ptot²=1; spin: T₁⁺; I=0
DDstarₛᵢ_A₁⁺1_T₁⁺_I0_nonlocal(i) = begin
    O = []

    for (p_x, p_y, p_z) in Iterators.product(0:1, 0:1, 0:1)
        p = [p_x, p_y, p_z]
        if p'*p == 1
            append!(O, DDstarₛᵢ_I0_nonlocal(i, p, -p))
        end
    end
    O = CIP.scalar_mul(O, 1/√6)
    return O
end


# %%#############
# Local Operators
#################

# Generic DD operator
#####################

DDstarₛᵢ_I0_local(i, p) = [
    Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2),
    Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>-1/√2)
]

DDstarₛᵢ_I1_local(i, p) = [
    Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2),
    Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>1/√2)
]


# T₁⁺ operators
###############
# Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=0
DDstarₛᵢ_A₁⁺p0_T₁⁺_I0_local(i) = DDstarₛᵢ_I0_local(i, [0, 0, 0])

# Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=1
DDstarₛᵢ_A₁⁺p0_T₁⁺_I1_local(i) = DDstarₛᵢ_I1_local(i, [0, 0, 0])
