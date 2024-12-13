# Useful definitions
γ = CIP.γ
e = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
δ(i, j) = Int(i==j)

# %%################
# Nonlocal Operators
####################

# Generic DD operator
#####################
# Subscript ₛ and ₐ denotes symmetric/antisymmetric under eschange of flavour and γ matrices
# Subscript ᵢ denotes that it is spin 1 (in direction i)

DDstarₐ_I0_nonlocal(p₁, p₂) = begin
    @assert p₁ != p₂ "operator is zero for p₁ == p₂"
    O = [
        Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
             "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
             "Γ₁"=>γ[5], "Γ₂"=>γ[5],"N"=>-1/√2, "trev"=>+1)
        ]
    return O
end

DDstarₛ_I1_nonlocal(p₁, p₂) = begin
    if p₁ == p₂
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1)
        ]
    else
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1)
        ]
    end
    return O
end

DDstarₛᵢ_I0_nonlocal(i, p₁, p₂) = begin
    if p₁ == p₂
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/√2, "trev"=>-1)
        ]
    else
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>-1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>1/2, "trev"=>-1)
        ]
    end
    return O
end

DstarDstarₛᵢ_S1_I0_nonlocal(i, p₁, p₂) = begin
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)

    if p₁ == p₂
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>1/√2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>-1/√2, "trev"=>+1)
        ]
    else
        O = [
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>1/2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>-1/2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip2], "Γ₂"=>γ[ip1], "N"=>-1/2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip2], "Γ₂"=>γ[ip1], "N"=>1/2, "trev"=>+1)
        ]
    end
    return O
end


# T₁⁺ operators (rest frame)
############################

# DD*; Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=0
DDstarₛᵢ_A₁⁺0_T₁⁺_I0_nonlocal(i) = DDstarₛᵢ_I0_nonlocal(i, [0, 0, 0], [0, 0, 0])

# Angular momentum: A₁⁺, p²=1; spin: T₁⁺; I=0
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

# DD*; Angular momentum: A₁⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_A₁⁺2_T₁⁺_I0_nonlocal(i) = begin
    O = []

    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            append!(O, DDstarₛᵢ_I0_nonlocal(i, p, -p))
            push!(used_p, p)
        end
    end
    O = CIP.scalar_mul(O, 1/√12)
    return O
end

# DD*; Angular momentum: E⁺, p²=1; spin: T₁⁺; I=0
DDstarₛᵢ_E⁺1_T₁⁺_I0_nonlocal(i) = begin
    O = []

    append!(O, DDstarₛᵢ_I0_nonlocal(i, e[i], -e[i]))
    for j in 1:3
        O_ = DDstarₛᵢ_I0_nonlocal(i, e[j], -e[j])
        append!(O, CIP.scalar_mul(O_, -1/3))
    end
    O = CIP.scalar_mul(O, √3/2)
    return O
end

# DD*; Angular momentum: J^P = 1⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_J1⁺2_T₁⁺_I0_nonlocal(i) = begin
    O = []

    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            for j in 1:3
                factor = p[i]*p[j] - δ(i, j)/3 * p'*p
                if factor != 0
                    append!(O, CIP.scalar_mul(DDstarₛᵢ_I0_nonlocal(j, p, -p), factor))
                end    
            end
            push!(used_p, p)
        end
    end
    
    O = CIP.scalar_mul(O, 1/√2)
    return O
end

# DD*; Angular momentum: J^P = 3⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_J3⁺2_T₁⁺_I0_nonlocal(i) = begin
    O = []

    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            factor = 3/5*(p[i]^2 - 1/3 * p'*p)
            if factor != 0
                append!(O, CIP.scalar_mul(DDstarₛᵢ_I0_nonlocal(i, p, -p), factor))
            end

            for j in 1:3
                j != i || continue

                factor = -2/5*p[i]*p[j]
                if factor != 0
                    append!(O, CIP.scalar_mul(DDstarₛᵢ_I0_nonlocal(j, p, -p), factor))
                end 
            end
            push!(used_p, p)
        end
    end
    
    O = CIP.scalar_mul(O, 1/√2)
    return O
end

# D*D*; Angular momentum: A₁⁺, p²=0; spin: S=1, T₁⁺; I=0; 
DstarDstarₛᵢ_A₁⁺0_S1_T₁⁺_I0_nonlocal(i) =
    DstarDstarₛᵢ_S1_I0_nonlocal(i, [0, 0, 0], [0, 0, 0])

# D*D*; Angular momentum: A₁⁺, p²=1; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_A₁⁺1_S1_T₁⁺_I0_nonlocal(i) = begin
    O = []

    for (p_x, p_y, p_z) in Iterators.product(0:1, 0:1, 0:1)
        p = [p_x, p_y, p_z]
        if p'*p == 1
            append!(O, DstarDstarₛᵢ_S1_I0_nonlocal(i, p, -p))
        end
    end
    O = CIP.scalar_mul(O, 1/√6)
    return O
end

# D*D*; Angular momentum: E⁺, p²=1; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_E⁺1_S1_T₁⁺_I0_nonlocal(i) = begin
    O = []

    append!(O, DstarDstarₛᵢ_S1_I0_nonlocal(i, e[i], -e[i]))
    for j in 1:3
        O_ = DstarDstarₛᵢ_S1_I0_nonlocal(i, e[j], -e[j])
        append!(O, CIP.scalar_mul(O_, -1/3))
    end
    O = CIP.scalar_mul(O, √3/2)
    return O
end


# %%#############
# Local Operators
#################

# Generic DD operator
#####################

DDstarₛᵢ_I0_local(i, p) = [
    Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2,
         "trev"=>-1),
    Dict("type"=>:DD_local, "flavour"=>:ūcd̄c, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/√2,
         "trev"=>-1)
]

DDstarₛᵢ_I1_local(i, p) = [
    Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2,
         "trev"=>-1),
    Dict("type"=>:DD_local, "flavour"=>:ūcd̄c, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2,
         "trev"=>-1)
]

DstarDstarₛᵢ_S1_I0_local(i, p) = begin
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)

    [
        Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2],
             "N"=>1/√2, "trev"=>1),
        Dict("type"=>:DD_local, "flavour"=>:ūcd̄c, "p"=>p, "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2],
             "N"=>-1/√2, "trev"=>1)
    ]
end


# T₁⁺ operators
###############
# DD*; Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=0
DDstarₛᵢ_A₁⁺_T₁⁺_I0_local(i) = DDstarₛᵢ_I0_local(i, [0, 0, 0])

# DD*; Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=1
DDstarₛᵢ_A₁⁺_T₁⁺_I1_local(i) = DDstarₛᵢ_I1_local(i, [0, 0, 0])

# D*D*; Angular momentum: A₁⁺, p²=0; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_A₁⁺_S1_T₁⁺_I0_local(i) = DstarDstarₛᵢ_S1_I0_local(i, [0, 0, 0])
