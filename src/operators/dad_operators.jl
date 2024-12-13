γ = CIP.γ

# Generic diquark-antidiquark operator
######################################

dadᵢ_I0_local(i, p) = [
    Dict("type"=>:dad_local, "flavour"=>:ccūd̄, "p"=>p, "Γ₁"=>"C"*γ[i], "Γ₂"=>"C"*γ[5],
         "N"=>1, "trev"=>-1)
]


# T₁⁺ operators
###############
# Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=0
dadᵢ_A₁⁺_T₁⁺_I0_local(i) = dadᵢ_I0_local(i, [0, 0, 0])
