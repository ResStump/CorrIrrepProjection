# %%########
# Tetraquark
############

function get_DD_nonlocal_corr(O_sink, O_src, raw_corr_dict)
    Γ_idx_permutation = [1, 2, 3, 4]

    raw_correlators = raw_corr_dict["DD_nonlocal"]["Correlators"]
    Γ_labels = raw_corr_dict["DD_nonlocal"]["Spin Structure"]

    # Total momentum
    p₁, p₂, p₃, p₄ = O_sink["p"]..., O_src["p"]...
    @assert p₁ + p₂ == p₃ + p₄ "Total momenta don't match."
    Ptot_str = "Ptot"*join(p₁ + p₂, ",")

    # Check if p₁ or p₂ is psink1 in raw correlator
    swap_ud = false
    psink1_str_arr = keys(raw_correlators[Ptot_str])
    if (p_str = "psink1_"*join(p₁, ",")) in psink1_str_arr
        psink1_str = p_str
    elseif (p_str = "psink1_"*join(p₂, ",")) in psink1_str_arr
        psink1_str = p_str
        
        # In that case also swap Γ₁<->Γ₂ (and possibly u<->d at source)
        permute!(Γ_idx_permutation, [2, 1, 3, 4])
        swap_ud = !swap_ud
    else
        throw(ArgumentError("none of the momenta p₁ and p₂ is in the correlator."))
    end

    # Check if p₃ or p₄ is psrc1 in raw correlator
    psrc1_str_arr = keys(raw_correlators[Ptot_str][psink1_str])
    if (p_str = "psrc1_"*join(p₃, ",")) in psrc1_str_arr
        psrc1_str = p_str
    elseif  (p_str = "psrc1_"*join(p₄, ",")) in psrc1_str_arr
        psrc1_str = p_str
        
        # In that case also swap Γ₃<->Γ₄ (and possibly u<->d at source)
        permute!(Γ_idx_permutation, [1, 2, 4, 3])
        swap_ud = !swap_ud
    else
        throw(ArgumentError("none of the momenta p₃ and p₄ is in the correlator."))
    end

    # Flavour content
    @assert O_sink["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of sink operator not valid."
    @assert O_src["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of source operator not valid."

    flavour = "ubar_c_dbar_c"
    if (O_sink["flavour"] == O_src["flavour"]) ⊻ swap_ud
        flavour *= "-cbar_u_cbar_d"
    else
        flavour *= "-cbar_d_cbar_u"
    end

    # Indices for gamma matrices
    Γ_idx = [
        Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_DD_1"]),
        Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_DD_2"]),
        Γ_str_to_idx(O_src["Γ₁"], Γ_labels["Gamma_DD_1"]),
        Γ_str_to_idx(O_src["Γ₂"], Γ_labels["Gamma_DD_2"])
    ]
    permute!(Γ_idx, Γ_idx_permutation)

    # Extract relevant part from raw correlators
    return raw_correlators[Ptot_str][psink1_str][psrc1_str][flavour][:, Γ_idx...]
end

function get_DD_local_corr(O_sink, O_src, raw_corr_dict)
    Γ_idx_permutation = [1, 2, 3, 4]

    raw_correlators = raw_corr_dict["DD_local"]["Correlators"]
    Γ_labels = raw_corr_dict["DD_local"]["Spin Structure"]

    # Momentum
    p′, p = O_sink["p"], O_src["p"]
    @assert p′ == p "Momenta don't match."
    p_str = "p"*join(p, ",")

    # Flavour
    @assert O_sink["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of sink operator not valid."
    @assert O_src["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of source operator not valid."
    if (O_sink["flavour"] != O_src["flavour"])
        # Swap Γ₃<->Γ₄
        permute!(Γ_idx_permutation, [1, 2, 4, 3])
    end

    # Indices for gamma matrices
    Γ_idx = [
        Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_DD_1"]),
        Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_DD_2"]),
        Γ_str_to_idx(O_src["Γ₁"], Γ_labels["Gamma_DD_1"]),
        Γ_str_to_idx(O_src["Γ₂"], Γ_labels["Gamma_DD_2"])
    ]
    permute!(Γ_idx, Γ_idx_permutation)
    
    # Extract relevant part from raw correlators
    return raw_correlators[p_str][:, Γ_idx...]
end

function get_DD_mixed_corr(O_sink, O_src, raw_corr_dict)
    Γ_idx_permutation = [1, 2, 3, 4]

    raw_correlators = raw_corr_dict["DD_mixed"]["Correlators"]
    Γ_labels = raw_corr_dict["DD_mixed"]["Spin Structure"]

    # Check if flavour content is valid
    @assert O_sink["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of sink operator not valid."
    @assert O_src["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of source operator not valid."

    # See which operator is at sink and which at source
    if O_sink["type"] == :DD_nonlocal && O_src["type"] == :DD_local
        # Order of operators
        op_order = "nonlocal-local"

        # Get momenta
        p₁, p₂ = O_sink["p"]
        Ptot = O_src["p"]

        # Flavour
        if (O_sink["flavour"] != O_src["flavour"])
            # Swap Γ₃<->Γ₄ (swap u<->d in the local operator)
            permute!(Γ_idx_permutation, [1, 2, 4, 3])
        end
    else
        # Order of operators
        op_order = "local-nonlocal"

        # Get momenta
        Ptot = O_sink["p"]
        p₁, p₂ = O_src["p"]

        # Flavour
        if (O_sink["flavour"] != O_src["flavour"])
            # Swap Γ₁<->Γ₂ (swap u<->d in the local operator)
            permute!(Γ_idx_permutation, [2, 1, 3, 4])
        end
    end

    # Total momentum
    @assert Ptot == p₁ + p₂ "Total momenta don't match."
    Ptot_str = "Ptot"*join(Ptot, ",")

    # Check if p₁ or p₂ is p_nonlocal1 in raw correlator
    p_nonlocal_str_arr = keys(raw_correlators[Ptot_str])
    if (p_str = "p_nonlocal1_"*join(p₁, ",")) in p_nonlocal_str_arr
        p_nonlocal1_str = p_str
    elseif (p_str = "p_nonlocal1_"*join(p₂, ",")) in p_nonlocal_str_arr
        p_nonlocal1_str = p_str
        
        # In that case also swap Γ₁<->Γ₂ and Γ₃<->Γ₄
        permute!(Γ_idx_permutation, [2, 1, 4, 3])
    else
        throw(ArgumentError("none of the chose momenta for the nonlocal interpolator " *
                            "is in the correlator."))
    end

    # Indices for gamma matrices
    Γ_idx = [
        Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_DD_1"]),
        Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_DD_2"]),
        Γ_str_to_idx(O_src["Γ₁"], Γ_labels["Gamma_DD_1"]),
        Γ_str_to_idx(O_src["Γ₂"], Γ_labels["Gamma_DD_2"])
    ]
    permute!(Γ_idx, Γ_idx_permutation)
    
    # Extract relevant part from raw correlators
    return raw_correlators[Ptot_str][p_nonlocal1_str][op_order][:, Γ_idx...]
end

function get_dad_local_corr(O_sink, O_src, raw_corr_dict)
    sign = 1

    raw_correlators = raw_corr_dict["dad_local"]["Correlators"]
    Γ_labels = raw_corr_dict["dad_local"]["Spin Structure"]

    # Momentum
    p′, p = O_sink["p"], O_src["p"]
    @assert p′ == p "Momenta don't match."
    p_str = "p"*join(p, ",")

    # Flavour
    @assert O_sink["flavour"] in [:ccūd̄, :ccd̄ū] "Flavour of sink operator not valid."
    @assert O_src["flavour"] in [:ccūd̄, :ccd̄ū] "Flavour of source operator not valid."
    if (O_sink["flavour"] == :ccd̄ū) ⊻ (O_src["flavour"] == :ccd̄ū)
        sign *= -1
    end

    # Indices for gamma matrices
    @assert O_sink["Γ₁"] == O_src["Γ₁"] "Diquark-antidiquark correlators with different " *
                                        "spin structure at sink and source not implemented."
    @assert O_sink["Γ₂"] == O_src["Γ₂"] "Diquark-antidiquark correlators with different " *
                                        "spin structure at sink and source not implemented."
    Γ_idx = [
        Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_dad_1"]),
        Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_dad_2"])
    ]
    
    # Extract relevant part from raw correlators (and multiply with sign)
    return sign * raw_correlators[p_str][:, Γ_idx...]
end

function get_dad_DD_local_mixed_corr(O_sink, O_src, raw_corr_dict)
    Γ_idx_permutation = [1, 2, 3, 4]
    sign = 1

    raw_correlators = raw_corr_dict["dad-DD_local"]["Correlators"]
    Γ_labels = raw_corr_dict["dad-DD_local"]["Spin Structure"]

    # Momentum
    p′, p = O_sink["p"], O_src["p"]
    @assert p′ == p "Momenta don't match."
    p_str = "p"*join(p, ",")

    # See which operator is at sink and which at source
    if O_sink["type"] == :dad_local && O_src["type"] == :DD_local
        # Order of operators
        op_order = "dad-DD"

        # Flavour content
        @assert O_sink["flavour"] in [:ccūd̄, :ccd̄ū] "Flavour of sink operator not valid."
        @assert O_src["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of source operator not valid."
        if O_sink["flavour"] == :ccd̄ū
            sign *= -1
        end
        if O_src["flavour"] == :d̄cūc
            # Swap Γ₃<->Γ₄
            permute!(Γ_idx_permutation, [1, 2, 4, 3])
        end

        # Indices for gamma matrices
        Γ_idx = [
            Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_dad_1"]),
            Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_dad_2"]),
            Γ_str_to_idx(O_src["Γ₁"], Γ_labels["Gamma_DD_1"]),
            Γ_str_to_idx(O_src["Γ₂"], Γ_labels["Gamma_DD_2"])
        ]
    else
        # Order of operators
        op_order = "DD-dad"

        # Flavour content
        @assert O_sink["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of sink operator not valid."
        @assert O_src["flavour"] in [:ccūd̄, :ccd̄ū] "Flavour of source operator not valid."
        if O_src["flavour"] == :ccd̄ū
            sign *= -1
        end
        if O_sink["flavour"] == :d̄cūc
            # Swap Γ₁<->Γ₂
            permute!(Γ_idx_permutation, [2, 1, 3, 4])
        end

        # Indices for gamma matrices
        Γ_idx = [
            Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_DD_1"]),
            Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_DD_2"]),
            Γ_str_to_idx(O_src["Γ₁"], Γ_labels["Gamma_dad_1"]),
            Γ_str_to_idx(O_src["Γ₂"], Γ_labels["Gamma_dad_2"])
        ]
    end
    
    permute!(Γ_idx, Γ_idx_permutation)

    # Extract relevant part from raw correlators (and multiply with sign)
    return sign * raw_correlators[p_str][op_order][:, Γ_idx...]
end

function get_dad_DD_local_nonlocal_mixed_corr(O_sink, O_src, raw_corr_dict)
    Γ_idx_permutation = [1, 2, 3, 4]
    sign = 1

    raw_correlators = raw_corr_dict["dad-DD_local-nonlocal"]["Correlators"]
    Γ_labels = raw_corr_dict["dad-DD_local-nonlocal"]["Spin Structure"]

    # See which operator is at sink and which at source
    if O_sink["type"] == :dad_local && O_src["type"] == :DD_nonlocal
        # Order of operators
        op_order = "local-nonlocal"

        # Get momenta
        Ptot = O_sink["p"]
        p₃, p₄ = O_src["p"]
        @assert Ptot == p₃ + p₄ "Total momenta don't match."
        Ptot_str = "Ptot"*join(Ptot, ",")

        # Check if p₃ or p₄ is p_nonlocal1 in raw correlator
        p_nonlocal_str_arr = keys(raw_correlators[Ptot_str])
        if (p_str = "p_nonlocal1_"*join(p₃, ",")) in p_nonlocal_str_arr
            p_nonlocal1_str = p_str
        elseif (p_str = "p_nonlocal1_"*join(p₄, ",")) in p_nonlocal_str_arr
            p_nonlocal1_str = p_str
            
            # Swap Γ₃<->Γ₄ and flip sign (since dad operator is I=0)
            permute!(Γ_idx_permutation, [1, 2, 4, 3])
            sign *= -1
        else
            throw(ArgumentError("none of the chose momenta for the nonlocal interpolator " *
                                "is in the correlator."))
        end

        # Flavour content
        @assert O_sink["flavour"] in [:ccūd̄, :ccd̄ū] "Flavour of sink operator not valid."
        @assert O_src["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of source operator not valid."
        if (O_sink["flavour"] == :ccd̄ū) ⊻ (O_src["flavour"] == :d̄cūc)
            # Swap u<->d in dad operator which flips sign since it is I=0
            sign *= -1
        end

        # Indices for gamma matrices
        Γ_idx = [
            Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_dad_1"]),
            Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_dad_2"]),
            Γ_str_to_idx(O_src["Γ₁"], Γ_labels["Gamma_DD_1"]),
            Γ_str_to_idx(O_src["Γ₂"], Γ_labels["Gamma_DD_2"])
        ]
    else
        # Order of operators
        op_order = "nonlocal-local"

        # Get momenta
        p₁, p₂ = O_sink["p"]
        Ptot = O_src["p"]
        @assert Ptot == p₁ + p₂ "Total momenta don't match."
        Ptot_str = "Ptot"*join(Ptot, ",")

        # Check if p₁ or p₂ is p_nonlocal1 in raw correlator
        p_nonlocal_str_arr = keys(raw_correlators[Ptot_str])
        if (p_str = "p_nonlocal1_"*join(p₁, ",")) in p_nonlocal_str_arr
            p_nonlocal1_str = p_str
        elseif (p_str = "p_nonlocal1_"*join(p₂, ",")) in p_nonlocal_str_arr
            p_nonlocal1_str = p_str
            
            # Swap Γ₁<->Γ₂ and flip sign (since dad operator is I=0)
            permute!(Γ_idx_permutation, [2, 1, 3, 4])
            sign *= -1
        else
            throw(ArgumentError("none of the chose momenta for the nonlocal interpolator " *
                                "is in the correlator."))
        end

        # Flavour
        @assert O_sink["flavour"] in [:ūcd̄c, :d̄cūc] "Flavour of sink operator not valid."
        @assert O_src["flavour"] in [:ccūd̄, :ccd̄ū] "Flavour of source operator not valid."
        if (O_sink["flavour"] == :d̄cūc) ⊻ (O_src["flavour"] == :ccd̄ū)
            # Swap u<->d in dad operator which flips sign since it is I=0
            sign *= -1
        end

        # Indices for gamma matrices
        Γ_idx = [
            Γ_str_to_idx(O_sink["Γ₁"], Γ_labels["Gamma_DD_1"]),
            Γ_str_to_idx(O_sink["Γ₂"], Γ_labels["Gamma_DD_2"]),
            Γ_str_to_idx(O_src["Γ₁"], Γ_labels["Gamma_dad_1"]),
            Γ_str_to_idx(O_src["Γ₂"], Γ_labels["Gamma_dad_2"])
        ]
    end
    
    permute!(Γ_idx, Γ_idx_permutation)
    
    # Extract relevant part from raw correlators (and multiply with sign)
    return sign * raw_correlators[Ptot_str][p_nonlocal1_str][op_order][:, Γ_idx...]
end

@doc raw"""
    get_tetraquark_corr(O_sink, O_src, raw_corr_dict) -> Cₜ

Get the correlator entries corresponding to Cₜ = <`O_sink`(t) `O_src`(0)^†> from the Dict
of raw correlators `raw_corr_dict`. This Dict is expected to contain the raw correlators
(key "Correlators") and the spin structure of it (key "Spin Structure).
"""
function get_tetraquark_corr(O_sink, O_src, raw_corr_dict)
    # Both nonlocal DD
    if O_sink["type"] == O_src["type"] == :DD_nonlocal
        Cₜ = get_DD_nonlocal_corr(O_sink, O_src, raw_corr_dict)

    # Both local DD
    elseif O_sink["type"] == O_src["type"] == :DD_local
        Cₜ = get_DD_local_corr(O_sink, O_src, raw_corr_dict)
    
    # Mixed nonlocal-local or local-nonlocal DD
    elseif (O_sink["type"] == :DD_nonlocal && O_src["type"] == :DD_local) ||
        (O_sink["type"] == :DD_local && O_src["type"] == :DD_nonlocal)
        Cₜ = get_DD_mixed_corr(O_sink, O_src, raw_corr_dict)

    # Local diquark-antidiquark
    elseif O_sink["type"] == O_src["type"] == :dad_local
        Cₜ = get_dad_local_corr(O_sink, O_src, raw_corr_dict)

    # Local diquark-antidiquark-DD and DD-diquark-antidiquark
    elseif (O_sink["type"] == :dad_local && O_src["type"] == :DD_local) ||
        (O_sink["type"] == :DD_local && O_src["type"] == :dad_local)
        Cₜ = get_dad_DD_local_mixed_corr(O_sink, O_src, raw_corr_dict)
    
    # Local-nonlocal diquark-antidiquark-DD and nonlocal-local DD-diquark-antidiquark
    elseif (O_sink["type"] == :dad_local && O_src["type"] == :DD_nonlocal) ||
        (O_sink["type"] == :DD_nonlocal && O_src["type"] == :dad_local)
        Cₜ = get_dad_DD_local_nonlocal_mixed_corr(O_sink, O_src, raw_corr_dict)
    
    else
        throw(ArgumentError("the combination of operators of type $(O_sink["type"]) and " *
                            "$(O_src["type"]) is not implemented."))
    end

    # Normalization factor
    Cₜ .*= O_sink["N"]*conj(O_src["N"])

    return Cₜ
end

@doc raw"""
    project_tetraquark_corr(O_sink_arr, O_src_arr, raw_corr_dict) -> Cₜ_fb

Compute the correlator `Cₜ_fb` which is the sum of correlators <O\_sink(t) O\_src(0)^†> for
all combinations of operators (O\_sink, O\_src) from the arrays of operators `O_sink_arr`
and `O_src_arr`. The correlator `Cₜ_fb` is in  forward/backward shape which means
it has the shape (Nₜ//2+1, 2) where the first column contains the forward correlator and the
second column the backward correlator.
The Dict `raw_corr_dict` is expected to contain the raw correlators
(key "Correlators") and the spin structure of it (key "Spin Structure).
"""
function project_tetraquark_corr(O_sink_arr, O_src_arr, raw_corr_dict)
    Nₜ = parms.Nₜ

    # Create Correlator (forward and backward shape)
    Cₜ_fb = zeros(ComplexF64, Nₜ÷2+1, 2)

    # Indices for forward and backward correlator
    fwd_indices = 1:Nₜ÷2+1
    bwd_indices = mod1.(Nₜ+1:-1:Nₜ÷2+1, Nₜ)

    # Get time reversal factors of sink/src operators
    trev_sink_arr = map(d -> d["trev"], O_sink_arr)
    trev_src_arr = map(d -> d["trev"], O_src_arr)

    # If time reversal factors are equal compute full correlator at once
    if allequal(trev_sink_arr) && allequal(trev_src_arr)
        Cₜ = zeros(ComplexF64, Nₜ)

        for (O_sink, O_src) in Iterators.product(O_sink_arr, O_src_arr)
            Cₜ .+= get_tetraquark_corr(O_sink, O_src, raw_corr_dict)
        end

        # Forward correlator
        Cₜ_fb[:, 1] = Cₜ[fwd_indices]

        # Backward correlator (multiply with time reversal factors)
        sign = trev_sink_arr[1]*trev_src_arr[1]
        Cₜ_fb[:, 2] = sign * Cₜ[bwd_indices]
    else
        for (O_sink, O_src) in Iterators.product(O_sink_arr, O_src_arr)
            Cₜ = get_tetraquark_corr(O_sink, O_src, raw_corr_dict)

            # Forward correlator
            Cₜ_fb[:, 1] += Cₜ[fwd_indices]

            # Backward correlator (multiply with time reversal factors)
            sign = O_sink["trev"]*O_src["trev"]
            Cₜ_fb[:, 2] += sign * Cₜ[bwd_indices]
        end
    end

    return Cₜ_fb
end