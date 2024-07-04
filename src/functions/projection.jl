# %%############
# DD, Tetraquark
################

@doc raw"""
    get_tetraquark_corr(O_sink, O_src, raw_corr_dict, Γ_indices) -> Cₜ

Get the correlator entries corresponding to Cₜ = <`O_sink`(t) `O_src`(0)^†> from the Dict
of raw correlators `raw_corr_dict`. The Dict `Γ_indices` contains the indices for the
(monomials) of gamma matrices used to compute the raw correlators.
"""
function get_tetraquark_corr(O_sink, O_src, raw_corr_dict, Γ_indices)
    Γ_idx_permutation = [1, 2, 3, 4]

    # Case: Both nonlocal DD
    ########################
    if O_sink["type"] == O_src["type"] == :DD_nonlocal
        # Total momentum
        p₁, p₂, p₃, p₄ = O_sink["p"]..., O_src["p"]...
        @assert p₁ + p₂ == p₃ + p₄ "Total momenta don't match."
        Ptot_str = "Ptot"*join(p₁ + p₂, ",")

        # Check if p₁ or p₂ is psink1 in raw correlator
        swap_ud = false
        psink1_str_arr = keys(raw_corr_dict["DD_nonlocal"][Ptot_str])
        if "psink1_"*join(p₁, ",") in psink1_str_arr
            psink1_str = "psink1_"*join(p₁, ",")
        else
            psink1_str = "psink1_"*join(p₂, ",")
            
            # In that case also swap Γ₁<->Γ₂ (and possibly u<->d at source)
            permute!(Γ_idx_permutation, [2, 1, 3, 4])
            swap_ud = !swap_ud
        end

        # Check if p₃ or p₄ is psrc1 in raw correlator
        psrc1_str_arr = keys(raw_corr_dict["DD_nonlocal"][Ptot_str][psink1_str])
        if "psrc1_"*join(p₃, ",") in psrc1_str_arr
            psrc1_str = "psrc1_"*join(p₃, ",")
        elseif  "psrc1_"*join(p₄, ",") in psrc1_str_arr
            psrc1_str = "psrc1_"*join(p₄, ",")
            
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

        # Get raw correlator
        raw_corr = raw_corr_dict["DD_nonlocal"][Ptot_str][psink1_str][psrc1_str][flavour]

    # Case: Both local DD
    #####################
    elseif O_sink["type"] == O_src["type"] == :DD_local
        # Momentum
        p′, p = O_sink["p"], O_src["p"]
        @assert p′ == p
        p_str = "p"*join(p, ",")

        # Flavour
        if (O_sink["flavour"] != O_src["flavour"])
            # Swap Γ₃<->Γ₄
            permute!(Γ_idx_permutation, [1, 2, 4, 3])
        end

        # Get war correlator
        raw_corr = raw_corr_dict["DD_local"][p_str]
    
    # Case: mixed nonlocal-local or local-nonlocal DD
    #################################################
    elseif (O_sink["type"] == :DD_nonlocal && O_src["type"] == :DD_local) ||
        (O_sink["type"] == :DD_local && O_src["type"] == :DD_nonlocal)
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
                # Swap Γ₃<->Γ₄ (swap u<->d in the local operator)
                permute!(Γ_idx_permutation, [2, 1, 3, 4])
            end
        end

        # Total momentum
        @assert Ptot == p₁ + p₂ "Total momenta don't match."
        Ptot_str = "Ptot"*join(Ptot, ",")

        # Check if p₁ or p₂ is p_nonlocal1 in raw correlator
        p_nonlocal_str_arr = keys(raw_corr_dict["DD_mixed"][Ptot_str])
        if "p_nonlocal1_"*join(p₁, ",") in p_nonlocal_str_arr
            p_nonlocal1_str = "p_nonlocal1_"*join(p₁, ",")
        elseif "p_nonlocal1_"*join(p₂, ",") in p_nonlocal_str_arr
            p_nonlocal1_str = "p_nonlocal1_"*join(p₂, ",")
            
            # In that case also swap Γ₁<->Γ₂ and Γ₃<->Γ₄
            permute!(Γ_idx_permutation, [2, 1, 4, 3])
        else
            throw(ArgumentError("none of the chose momenta for the nonlocal interpolator " *
                                "is in the correlator."))
        end

        # Get raw correlator
        raw_corr = raw_corr_dict["DD_mixed"][Ptot_str][p_nonlocal1_str][op_order]
    end

    # Indices for gamma matrices
    Γ_idx = [Γ_indices[O_sink["Γ₁"]], Γ_indices[O_sink["Γ₂"]],
             Γ_indices[O_src["Γ₁"]], Γ_indices[O_src["Γ₂"]]]
    permute!(Γ_idx, Γ_idx_permutation)
    
    Cₜ = raw_corr[:, Γ_idx...]

    # Normalization factor
    Cₜ .*= O_sink["N"]*conj(O_src["N"])

    return Cₜ
end

@doc raw"""
    project_tetraquark_corr(O_sink_arr, O_src_arr, raw_corr_dict, Γ_indices) -> Cₜ

Compute the correlator Cₜ which is the sum of correlators <O\_sink(t) O\_src(0)^†> for all
combinations of operators (O\_sink, O\_src) from the arrays of operators `O_sink_arr` and
`O_src_arr`. The raw correlators are stored in the Dict `raw_corr_dict`.  The Dict 
`Γ_indices` contains the indices for the (monomials) of gamma matrices used to compute the
raw correlators.
"""
function project_tetraquark_corr(O_sink_arr, O_src_arr, raw_corr_dict, Γ_indices)
    # Create Correlator
    Cₜ = zeros(ComplexF64, parms.Nₜ)

    for (O_sink, O_src) in Iterators.product(O_sink_arr, O_src_arr)
        Cₜ .+= get_tetraquark_corr(O_sink, O_src, raw_corr_dict, Γ_indices)
    end

    return Cₜ
end