struct DiscreteMapping_LR{tT,mxT,vT,sysT}
    ts::Vector{tT}
    LmappingMX::mxT
    RmappingMX::mxT
    mappingVs::Vector{vT}
    A_fixpoint::sysT # Pre-assembled (L - R) matrix for fast fixed point solve
end

###############################################################################
############################## Discrete Mapping ###############################
###############################################################################

#Direct Sparse fill (in one step)
function DiscreteMapping_LR(LDDEP::AbstractLDDEProblem, method::DiscretizationMethod, DiscretizationLength::Real; n_steps::Int64=nStepOfLength(DiscretizationLength, method.Δt), calculate_additive::Bool=false)
    result = calculateResults(LDDEP, method, DiscretizationLength, n_steps=n_steps, calculate_additive=calculate_additive)
    DiscreteMapping_LR(DiscreteMappingSteps_LR(result)...)
end

function rangeshift_LR!(rst::AbstractResult{d}) where {d}
    #for (iIPR, smx_IPR) in enumerate(rst.subMXs)#P,R1,R2....
    for smx_IPR in rst.subMXs#P,R1,R2....
        for (it, smx) in enumerate(smx_IPR) #each time
            for iloc in eachindex(smx.ranges, smx.MXs)
                smx.ranges[iloc] =
                    (smx.ranges[iloc][1] .- (it - 1 - (rst.n_steps - 1)) * d,
                        smx.ranges[iloc][2] .- (it - 1 - (rst.n_steps)) * d)
            end
        end
    end
end


function DiscreteMappingSteps_LR(rst::AbstractResult{d}) where {d}
    rangeshift_LR!(rst)
    
    p = rst.n_steps
    r = rst.n ÷ d
    rhat = maximum([r, p - 1])
    n_full = (rhat + 1) * d
    
    # Pre-calculate counts for L and R parts
    nelements_L_core = 0
    nelements_R_core = 0
    for smx_IPR in rst.subMXs
        for smx in smx_IPR
            for range_pair in smx.ranges
                cols = range_pair[2]
                if cols[end] <= p * d
                    nelements_L_core += d * d
                else
                    for col_idx in cols
                        if col_idx <= p * d
                            nelements_L_core += d
                        else
                            nelements_R_core += d
                        end
                    end
                end
            end
        end
    end

    n_extra = (rhat + 1 - p) * d
    
    # Coordinates for L
    iL = Vector{Int}(undef, nelements_L_core + p * d + n_extra)
    jL = Vector{Int}(undef, nelements_L_core + p * d + n_extra)
    vL = Vector{typeof(rst.subMXs[1][1].MXs[1][1])}(undef, nelements_L_core + p * d + n_extra)
    
    # Coordinates for R
    iR = Vector{Int}(undef, nelements_R_core + n_extra)
    jR = Vector{Int}(undef, nelements_R_core + n_extra)
    vR = Vector{typeof(rst.subMXs[1][1].MXs[1][1])}(undef, nelements_R_core + n_extra)
    
    # Coordinates for A = L - R
    iA = Vector{Int}(undef, length(iL) + length(iR))
    jA = Vector{Int}(undef, length(jL) + length(jR))
    vA = Vector{typeof(rst.subMXs[1][1].MXs[1][1])}(undef, length(vL) + length(vR))

    currL = 1
    currR = 1
    currA = 1
    for smx_IPR in rst.subMXs
        for smx in smx_IPR
            for (range_pair, mx) in zip(smx.ranges, smx.MXs)
                rows = range_pair[1]
                cols = range_pair[2]
                for (c_idx, col_idx) in enumerate(cols)
                    if col_idx <= p * d
                        # Unroll for d=2 if possible
                        if d == 2
                            # row 1
                            val1 = -mx[1, c_idx]
                            row1 = rows[1]
                            iL[currL] = row1; jL[currL] = col_idx; vL[currL] = val1; currL += 1
                            iA[currA] = row1; jA[currA] = col_idx; vA[currA] = val1; currA += 1
                            # row 2
                            val2 = -mx[2, c_idx]
                            row2 = rows[2]
                            iL[currL] = row2; jL[currL] = col_idx; vL[currL] = val2; currL += 1
                            iA[currA] = row2; jA[currA] = col_idx; vA[currA] = val2; currA += 1
                        else
                            for r_idx in 1:d
                                val = -mx[r_idx, c_idx]
                                row = rows[r_idx]
                                iL[currL] = row; jL[currL] = col_idx; vL[currL] = val; currL += 1
                                iA[currA] = row; jA[currA] = col_idx; vA[currA] = val; currA += 1
                            end
                        end
                    else
                        col = col_idx - p * d
                        if d == 2
                            # row 1
                            val1 = mx[1, c_idx]
                            row1 = rows[1]
                            iR[currR] = row1; jR[currR] = col; vR[currR] = val1; currR += 1
                            iA[currA] = row1; jA[currA] = col; vA[currA] = -val1; currA += 1
                            # row 2
                            val2 = mx[2, c_idx]
                            row2 = rows[2]
                            iR[currR] = row2; jR[currR] = col; vR[currR] = val2; currR += 1
                            iA[currA] = row2; jA[currA] = col; vA[currA] = -val2; currA += 1
                        else
                            for r_idx in 1:d
                                val = mx[r_idx, c_idx]
                                row = rows[r_idx]
                                iR[currR] = row; jR[currR] = col; vR[currR] = val; currR += 1
                                iA[currA] = row; jA[currA] = col; vA[currA] = -val; currA += 1
                            end
                        end
                    end
                end
            end
        end
    end

    # PHIL diagonal identity
    for k in 1:(p * d)
        iL[currL] = k; jL[currL] = k; vL[currL] = 1.0; currL += 1
        iA[currA] = k; jA[currA] = k; vA[currA] = 1.0; currA += 1
    end
    
    # Extra parts for PHILL and PHIRR
    for k in 1:n_extra
        idx = p * d + k
        iL[currL] = idx; jL[currL] = idx; vL[currL] = 1.0; currL += 1
        iA[currA] = idx; jA[currA] = idx; vA[currA] = 1.0; currA += 1
    end
    
    for k in 1:n_extra
        row = p * d + k
        col = k
        iR[currR] = row; jR[currR] = col; vR[currR] = 1.0; currR += 1
        iA[currA] = row; jA[currA] = col; vA[currA] = -1.0; currA += 1
    end

    PHILL = sparse(iL, jL, vL, n_full, n_full)
    PHIRR = sparse(iR, jR, vR, n_full, n_full)
    A_fix = sparse(iA, jA, vA, n_full, n_full)
    
    mappingV = zeros(eltype(vL), n_full)
    for (it, subV) in enumerate(rst.subVs)
        target_block = p - it + 1
        pos = (target_block - 1) * d
        mappingV[(1+pos):(d+pos)] .+= subV.V
    end
    mappingVs = [mappingV]
    
    ([rst.ts[1], rst.ts[end]], PHILL, PHIRR, mappingVs, A_fix)
end

function spectralRadiusOfMapping(mappLR::DiscreteMapping_LR{tT,mxT,vT}; useKrylovKit=true, nev=1, tol=1e-6, args...)::mxT.parameters[1] where {tT,mxT,vT}
    if useKrylovKit
        # Use KrylovKit on the operator L \ (R * x) to find the largest magnitude eigenvalue
        # This works for non-symmetric matrices and avoids geneigsolve's symmetric requirement
        vals, vecs, info = eigsolve(x -> mappLR.LmappingMX \ (mappLR.RmappingMX * x), rand(eltype(mappLR.RmappingMX), size(mappLR.RmappingMX, 1)), nev, :LM; tol=tol, args...)
        return abs(vals[1])::mxT.parameters[1]
    else
        # Use Arpack for generalized eigenvalue problem: R*x = lambda*L*x
        return abs(Arpack.eigs(mappLR.RmappingMX, mappLR.LmappingMX; nev=nev, tol=tol, args...)[1][1])::mxT.parameters[1]
    end
end

function fixPointOfMapping(mappLR::DiscreteMapping_LR{tT,mxT,vT}) where {tT,mxT,vT}
    # Use the pre-assembled system matrix for fix point solve
    return mappLR.A_fixpoint \ mappLR.mappingVs[1]
end
