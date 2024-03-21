struct DiscreteMapping_LR{tT,mxT,vT}
    ts::Vector{tT}
    LmappingMX::mxT
    RmappingMX::mxT
    mappingVs::Vector{vT}
end

###############################################################################
############################## Discrete Mapping ###############################
###############################################################################
function DiscreteMapping_LR(LDDEP::AbstractLDDEProblem, method::DiscretizationMethod, DiscretizationLength::Real; n_steps::Int64=nStepOfLength(DiscretizationLength, method.Δt), calculate_additive::Bool=false)
    result = calculateResults(LDDEP, method, DiscretizationLength, n_steps=n_steps, calculate_additive=calculate_additive)
    DiscreteMapping_LR(DiscreteMappingSteps_LR(result)...)
end

function DiscreteMapping_LR_directSparse(LDDEP::AbstractLDDEProblem, method::DiscretizationMethod, DiscretizationLength::Real; n_steps::Int64=nStepOfLength(DiscretizationLength, method.Δt), calculate_additive::Bool=false)
    result = calculateResults(LDDEP, method, DiscretizationLength, n_steps=n_steps, calculate_additive=calculate_additive)
    DiscreteMapping_LR(DiscreteMappingSteps_LR_directSparse(result)...)
end
function rangeshift_LR!(rst::AbstractResult{d}) where {d}
    for (iIPR, smx_IPR) in enumerate(rst.subMXs)#P,R1,R2....
        for (it, smx) in enumerate(smx_IPR) #each time
            for iloc in eachindex(smx.ranges, smx.MXs)
                smx.ranges[iloc] =
                    (smx.ranges[iloc][1] .- (it - 1 - (rst.n_steps - 1)) * d,
                        smx.ranges[iloc][2] .- (it - 1 - (rst.n_steps)) * d)
            end
        end
    end
end


function DiscreteMappingSteps_LR_directSparse(rst::AbstractResult{d}) where {d}
    rangeshift_LR!(rst)

    i = zeros(Int, 0)
    j = zeros(Int, 0)
    v = zeros(typeof(rst.subMXs[1][1].MXs[1][1]), 0)

    for (iIPR, smx_IPR) in enumerate(rst.subMXs)#P,R1,R2....
        for (it, smx) in enumerate(smx_IPR) #each time
            for iloc in eachindex(smx.ranges, smx.MXs)

                #println("----")
                #println([iloc, it, iIPR])
                push!(i, collect(smx.ranges[iloc][1]) .* ones(1, d)...)
                push!(j, ones(d, 1) .* collect(smx.ranges[iloc][2])'...)
                push!(v, -smx.MXs[iloc]...)
                #i[:,:,iloc,it,iIPR] .= getindex.(collect(eachindex(IndexCartesian(),smloc)),1)
                #j[:,:,iloc,it,iIPR] .= getindex.(collect(eachindex(IndexCartesian(),smloc)),2)
                #v[:,:,iloc,it,iIPR] .= smloc
            end
        end
    end

    push!(i, (1:rst.n_steps*d)...)
    push!(j, (1:rst.n_steps*d)...)
    push!(v, 1.0 .* ones(rst.n_steps * d, 1)...)
    #    sparse(i, j, v,[ m, n, combine])

    r = rst.n ÷ d
    p = rst.n_steps
    rhat = maximum([r, p - 1])

    PHI = sparse(i, j, v, p * d, (rhat + p + 1) * d)
    dropzeros!(PHI)

    #PHIL = -PHI[1:(rhat+1)*d, 1:(rhat+1)*d]
    #PHIR = PHI[1:(rhat+1)*d, (rhat+1)*d+1:end]
    PHIL = -PHI[1:(p)*d, 1:(p)*d]
    PHIR = PHI[1:(p)*d, (p)*d+1:end]

    PHILL = vcat(
        hcat(PHIL, spzeros((p) * d, (rhat + 1 - p) * d)),
        hcat(spzeros((rhat + 1 - p) * d, (p) * d), spdiagm(0 => ones((rhat + 1 - p) * d)))
    )
    PHIRR = vcat(
        PHIR,
        hcat(spdiagm(0 => ones((rhat + 1 - p) * d)), spzeros((rhat + 1 - p) * d, (p) * d))
    )

    mappingVs = [spzeros(size(PHILL, 1))]
    for (it, subV) in enumerate(rst.subVs)
        # pos=(it-1)*d
        # mappingVs[1][(1+pos) : (d+pos)] .+= subV.V
        mappingVs[1][(1+(it-1)*d):(d+(it-1)*d)] .+= subV.V
    end

    ([rst.ts[1], rst.ts[end]], PHILL, PHIRR, mappingVs)
end

function DiscreteMappingSteps_LR(rst::AbstractResult{d}) where {d}
    NT = size(rst.subMXs[1])[1]
    Nτ = rst.n ÷ d
    Nmaximal = max(Nτ, NT)

    ΦL = spdiagm(0 => -ones(d * Nmaximal))
    if Nτ > NT
        ΦR = spdiagm(-d * (NT) => -ones(d * (Nτ - NT)))
    else
        ΦR = spzeros(d * Nmaximal, d * Nmaximal)
    end

    for smx_t in rst.subMXs
        for (it, smx) in enumerate(smx_t)
            # pos=(it-1)*d
            # addSubmatrixToL!(ΦL,smx,pos,pos,NT*d)
            # addSubmatrixToR!(ΦR,smx,pos,pos-NT*d,Nmaximal*d)
            addSubmatrixToL!(ΦL, smx, (it - 1) * d, it * d, NT * d)
            addSubmatrixToR!(ΦR, smx, (it - 1) * d, (it - NT) * d, Nmaximal * d)
        end
    end

    mappingVs = [spzeros(size(ΦL, 1))]
    for (it, subV) in enumerate(rst.subVs)
        # pos=(it-1)*d
        # mappingVs[1][(1+pos) : (d+pos)] .+= subV.V
        mappingVs[1][(1+(it-1)*d):(d+(it-1)*d)] .+= subV.V
    end
    #mappingVs=vcat([subV.V for (it, subV) in enumerate(rst.subVs)]...)
    ([rst.ts[1], rst.ts[end]], ΦL, ΦR, mappingVs)
end

function addSubmatrixToL!(Φ::AbstractArray, subMX::SubMX, it_shift1, it_shift2, nlimit)
    for i in eachindex(subMX.ranges, subMX.MXs)
        positionrange = LR_positionshift(subMX.ranges[i], it_shift1, it_shift2)
        if (positionrange[2][end]) <= nlimit && (positionrange[2][1]) > 0
            Φ[positionrange...] .+= subMX.MXs[i]
        end
    end
end

function addSubmatrixToR!(Φ::AbstractArray, subMX::SubMX, it_shift1, it_shift2, nlimit)
    for i in eachindex(subMX.ranges, subMX.MXs)
        positionrange = LR_positionshift(subMX.ranges[i], it_shift1, it_shift2)
        if (positionrange[2][1]) > 0
            Φ[positionrange...] .-= subMX.MXs[i]
        end
    end
end


function LR_positionshift(range::Tuple{UnitRange{Int64},UnitRange{Int64}}, d)
    ((range[1][1]+d):(range[1][end]+d), (range[2][1]+d):(range[2][end]+d))
end

function LR_positionshift(range::Tuple{UnitRange{Int64},UnitRange{Int64}}, d1, d2)
    ((range[1][1]+d1):(range[1][end]+d1), (range[2][1]+d2):(range[2][end]+d2))
end


function spectralRadiusOfMapping(mappLR::DiscreteMapping_LR; args...)
    return abs(eigs(mappLR.RmappingMX, mappLR.LmappingMX; args...)[1][1])
end

function fixPointOfMapping(mappLR::DiscreteMapping_LR)
    (mappLR.LmappingMX - mappLR.RmappingMX) \ Vector(mappLR.mappingVs[1])
end
