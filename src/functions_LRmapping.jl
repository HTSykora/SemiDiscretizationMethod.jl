struct DiscreteMapping_LR{tT,mxT,vT}
    ts::Vector{tT}
    LmappingMX::mxT
    RmappingMX::mxT
    mappingVs::Vector{vT}
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
    i = zeros(Int, 0)
    j = zeros(Int, 0)
    v = zeros(typeof(rst.subMXs[1][1].MXs[1][1]), 0)


    Nelements_approx = ((methodorder(rst.method) + 1) * (size(rst.subMXs, 1) - 1) + 1) * rst.n_steps * d * d+ rst.n_steps * d 
  # println(Nelements_approx)
    sizehint!(i, Nelements_approx)
    sizehint!(j, Nelements_approx)
    sizehint!(v, Nelements_approx)

    for (iIPR, smx_IPR) in enumerate(rst.subMXs)#P,R1,R2....
        for (it, smx) in enumerate(smx_IPR) #each time
            for iloc in eachindex(smx.ranges, smx.MXs)

                #println("----")
                #println([iloc, it, iIPR])
                push!(i, collect(smx.ranges[iloc][1]) .* ones(Int, 1, d)...)
                push!(j, ones(Int, d, 1) .* transpose(collect(smx.ranges[iloc][2]))...)
                #push!(j, collect(smx.ranges[iloc][2]) .* ones(Int, 1, d)...)
                push!(v, -smx.MXs[iloc]...)
                #i[:,:,iloc,it,iIPR] .= getindex.(collect(eachindex(IndexCartesian(),smloc)),1)
                #j[:,:,iloc,it,iIPR] .= getindex.(collect(eachindex(IndexCartesian(),smloc)),2)
                #v[:,:,iloc,it,iIPR] .= smloc
            end
        end
    end
   # println(size(i))
    push!(i, (1:rst.n_steps*d)...)
    push!(j, (1:rst.n_steps*d)...)
    push!(v, 1.0 .* ones(rst.n_steps * d, 1)...)
    #    sparse(i, j, v,[ m, n, combine])

   # println(size(i)) 
   # println("--------")
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


function spectralRadiusOfMapping(mappLR::DiscreteMapping_LR{tT,mxT,vT}; args...)::mxT.parameters[1] where {tT,mxT,vT}
   # if size(mappLR.LmappingMX,1)>30
        #println("Full")
    #    return abs(eigen(collect(mappLR.RmappingMX),collect(mappLR.LmappingMX),sortby=abs).values[end])
    #else
    #    #println("SP")
    return abs(eigs(mappLR.RmappingMX, mappLR.LmappingMX; args...)[1][1])::mxT.parameters[1]
    #end
end

function fixPointOfMapping(mappLR::DiscreteMapping_LR{tT,mxT,vT}) where {tT,mxT,vT}
    (mappLR.LmappingMX - mappLR.RmappingMX) \ Vector(mappLR.mappingVs[1])
end
