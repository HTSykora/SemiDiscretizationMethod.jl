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
    #rhat = maximum([r, p - 1])+10
    #rhat = maximum([maximum(j), p - 1])
    #@show (maximum(i), maximum(j), p , (rhat + p + 1) )
    #@show (maximum(i), maximum(j), p *d, (rhat + p + 1)*d )
    #@show i
    #@show j
    PHI = sparse(i, j, v, p * d, (rhat + p + 1) * d)
    dropzeros!(PHI)

    #PHIL = -PHI[1:(rhat+1)*d, 1:(rhat+1)*d]
    #PHIR = PHI[1:(rhat+1)*d, (rhat+1)*d+1:end]

    #Note: sign is changed according to the Journal Publication
    PHIL = PHI[1:(p)*d, 1:(p)*d]
    PHIR = -PHI[1:(p)*d, (p)*d+1:end]

    PHILL = vcat(
        hcat(PHIL, spzeros((p) * d, (rhat + 1 - p) * d)),
        hcat(spzeros((rhat + 1 - p) * d, (p) * d), spdiagm(0 => ones((rhat + 1 - p) * d)))
    )
    PHIRR = vcat(
        PHIR,
        hcat(spdiagm(0 => ones((rhat + 1 - p) * d)), spzeros((rhat + 1 - p) * d, (p) * d))
    )
    # Wrong order of the blocks!
    # mappingVs = [spzeros(size(PHILL, 1))]
    # for (it, subV) in enumerate(rst.subVs)
    #     # pos=(it-1)*d
    #     # mappingVs[1][(1+pos) : (d+pos)] .+= subV.V
    #     mappingVs[1][(1+(it-1)*d):(d+(it-1)*d)] .+= subV.V
    # end

    # Direction of the excitation vector is reversed accordingly
    mappingVs = [spzeros(size(PHILL, 1))]
    p_steps = rst.n_steps
    for (it, subV) in enumerate(rst.subVs)
        target_block = p_steps - it + 1
        pos = (target_block - 1) * d
        mappingVs[1][(1+pos):(d+pos)] .+= subV.V
    end
    ([rst.ts[1], rst.ts[end]], PHILL, PHIRR, mappingVs)
end

function spectralRadiusOfMapping(mappLR::DiscreteMapping_LR{tT,mxT,vT}; args...)::mxT.parameters[1] where {tT,mxT,vT}
    return abs(eigs(mappLR.RmappingMX, mappLR.LmappingMX; args...)[1][1])::mxT.parameters[1]
end

function fixPointOfMapping(mappLR::DiscreteMapping_LR{tT,mxT,vT}) where {tT,mxT,vT}
    (mappLR.LmappingMX - mappLR.RmappingMX) \ Vector(mappLR.mappingVs[1])
end



"""
    extract_SDM_system(dde_rhs, p, ::Val{N}; t_test = 0.0) where {N}

Automatically extracts the necessary matrices and additive term for 
SemiDiscretizationMethod.jl from a black-box DDE "right-hand-side" function.
This is an experimetal tool, only for finite number of delays (definitly different at all times).
It migh fail for complex systems.
"""
function extract_SDM_system(dde_rhs, p, ::Val{N}; t_test = 0.0) where {N}
    
    u_eq = @SVector zeros(Float64, N)
    
    # --- Step 1: Detect delays using a "Spy function" ---
    requested_times = Float64[]
    h_spy(p_arg, t_req; kwargs...) = begin
        push!(requested_times, t_req)
        return u_eq
    end
    
    # Run at the test time to catch the requested past times
    dde_rhs(u_eq, h_spy, p, t_test)
    
    # Calculate tau values, round them to avoid minor numerical errors, and filter
    taus = sort(unique(round.(t_test .- requested_times, digits=10)))
    filter!(tau -> tau > 1e-10, taus) # Filter out tau=0 requests, if any
    
    # --- Step 2: Set up ForwardDiff Duals for Jacobian calculation ---
    # N-dimensional perturbation for gradient/Jacobian extraction
    partials = ntuple(i -> ForwardDiff.Partials{N, Float64}(ntuple(j -> i==j ? 1.0 : 0.0, N)), N)
    A_duals = SVector{N}(ntuple(i -> ForwardDiff.Dual{nothing, Float64, N}(0.0, partials[i]), N))
    
    u_dual = u_eq .+ A_duals
    
    # Helper function to construct SMatrix from the Dual result
    function get_jac(du_dual)
        SMatrix{N, N, Float64, N*N}(ntuple(k -> begin
            i = (k - 1) % N + 1
            j = (k - 1) ÷ N + 1
            ForwardDiff.partials(du_dual[i], j)
        end, N*N))
    end

    # --- Step 3: Extract A(t) - Present time Jacobian matrix ---
    A_func = t -> begin
        # Freeze the past: no dual perturbation in the past
        h_const(p_arg, t_req; kwargs...) = u_eq 
        du_out = dde_rhs(u_dual, h_const, p, t)
        return get_jac(du_out)
    end
    AMx = ProportionalMX(A_func)

    # --- Step 4: Extract B_i(t) - Delayed Jacobian matrices ---
    # Fix: Strongly typed array for DelayMX to avoid MethodError in LDDEProblem
    delay_matrices = DelayMX{N}[]
    for tau in taus
        B_func = t -> begin
            # The past receives a perturbed (Dual) value only at the specific tau time
            h_tau(p_arg, t_req; kwargs...) = begin
                if isapprox(t - t_req, tau, atol=1e-8)
                    return u_dual
                else
                    return u_eq
                end
            end
            du_out = dde_rhs(u_eq, h_tau, p, t)
            return get_jac(du_out)
        end
        # Create the DelayMX object for the given delay
        push!(delay_matrices, DelayMX(t -> tau, B_func))
    end

    # --- Step 5: Extract c(t) - Additive term ---
    c_func = t -> begin
        # Since u_eq = 0, evaluating the rhs alone yields the additive c(t) term
        h_zero(p_arg, t_req; kwargs...) = u_eq
        du_out = dde_rhs(u_eq, h_zero, p, t)
        return SVector{N, Float64}(ForwardDiff.value.(du_out))
    end
    cVec = Additive(c_func)

    # Return the constructed system
    return LDDEProblem(AMx, delay_matrices, cVec)
end

