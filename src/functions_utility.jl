@doc """
    rOfDelay(τ::Float64,dt::Float64,ord::Int64)  
    rOfDelay(τ::Float64,dt::Float64,ord::ORDER)

Calculate the delay resolution for given delay `τ`, time resolution `dt` and order `ord` defined by an integer or `ORDER` type.

## Arguments
- `τ::Float64`: time delay to discretise
- `dt::Float64`: time resolution (length)
- `ord::Int64`/`ord::ORDER`: order of the discretisation
""" rOfDelay
function rOfDelay(τ::Real, dt::Real, ord::Integer)
    Int64((τ + 100eps(τ)) ÷ dt + ord ÷ 2)
  # trunc(Int64,τ/dt+ord/2.+1000*eps())
end
function rOfDelay(τ::Real, ord::Union{FullDiscretization,SemiDiscretization})
    rOfDelay(τ, ord.Δt, methodorder(ord))
end
function nStepOfLength(τ::Real,dt::Real)
    Int64((τ + 100eps(τ)) ÷ dt)
end
@doc """
    subMXrange(i::Int64,d::Int64)

Return the index range of the submatrix in the large discretisation matrix (used for constructing the deterministic matrix)

## Arguments
- `i::Int64`: the index of the state space vector the submatrix is multiplied with
- `d::Int64`: dimension of the state space
""" subMxRange
function subMxRange(i::Int64, d::Int64)
    (1:d, i * d + 1:(i + 1) * d)
end
@doc """
    subMXArray(i::Int64,d::Int64)

Return the index range of the submatrix in the large discretisation matrix as a 1 dimensional array (used for constructing the stochastic matrix)

## Arguments
- `i::Int64`: the index of the state space vector the submatrix is multiplied with
- `d::Int64`: dimension of the state space
""" subMXArray
function subMXArray(i::Int64, d::Int64)
    sMXrng = subMxRange(i, d);
    [(i, j) for i in sMXrng[1],j in sMXrng[2]]
end

@doc """
    subVecRange(d::Int64)

Return the index range of the subvector in the large discretised additive vector (it is always added to the present state)

## Arguments
- `d::Int64`: dimension of the state space
""" subVecRange
function subVecRange(d::Int64)
    (1:nDim, 1:1)
end
@doc """
    subVecArray(d::Int64)

Return the index range of the subvector in the large discretised additive vector as a 1 dimensional array (it is always added to the present state)

## Arguments
- `d::Int64`: dimension of the state space
""" subVecArray
function subVecArray(d::Int64)
    sMXrng = subVecRange(d);
    [(i, j) for i in sMXrng[1],j in sMXrng[2]]
end

@doc """
    addSubmatrixToResult!(result::AbstractArray,subMX::detSubMX)

Add the deterministic submatrix `subMX` to the deterministic mapping matrix
Return the index range of the subvector in the large discretised additive vector as a 1 dimensional array (it is always added to the present state)

## Arguments
- `result::AbstractArray`: mapping matrix to add the subMX to
- `subMX::detSubMX`: submatrix to add to the mapping matrix
""" addSubmatrixToResult!
function addSubmatrixToResult!(result::AbstractArray, subMX::SubMX)
    for i in eachindex(subMX.ranges, subMX.MXs)
        result[subMX.ranges[i]...] .+= subMX.MXs[i]
    end
end
function addSubvectorToResults!(results::AbstractArray, subV::SubV)
    results[1:size(subV.V, 1)] .+= subV.V
end
# function addSubvectorToResults!(results::AbstractArray, subV::SubV)
#     d = size(subV.Vs[1], 1)
#     for j in eachindex(subV.Vs, results)
#         for i in eachindex(subV.Vs[j])
#             results[j][1:d] .+= subV.Vs[j][:]
#         end
#     end
# end

# Towards higher order ############################################################
# Building the Lagrangian polynomial for higher order methods
function lagr_atom(i::Real, k::Real, l::Real, τ::Real, r::Real, h::Real, t::Real)
    ((t - τ) - (i - r) * h - l * h) / ((k - l) * h)
end

function lagr_el(q::Real, i::Real, k::Real, τ::Real, r::Real, h::Real, t::Real)
    range = [0:k - 1; k + 1:q]
    prod([lagr_atom(i, k, l, τ, r, h, t) for l in range])
end

function lagr_atom0(k::Integer, l::Integer, τerr::Real, h::Real, t::Real)
    return (t - l * h - τerr) / ((k - l) * h)
end

function lagr_el0(q::Integer, k::Integer, τerr::Real, h::Real, t::Real)
    if q == 1
        if k == 0
            return (h + τerr - t) / h
        else
            return (t - τerr) / h
        end
    elseif q == 0
        return 1.0
    else
        val = 1.0
        for l in 0:q
            if l != k
                val *= (t - l * h - τerr) / ((k - l) * h)
            end
        end
        return val
    end
end

# Handling multiple time-steps
# Reduction of the mapping matrices and -vectors on a time interval 
function prodl(mappingMXs::Vector{TM}) where TM<:AbstractMatrix
    mx0 = mappingMXs[1]
    for i in 2:length(mappingMXs)
        mx0 = mappingMXs[i] * mx0
    end
    return mx0
end
function prodl(mappingMXs::Vector{TM},idxs::AbstractVector{<:Integer}) where TM<:AbstractMatrix
    mx0 = mappingMXs[1][idxs,idxs]
    for i in 2:length(mappingMXs)
        mx0 = mappingMXs[i][idxs,idxs] * mx0
    end
    return mx0
end
# function prodl(mappingMXs::Vector{TM}) where TM<:AbstractMatrix
#     prod(reverse(mappingMXs))
# end
function reduce_additive(mappingMXs::Vector{TM}, mappingVs::Vector{TV}) where {TM<:AbstractMatrix,TV<:AbstractVector}
    mappingV = mappingVs[1];
    for i in 2:length(mappingVs)
        mappingV = mappingMXs[i] * mappingV + mappingVs[i]
    end
    return mappingV
end
function reduce_additive(mappingMXs::Vector{TM}, mappingVs::Vector{TV}, idxs::AbstractVector{<:Integer}) where {TM<:AbstractMatrix,TV<:AbstractVector}
    mappingV = mappingVs[1][idxs];
    for i in 2:length(mappingVs)
        mappingV = mappingMXs[i][idxs,idxs] * mappingV + mappingVs[i][idxs]
    end
    return mappingV
end


function reduce_additive(mappingMX0::TM, mappingV0::TV, q::Int) where {TM<:AbstractMatrix,TV<:AbstractVector}
    mappingV = mappingV0;
    for i in 2:q
        mappingV = mappingMX0 * mappingV + mappingV0
    end
    return mappingV
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