struct CyclicVector{T} <: AbstractVector{T}
    s::Int64 # period
    n::Int64 # total length
    V::Vector{T}
end
CyclicVector(V::Vector{T}, n::Int64) where T = CyclicVector{T}(length(V), n, V)
CyclicVector(V::Vector{T}) where T = CyclicVector{T}(length(V), length(V), V)
Base.getindex(cV::CyclicVector, idx) = cV.V[((idx .- 1) .% cV.s .+ 1)]
Base.size(LV::CyclicVector) = (LV.n,)
Base.size(LV::CyclicVector, ::Integer) = LV.n

abstract type subArray{mT} end
struct SubMX{mT} <: subArray{mT}
    ranges::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}}
    MXs::Vector{mT}
end
SubMX(range::Tuple{<:AbstractVector,<:AbstractVector}, MX::SMatrix{<:Real}) = SubMX([range], [MX])
SubMX(range::Tuple{<:AbstractVector,<:AbstractVector}, MX::AbstractMatrix{<:Real}) = SubMX([range], [MX])


struct SubV{vT} <: subArray{vT}
    V::vT
end

abstract type AbstractResult{d} end
struct Result{d,lddepT<:AbstractLDDEProblem{d},mT,submxT,subvT,tT,AavgT,expAT,expAPtsT} <: AbstractResult{d}
    LDDEP::lddepT
    method::mT #::DiscretizationMethod
    subMXs::Vector{Vector{submxT}} # [[A(t1),A(t2),...],[B1(t1),B1(t2),...],[B2(t1),B2(t2),...],...]
    subVs::Vector{subvT} # [c(t1),c(t2),...]
    #
    A_avgs::AavgT
    expA_hs::expAT # Pre-calculated matrix exponentials exp(A_avgs[i] * h)
    expA_pts::expAPtsT # Pre-calculated matrix exponentials at Gauss points for this step
    ts::Vector{tT} # [0,t1,t2,...]
    n_steps::Int64 # number of time steps
    n::Int64 # Large discretisation matrix size
    calculate_additive::Bool
    # d::Integer # state space dimension
end

function calculate_Aavgs(A::ProportionalMX{d,<:Function}, ts::AbstractVector{<:Real}, Δt::Real, n_steps::Int64) where d
    if A.T > 0 && n_steps > Int(round(A.T / Δt))
        n_period = Int(round(A.T / Δt))
        avgs = [A(ts[i] + Δt/2) for i in 1:n_period]
        return CyclicVector(avgs, n_steps)
    else
        return [A(ts[i] + Δt/2) for i in 1:n_steps]
    end
end
calculate_Aavgs(A::ProportionalMX{d,<:mT}, ts::AbstractVector{<:Real},Δt::Real, n_steps::Int64)  where d where mT <: AbstractMatrix{T} where T = CyclicVector([A.MX], n_steps)

function Result(LDDEP::LDDEProblem{d,AT,BT, cT}, method::DiscretizationMethod{fT}, DiscretizationLength::Real; n_steps::Int64=nStepOfLength(DiscretizationLength, method.Δt), calculate_additive::Bool=false) where {d,AT,BT,cT,fT}
    # DiscretizationLength discretisated time interval length
    # n_steps: how many mapping matrix to calculate
    ts = collect(range(zero(fT), length=n_steps + 1, step=method.Δt))
    n = (rOfDelay(DiscretizationLength, method) + 1) * d
    A_avgs = calculate_Aavgs(LDDEP.A, ts, method.Δt, n_steps)
    
    h = method.Δt
    # 3-point Gauss-Legendre quadrature points for [0, h]
    pts = (h/2 .* (1 .- [sqrt(3/5), 0.0, -sqrt(3/5)]))

    # Pre-calculate matrix exponentials
    if A_avgs isa CyclicVector
        expA_hs = CyclicVector([exp(A * h) for A in A_avgs.V], n_steps)
        expA_pts = CyclicVector([[exp(A * (h - t)) for t in pts] for A in A_avgs.V], n_steps)
    else
        expA_hs = [exp(A * h) for A in A_avgs]
        expA_pts = [[exp(A * (h - t)) for t in pts] for A in A_avgs]
    end

    subMXs = [Vector{SubMX{eltype(A_avgs)}}(undef, n_steps) for i in 1:(length(LDDEP.Bs) + 1)] # []
    if calculate_additive
        # subVs = Vector{SubV{eltype(LDDEP.cT.V)}}(undef, n_steps) # []
        subVs = Vector{SubV{SVector{d,eltype(eltype(A_avgs))}}}(undef, n_steps) # []
    else
        # subVs = Vector{SubV{eltype(LDDEP.cT.V)}}(undef, 0)
        subVs = Vector{SubV{SVector{d,eltype(eltype(A_avgs))}}}(undef, 0)
    end
    
    Result(LDDEP, method, subMXs, subVs, A_avgs, expA_hs, expA_pts, ts, n_steps, n, calculate_additive)
end

struct DiscreteMapping{tT,mxT,vT}
    ts::Vector{tT}
    mappingMXs::Vector{mxT}
    mappingVs::Vector{vT}
end

function initializeOneStepMapping(dm::DiscreteMapping{tT,mxT,vT}) where {tT,mxT,vT}
    DiscreteMapping(Vector{tT}(undef,2),Vector{mxT}(undef,1),Vector{vT}(undef,1))
end
