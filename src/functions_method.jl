###############################################################################
###################### Submatrix from the delayed terms #######################
###############################################################################
#τ discretisation
function (method::SemiDiscretization)(τ::Real, rst::AbstractResult{d}) where d
    r = rOfDelay(τ, method)
    τerr = τ - r * method.Δt
    ranges = [subMxRange(r - k, d) for k in 0:methodorder(method)]
    return (τerr, ranges)
end
function (method::SemiDiscretization)(τs::Vector{<:Real}, rst::AbstractResult{d}) where d
    rs = rOfDelay.(τs, Ref(method))
    τerrs = τs .- rs .* method.Δt
    rangess = [[subMxRange(r - k, d) for k in 0:methodorder(method)] for r in rs]
    return (τerrs, rangess)
end

###############################################################################
# To distinguish between B as a constant matrix vs B as a function matrix

function (method::SemiDiscretization{NumericSD{N}})(τ::Real, B::DelayMX{d,dT,<:AbstractMatrix}, rst::AbstractResult) where {d,dT,N}
    (τerr, ranges) = method(τ, rst)
    h = method.Δt
    ord = methodorder(method)
    
    # 3-point Gauss-Legendre weights
    wts = (h/2 .* [5/9, 8/9, 5/9])
    pts_std = (1 .- [sqrt(3/5), 0.0, -sqrt(3/5)])

    # Check for periodicity
    if rst.expA_pts isa CyclicVector && B.T > 0 && isapprox(B.T, rst.LDDEP.A.T, atol=1e-8)
        n_period = rst.expA_pts.s
        MXs = [Vector{typeof(B.MX)}(undef, ord + 1) for _ in 1:n_period]
        for i in 1:n_period
            Es = rst.expA_pts[i]
            for k in 0:ord
                val = zero(B.MX)
                for j in 1:length(wts)
                    l_val = lagr_el0(ord, k, τerr, h, (h/2 * pts_std[j]))
                    val += wts[j] * (Es[j] * l_val)
                end
                MXs[i][k + 1] = val * B.MX
            end
        end
        return CyclicVector([SubMX(ranges, m) for m in MXs], rst.n_steps)
    end

    # Generic loop
    MXs = [Vector{typeof(B.MX)}(undef, ord + 1) for _ in 1:rst.n_steps]
    for i in 1:rst.n_steps
        Es = rst.expA_pts[i]
        for k in 0:ord
            val = zero(B.MX)
            for j in 1:length(wts)
                l_val = lagr_el0(ord, k, τerr, h, (h/2 * pts_std[j]))
                val += wts[j] * (Es[j] * l_val)
            end
            MXs[i][k + 1] = val * B.MX
        end
    end
    return SubMX.(Ref(ranges), MXs)
end

function (method::SemiDiscretization{NumericSD{N}})(τs::Vector{<:Real}, B::DelayMX{d,dT,<:AbstractMatrix}, rst::AbstractResult) where {d,dT,N}
    (τerrs, rangess) = method(τs, rst)
    h = method.Δt
    ord = methodorder(method)
    wts = (h/2 .* [5/9, 8/9, 5/9])
    pts_std = (1 .- [sqrt(3/5), 0.0, -sqrt(3/5)])

    results = Vector{SubMX{typeof(B.MX)}}(undef, rst.n_steps)
    for i in 1:rst.n_steps
        τerr = τerrs[i]
        mxs_i = Vector{typeof(B.MX)}(undef, ord + 1)
        Es = rst.expA_pts[i]
        for k in 0:ord
            val = zero(B.MX)
            for j in 1:length(wts)
                l_val = lagr_el0(ord, k, τerr, h, (h/2 * pts_std[j]))
                val += wts[j] * (Es[j] * l_val)
            end
            mxs_i[k + 1] = val * B.MX
        end
        results[i] = SubMX(rangess[i], mxs_i)
    end
    return results
end

function (method::SemiDiscretization{<:NumericSD})(τ::Real, B::DelayMX{d,dT,<:Function}, rst::AbstractResult) where {d,dT}
    (τerr, ranges) = method(τ, rst)
    MXs = [[QuadGK.quadgk((t -> (exp(rst.A_avgs[i] * (rst.ts[i + 1] - t)) * B(t) * lagr_el0(methodorder(method), k, τerr, method.Δt, t-rst.ts[i]))), rst.ts[i], rst.ts[i + 1])[1] for k in 0:methodorder(method)] for i in 1:rst.n_steps]
    SubMX.(Ref(ranges), MXs)
end
function (method::SemiDiscretization{<:NumericSD})(τs::Vector{<:Real}, B::DelayMX{d,dT,<:Function}, rst::AbstractResult) where {d,dT}
    (τerrs, rangess) = method(τs, rst)
    MXs = [[QuadGK.quadgk((t -> (exp(rst.A_avgs[i] * (rst.ts[i + 1] - t)) * B(t) * lagr_el0(methodorder(method), k, τerrs[i], method.Δt, t-rst.ts[i]))), rst.ts[i], rst.ts[i + 1])[1] for k in 0:methodorder(method)] for i in 1:rst.n_steps]
    SubMX.(rangess, MXs)
end

###############################################################################
# To distinguish between τ as constant vs τ as function
function (method::SemiDiscretization{<:NumericSD})(B::DelayMX{d,<:Real,mT}, rst::AbstractResult) where {d,mT}
    method(B.τ.τ, B, rst)
end

function (method::SemiDiscretization{<:NumericSD})(B::DelayMX{d,<:Function,mT}, rst::AbstractResult) where {d,mT}
    τis = [quadgk(B.τ, rst.ts[i], rst.ts[i + 1])[1] / method.Δt for i in 1:rst.n_steps]
    method(τis, B, rst)
end

###############################################################################
###################### Submatrix from the additive terms ######################
###############################################################################
function (method::SemiDiscretization)(c::Additive{d,<:AbstractArray{<:Real}}, rst::AbstractResult{d}) where d
    h = method.Δt
    wts = (h/2 .* [5/9, 8/9, 5/9])
    
    # Check for periodicity
    if rst.expA_pts isa CyclicVector && c.T > 0 && isapprox(c.T, rst.LDDEP.A.T, atol=1e-8)
        n_period = rst.expA_pts.s
        Vs = Vector{SubV{SVector{d,eltype(rst.A_avgs[1])}}}(undef, n_period)
        for i in 1:n_period
            Es = rst.expA_pts[i]
            val = zeros(SVector{d,eltype(rst.A_avgs[1])})
            for j in 1:length(wts)
                val += wts[j] * (Es[j] * c.V)
            end
            Vs[i] = SubV(val)
        end
        return CyclicVector(Vs, rst.n_steps)
    end

    Vs = Vector{SubV{SVector{d,eltype(rst.A_avgs[1])}}}(undef, rst.n_steps)
    for i in 1:rst.n_steps
        Es = rst.expA_pts[i]
        val = zeros(SVector{d,eltype(rst.A_avgs[1])})
        for j in 1:length(wts)
            val += wts[j] * (Es[j] * c.V)
        end
        Vs[i] = SubV(val)
    end
    return Vs
end
function (method::SemiDiscretization)(c::Additive{d,<:Function}, rst::AbstractResult{d}) where d
    Vs = SubV.([SVector{d}(vec(quadgk(t -> exp(rst.A_avgs[i] * (rst.ts[i + 1] - t)) * c(t), rst.ts[i], rst.ts[i + 1])[1])) for i in 1:rst.n_steps])
end
function (method::SemiDiscretization)(c::Additive{d,<:Array{<:Nothing}}, rst::AbstractResult) where d
    Vs = fill(SubV(@SArray zeros(d)),length(rst.ts)-1)
end
