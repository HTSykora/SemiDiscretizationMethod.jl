# Independent verification of DiscreteMapping_LR against DiscreteMapping and package-free references:
#   S1  LTI 1DOF turning, constant τ::Real, T = τ, 2τ, 3τ: ρ vs rightmost characteristic root,
#       fixpoint vs the exact forced periodic response and vs a long RK4 simulation
#   S2  delayed Mathieu, principal period T = 2τ, 3τ: ρ vs an RK4 monodromy matrix,
#       fixpoint vs a long forced RK4 simulation
#   S3  time-periodic delay without time-reversal symmetry (single point), SD order 1 and 2
# Usage: julia --project=<env with SemiDiscretizationMethod> examples/verify_LR_mapping.jl [s1] [s2] [s3]
using SemiDiscretizationMethod, LinearAlgebra, SparseArrays, StaticArrays, Printf
const SDM = SemiDiscretizationMethod

# ----------------------------------------------------------------- SD wrappers (dense eig)
ρ_classic(dm) = maximum(abs, eigvals(Matrix(SDM.prodl(dm.mappingMXs))))
ρ_LR(mp) = maximum(abs, eigvals(Matrix(mp.LmappingMX) \ Matrix(mp.RmappingMX)))
# both fixpoints are ordered [x(0); x(-Δt); x(-2Δt); ...]
fix_classic(dm) = Vector(fixPointOfMapping(dm))
fix_LR(mp) = Vector(fixPointOfMapping(mp))
blk(v, d, k) = v[k*d+1:(k+1)*d]

function run_sd(prob, order, Δt, τmax, p)
    method = SemiDiscretization(order, Δt)
    dm = DiscreteMapping(prob, method, τmax; n_steps = p, calculate_additive = true)
    out = (ρc = ρ_classic(dm), zc = fix_classic(dm), ρl = NaN, zl = Float64[], err = "")
    try
        mp = DiscreteMapping_LR(prob, method, τmax; n_steps = p, calculate_additive = true)
        out = merge(out, (ρl = ρ_LR(mp), zl = fix_LR(mp)))
    catch e
        out = merge(out, (err = first(split(sprint(showerror, e), '\n')),))
    end
    out
end
fixdiff(r, d, K) = isempty(r.zl) ? NaN : maximum(abs, r.zl[1:K*d] - r.zc[1:K*d])

# ----------------------------------------------------------------- reference 1: RK4 DDE stepper
# Fixed step h, pointwise delay τ(t), delayed state by 6-point Lagrange interpolation of the stored
# grid (stencil checked to lie in the known past). Works on d×M matrix states: M = 1 for a forced
# simulation, M = d(K+1) with a unit-history basis for the monodromy matrix.
struct DDE{FA,FB,FT,FC}
    A::FA
    B::FB
    τ::FT
    c::FC
end

function rk4_dde(sys::DDE, h, K, nsteps, hist::Vector{Matrix{Float64}}; forced = true)
    X = Vector{Matrix{Float64}}(undef, K + 1 + nsteps)
    for k in 0:K
        X[K+1-k] = hist[k+1]
    end
    I0 = K + 1                                   # X[I0 + j] = X(j h)
    function delayed(t, nknown)
        s = t - sys.τ(t)
        j0 = floor(Int, s / h)
        θ = s / h - j0
        (j0 - 2 >= -K && j0 + 3 <= nknown) || error("interpolation stencil out of range")
        acc = zero(X[I0])
        for m in -2:3
            w = 1.0
            for l in -2:3
                l == m || (w *= (θ - l) / (m - l))
            end
            acc .+= w .* X[I0+j0+m]
        end
        acc
    end
    f(t, Y, n) = forced ? (sys.A(t) * Y + sys.B(t) * delayed(t, n) .+ sys.c(t)) : (sys.A(t) * Y + sys.B(t) * delayed(t, n))
    for n in 0:nsteps-1
        t, Y = n * h, X[I0+n]
        k1 = f(t, Y, n)
        k2 = f(t + h / 2, Y + h / 2 * k1, n)
        k3 = f(t + h / 2, Y + h / 2 * k2, n)
        k4 = f(t + h, Y + h * k3, n)
        X[I0+n+1] = Y + h / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    X, I0
end

function ρ_reference(sys, d, T, nT, τmax)
    h = T / nT
    K = ceil(Int, τmax / h) + 6
    M = d * (K + 1)
    Id = Matrix{Float64}(I, M, M)
    hist = [Id[k*d+1:(k+1)*d, :] for k in 0:K]
    X, I0 = rk4_dde(sys, h, K, nT, hist; forced = false)
    U = reduce(vcat, [X[I0+nT-k] for k in 0:K])
    maximum(abs, eigvals(U))
end

# long forced simulation from zero history; returns x(t_end - k Δt), k = 0..Kout-1
function longsim(sys, d, T, Δt, m, nper, τmax, Kout)
    h = Δt / m
    nT = round(Int, T / h)
    K = ceil(Int, τmax / h) + 6
    X, I0 = rk4_dde(sys, h, K, nper * nT, [zeros(d, 1) for _ in 0:K])
    N = nper * nT
    reduce(vcat, [vec(X[I0+N-k*m]) for k in 0:Kout-1]), reduce(vcat, [vec(X[I0+N-nT-k*m]) for k in 0:Kout-1])
end

# ----------------------------------------------------------------- S1: LTI 1DOF turning
# x'' + 2ζx' + (1+w)x = w x(t-τ) + f(t),   f(t) = f0 + f1 cos(2πt/τ) + f2 sin(2πt/T + φ)
function s1(; ζ = 0.05, w = 0.2, τ = 2π, f0 = 1.0, f1 = 0.5, f2 = 0.3, φ = 0.3)
    println("="^100)
    println("S1  LTI turning, constant τ::Real, T = kτ, excitation with periods τ and T")
    A = @SMatrix [0.0 1.0; -(1 + w) -2ζ]
    B = @SMatrix [0.0 0.0; w 0.0]
    # rightmost characteristic root of λ² + 2ζλ + 1 + w − w e^{−λτ} = 0 (Newton from a seed grid)
    D(λ) = λ^2 + 2ζ * λ + 1 + w - w * exp(-λ * τ)
    dD(λ) = 2λ + 2ζ + w * τ * exp(-λ * τ)
    roots = ComplexF64[]
    for re in -0.6:0.1:0.3, im in 0.0:0.05:6.0
        λ = complex(re, im)
        for _ in 1:60
            λ -= D(λ) / dD(λ)
        end
        abs(D(λ)) < 1e-12 && push!(roots, λ)
    end
    λmax = roots[argmax(real.(roots))]
    @printf("rightmost characteristic root λ = %.12f %+.12fi\n", real(λmax), imag(λmax))
    H(s) = inv(s * I - Matrix(A) - Matrix(B) * exp(-s * τ))
    cT(T) = t -> @SVector [0.0, f0 + f1 * cos(2π * t / τ) + f2 * sin(2π * t / T + φ)]
    # exact periodic response (superposition of the three harmonics)
    xp(t, T) = real(-(Matrix(A + B)) \ [0.0, f0] + H(2π / τ * im) * [0.0, f1] * exp(2π / τ * im * t) +
                    H(2π / T * im) * [0.0, -im * f2 * exp(im * φ)] * exp(2π / T * im * t))
    for k in (1, 2, 3)
        T = k * τ
        prob = LDDEProblem(ProportionalMX(A), [DelayMX(τ, B)], Additive(cT(T)))
        ρex = exp(real(λmax) * T)
        for N in (100, 200)
            Δt = τ / N
            r = run_sd(prob, 1, Δt, τ, k * N)
            Kc = N + 1
            zex = reduce(vcat, [xp(-j * Δt, T) for j in 0:Kc-1])
            @printf("T=%dτ Δt=τ/%d | ρ: exact %.8f classic %.8f LR %s | rel.err %.1e | fixpt |classic-exact| %.1e |LR-classic| %s %s\n",
                k, N, ρex, r.ρc, isnan(r.ρl) ? "  ERROR   " : @sprintf("%.8f", r.ρl), abs(r.ρc - ρex) / ρex,
                maximum(abs, r.zc[1:2Kc] - zex), isnan(r.ρl) ? "-" : @sprintf("%.1e", fixdiff(r, 2, Kc)), r.err)
        end
    end
    # long simulation check of the T = 2τ fixpoint
    T, N = 2τ, 100
    sys = DDE(t -> A, t -> B, t -> τ, cT(T))
    xs, xs_prev = longsim(sys, 2, T, τ / N, 4, 120, τ, N + 1)
    zex = reduce(vcat, [xp(-j * τ / N, T) for j in 0:N])
    r = run_sd(LDDEProblem(ProportionalMX(A), [DelayMX(τ, B)], Additive(cT(T))), 1, τ / N, τ, 2N)
    @printf("long sim (RK4 h=Δt/4, 120 periods of T=2τ): |sim-exact| %.1e, period-to-period drift %.1e, |SD classic-sim| %.1e, |SD LR-sim| %s\n",
        maximum(abs, xs - zex), maximum(abs, xs - xs_prev), maximum(abs, r.zc[1:length(xs)] - xs),
        isempty(r.zl) ? "ERROR" : @sprintf("%.1e", maximum(abs, r.zl[1:length(xs)] - xs)))
end

# ----------------------------------------------------------------- S2: delayed Mathieu, T = kτ
# x'' + κx' + (δ + ε cos(2πt/T))x = b x(t-τ) + sin(2πt/T + 0.4) + 0.3
function s2(; δ = 3.0, ε = 2.0, b = -0.15, κ = 0.1, τ = 2π)
    println("="^100)
    println("S2  delayed Mathieu, constant τ::Real, principal period T = kτ, excitation with period T")
    B = @SMatrix [0.0 0.0; b 0.0]
    for k in (2, 3)
        T = k * τ
        Af = t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π * t / T) -κ]
        c = t -> @SVector [0.0, sin(2π * t / T + 0.4) + 0.3]
        prob = LDDEProblem(ProportionalMX(Af), [DelayMX(τ, B)], Additive(c))
        sys = DDE(Af, t -> B, t -> τ, c)
        ρref = ρ_reference(sys, 2, T, k * 400, τ)
        ρref2 = ρ_reference(sys, 2, T, k * 800, τ)
        @printf("T=%dτ reference ρ (RK4 monodromy, h=τ/400 | τ/800): %.10f | %.10f\n", k, ρref, ρref2)
        for N in (100, 200)
            Δt = τ / N
            r = run_sd(prob, 1, Δt, τ, k * N)
            xs, xs_prev = longsim(sys, 2, T, Δt, 4, 60, τ, N + 1)
            @printf("  Δt=τ/%d | ρ classic %.8f LR %s rel.err %.1e | fixpt |classic-sim| %.1e |LR-classic| %s (sim drift %.0e) %s\n",
                N, r.ρc, isnan(r.ρl) ? "  ERROR   " : @sprintf("%.8f", r.ρl), abs(r.ρc - ρref2) / ρref2,
                maximum(abs, r.zc[1:length(xs)] - xs), isnan(r.ρl) ? "-" : @sprintf("%.1e", fixdiff(r, 2, N + 1)),
                maximum(abs, xs - xs_prev), r.err)
        end
    end
end

# ----------------------------------------------------------------- S3: time-periodic delay, no time symmetry
# τ(t) = τ0 (1 + a (sin νt + 0.5 sin(2νt + 1))), ν = 2π/T, T = 2.5 τ0 (non-integer T/τ0)
# The two harmonics with phase 1 rad have no common symmetry axis or centre, so τ(t) is neither
# mirror symmetric (τ(c+t) = τ(c−t)) nor point symmetric; the forcing likewise.
function s3(; ζ = 0.05, w = 0.15, τ0 = 2π, a = 0.15)
    println("="^100)
    println("S3  turning with time-periodic, time-asymmetric delay τ(t), T = 2.5 τ0, single parameter point")
    T = 2.5τ0
    ν = 2π / T
    τf = t -> τ0 * (1 + a * (sin(ν * t) + 0.5 * sin(2ν * t + 1.0)))
    tg = range(0, T, length = 10001)
    τmax = maximum(τf, tg) * (1 + 1e-9)
    A = @SMatrix [0.0 1.0; -(1 + w) -2ζ]
    B = @SMatrix [0.0 0.0; w 0.0]
    c = t -> @SVector [0.0, 1.0 + sin(ν * t) + 0.4 * cos(2ν * t + 0.7)]
    prob = LDDEProblem(ProportionalMX(A), [DelayMX(τf, B)], Additive(c))
    sys = DDE(t -> A, t -> B, τf, c)
    asym = minimum(maximum(abs(τf(t0 + s) - τf(t0 - s)) for s in range(0, T / 2, length = 200)) for t0 in range(0, T, length = 2001))
    @printf("τ(t) ∈ [%.4f, %.4f] s; min over axes t0 of max_s |τ(t0+s) − τ(t0−s)| = %.3f s (> 0: no mirror symmetry)\n",
        minimum(τf, tg), τmax, asym)
    ρref = ρ_reference(sys, 2, T, 1000, τmax)
    ρref2 = ρ_reference(sys, 2, T, 2000, τmax)
    @printf("reference ρ (RK4 monodromy, h=T/1000 | T/2000): %.10f | %.10f\n", ρref, ρref2)
    for order in (1, 2), p in (250, 500)
        Δt = T / p
        r = run_sd(prob, order, Δt, τmax, p)
        Kc = SDM.rOfDelay(τmax, Δt, order) + 1
        xs, xs_prev = longsim(sys, 2, T, Δt, 4, 80, τmax, min(p, Kc))
        @printf("  order %d Δt=T/%d | ρ classic %.8f LR %s rel.err %.1e | fixpt |classic-sim| %.1e |LR-classic| %s (sim drift %.0e) %s\n",
            order, p, r.ρc, isnan(r.ρl) ? "  ERROR   " : @sprintf("%.8f", r.ρl), abs(r.ρc - ρref2) / ρref2,
            maximum(abs, r.zc[1:length(xs)] - xs), isnan(r.ρl) ? "-" : @sprintf("%.1e", fixdiff(r, 2, min(p, Kc))),
            maximum(abs, xs - xs_prev), r.err)
    end
end

const SECTIONS = isempty(ARGS) ? ["s1", "s2", "s3"] : ARGS
"s1" in SECTIONS && s1()
"s2" in SECTIONS && s2()
"s3" in SECTIONS && s3()
