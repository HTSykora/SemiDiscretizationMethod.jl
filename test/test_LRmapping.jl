using SemiDiscretizationMethod
using Test
const SDM = SemiDiscretizationMethod

# Dense eigenvalues, so that the comparisons do not depend on the iterative solver tolerances.
ρ_classic(dm) = maximum(abs, eigvals(Matrix(SDM.prodl(dm.mappingMXs))))
ρ_LR(mp) = maximum(abs, eigvals(Matrix(mp.LmappingMX) \ Matrix(mp.RmappingMX)))

# 1DOF turning: x'' + 2ζx' + (1+w)x = w x(t-τ) + f(t)
turningA(w, ζ) = @SMatrix [0.0 1.0; -(1 + w) -2ζ]
turningB(w) = @SMatrix [0.0 0.0; w 0.0]

# Package-independent reference: RK4 with fixed step h, pointwise delay τ(t) and 6-point Lagrange
# interpolation of the stored grid. State is d×M (M = 1: forced simulation, M = d(K+1): monodromy).
function rk4_dde(Af, Bf, τf, cf, h, K, nsteps, hist; forced = true)
    X = Vector{Matrix{Float64}}(undef, K + 1 + nsteps)
    for k in 0:K
        X[K+1-k] = hist[k+1]
    end
    I0 = K + 1
    function delayed(t, n)
        s = t - τf(t)
        j0 = floor(Int, s / h)
        θ = s / h - j0
        @assert j0 - 2 >= -K && j0 + 3 <= n
        acc = zero(X[I0])
        for m in -2:3
            wgt = 1.0
            for l in -2:3
                l == m || (wgt *= (θ - l) / (m - l))
            end
            acc .+= wgt .* X[I0+j0+m]
        end
        acc
    end
    f(t, Y, n) = forced ? Af(t) * Y + Bf(t) * delayed(t, n) .+ cf(t) : Af(t) * Y + Bf(t) * delayed(t, n)
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

@testset "DiscreteMapping_LR" begin
    ζ, w, τ = 0.05, 0.2, 2π
    A, B = turningA(w, ζ), turningB(w)

    @testset "constant (Real) delay, T = $k τ, LTI vs characteristic root and exact forced response" for k in (1, 2, 3)
        N = 100
        T, Δt = k * τ, τ / N
        f0, f1, f2, φ = 1.0, 0.5, 0.3, 0.3
        c = t -> @SVector [0.0, f0 + f1 * cos(2π * t / τ) + f2 * sin(2π * t / T + φ)]
        prob = LDDEProblem(ProportionalMX(A), [DelayMX(τ, B)], Additive(c))
        method = SemiDiscretization(1, Δt)
        dm = DiscreteMapping(prob, method, τ; n_steps = k * N, calculate_additive = true)
        mp = DiscreteMapping_LR(prob, method, τ; n_steps = k * N, calculate_additive = true)

        @test ρ_LR(mp) ≈ ρ_classic(dm) rtol = 1e-10
        zc, zl = fixPointOfMapping(dm), fixPointOfMapping(mp)
        nc = length(zc)
        @test maximum(abs, zl[1:nc] - zc) < 1e-10

        # rightmost root of λ² + 2ζλ + 1 + w − w e^{−λτ} = 0
        D(λ) = λ^2 + 2ζ * λ + 1 + w - w * exp(-λ * τ)
        dD(λ) = 2λ + 2ζ + w * τ * exp(-λ * τ)
        λ = complex(0.0, 1.0)
        for _ in 1:50
            λ -= D(λ) / dD(λ)
        end
        @test abs(D(λ)) < 1e-12
        @test ρ_LR(mp) ≈ exp(real(λ) * T) rtol = 1e-3     # O(Δt²): 1.3e-4 … 4e-4 here

        H(s) = inv(s * I - Matrix(A) - Matrix(B) * exp(-s * τ))
        xp(t) = real(-(Matrix(A + B)) \ [0.0, f0] + H(2π / τ * im) * [0.0, f1] * exp(2π / τ * im * t) +
                     H(2π / T * im) * [0.0, -im * f2 * exp(im * φ)] * exp(2π / T * im * t))
        zex = reduce(vcat, [xp(-j * Δt) for j in 0:nc÷2-1])
        @test maximum(abs, zl[1:nc] - zex) < 1e-2         # O(Δt²): 3.3e-3 … 4.3e-3 here
    end

    @testset "DiscreteMappingSteps_LR does not mutate the shared Result" begin
        prob = LDDEProblem(ProportionalMX(A), [DelayMX(τ, B)])
        rst = SDM.calculateResults(prob, SemiDiscretization(1, τ / 40), τ; n_steps = 80)
        L1, R1 = SDM.DiscreteMappingSteps_LR(rst)[2:3]
        L2, R2 = SDM.DiscreteMappingSteps_LR(rst)[2:3]
        @test L1 == L2 && R1 == R2
        @test SDM.DiscreteMappingSteps(rst)[2] == DiscreteMapping(prob, SemiDiscretization(1, τ / 40), τ; n_steps = 80).mappingMXs
    end

    @testset "constant delay given as τ::Function equals τ::Real (N = $N, T = $k τ)" for N in (100, 101), k in (2, 3)
        method = SemiDiscretization(1, τ / N)
        probR = LDDEProblem(ProportionalMX(A), [DelayMX(τ, B)])
        probF = LDDEProblem(ProportionalMX(A), [DelayMX(t -> τ, B)])
        ρR = ρ_classic(DiscreteMapping(probR, method, τ; n_steps = k * N))
        @test ρ_classic(DiscreteMapping(probF, method, τ; n_steps = k * N)) ≈ ρR rtol = 1e-12
        @test ρ_LR(DiscreteMapping_LR(probF, method, τ; n_steps = k * N)) ≈ ρR rtol = 1e-12
        @test ρ_LR(DiscreteMapping_LR(probR, method, τ; n_steps = k * N)) ≈ ρR rtol = 1e-12
    end

    @testset "periodic coefficients with CyclicVector reuse (n_steps = 2 principal periods)" begin
        δ, ε, b0, κ, T = 3.0, 2.0, -0.15, 0.1, 2π
        Am = ProportionalMX(t -> @SMatrix([0.0 1.0; -δ-ε*cos(2π / T * t) -κ]); T = T)
        prob = LDDEProblem(Am, [DelayMX(τ, @SMatrix([0.0 0.0; b0 0.0]); T = T)], Additive([0.0, 1.0]; T = T))
        method = SemiDiscretization(1, T / 60)
        dm1 = DiscreteMapping(prob, method, τ; n_steps = 60, calculate_additive = true)
        dm2 = DiscreteMapping(prob, method, τ; n_steps = 120, calculate_additive = true)
        mp2 = DiscreteMapping_LR(prob, method, τ; n_steps = 120, calculate_additive = true)
        @test ρ_LR(mp2) ≈ ρ_classic(dm2) rtol = 1e-10
        @test ρ_LR(mp2) ≈ ρ_classic(dm1)^2 rtol = 1e-10
        z1, z2 = fixPointOfMapping(dm1), fixPointOfMapping(mp2)
        @test maximum(abs, z2[1:length(z1)] - z1) < 1e-10
    end

    @testset "time-periodic delay without time-reversal symmetry, order $order" for order in (1, 2)
        # τ(t) has no mirror axis and no centre of symmetry; T/τ0 = 2.5 is not an integer
        ζ3, w3, τ0, a = 0.05, 0.15, 2π, 0.15
        T = 2.5τ0
        ν = 2π / T
        τf = t -> τ0 * (1 + a * (sin(ν * t) + 0.5 * sin(2ν * t + 1.0)))
        τmax = maximum(τf, range(0, T, length = 10001)) * (1 + 1e-9)
        A3, B3 = turningA(w3, ζ3), turningB(w3)
        c = t -> @SVector [0.0, 1.0 + sin(ν * t) + 0.4 * cos(2ν * t + 0.7)]
        prob = LDDEProblem(ProportionalMX(A3), [DelayMX(τf, B3)], Additive(c))
        p = 250
        Δt = T / p
        method = SemiDiscretization(order, Δt)
        dm = DiscreteMapping(prob, method, τmax; n_steps = p, calculate_additive = true)
        mp = DiscreteMapping_LR(prob, method, τmax; n_steps = p, calculate_additive = true)
        ρl = ρ_LR(mp)
        @test ρl ≈ ρ_classic(dm) rtol = 1e-10
        zc, zl = fixPointOfMapping(dm), fixPointOfMapping(mp)
        @test maximum(abs, zl[1:length(zc)] - zc) < 1e-10

        # reference multiplier: RK4 monodromy on a history grid (converged to ~1e-8 at this h)
        h = T / 500
        K = ceil(Int, τmax / h) + 6
        M = 2 * (K + 1)
        Id = Matrix{Float64}(I, M, M)
        X, I0 = rk4_dde(t -> A3, t -> B3, τf, c, h, K, 500, [Id[2k+1:2k+2, :] for k in 0:K]; forced = false)
        ρref = maximum(abs, eigvals(reduce(vcat, [X[I0+500-k] for k in 0:K])))
        @test ρl ≈ ρref rtol = (order == 1 ? 5e-4 : 1e-4)   # 1.7e-4 / 2.7e-5 here

        # fixpoint vs long forced simulation from zero history (ρ ≈ 0.4: 40 periods suffice)
        m, nper = 2, 40
        hs = Δt / m
        Ks = ceil(Int, τmax / hs) + 6
        Xs, I0s = rk4_dde(t -> A3, t -> B3, τf, c, hs, Ks, nper * p * m, [zeros(2, 1) for _ in 0:Ks])
        nb = min(p, length(zc) ÷ 2)
        xs = reduce(vcat, [vec(Xs[I0s+nper*p*m-j*m]) for j in 0:nb-1])
        @test maximum(abs, zl[1:2nb] - xs) < (order == 1 ? 5e-4 : 2e-4)   # 1.1e-4 / 4.2e-5 here
    end
end
