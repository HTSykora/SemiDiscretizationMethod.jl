using SemiDiscretizationMethod
using Arpack
using KrylovKit
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using StaticArrays

# System 1: Standard Mathieu (d=2)
function get_mathieu_mapping(N)
    T = 20π
    τ = 2π
    AMx = ProportionalMX(t -> SMatrix{2,2,Float64}(0.0, -3.0 - 1.5 * cos(2π / T * t), 1.0, -0.2), T=T)
    BMx1 = DelayMX(t -> τ, SMatrix{2,2,Float64}(0.0, -0.5, 0.0, 0.0))
    prob = LDDEProblem(AMx, BMx1)
    method = SemiDiscretization(1, T / N)
    return prob, method, τ
end

function verify_logic(N::Int=500)
   
    println(" ------------------------------ Multiplcation free N=$N ---------------------")
    p, m, tau = get_mathieu_mapping(N)
    mapp = DiscreteMapping_LR(p, m, tau)

    println("Testing DiscreteMapping_LR with useKrylovKit=true (default)")
    @time mu1 = spectralRadiusOfMapping(mapp, useKrylovKit=true)
    @show mu1

    println("Testing DiscreteMapping_LR with useKrylovKit=false")
    @time mu2 = spectralRadiusOfMapping(mapp, useKrylovKit=false)
    @show mu2

    @assert isapprox(mu1, mu2, atol=1e-8)
    println("SUCCESS: Both solvers return the same spectral radius.")

    println(" ------------------------------ Traditional N=$N ---------------------")

    p, m, tau = get_mathieu_mapping(N)
    mapp = DiscreteMapping(p, m, tau)

    println("Testing DiscreteMapping with useKrylovKit=false (default)")
    @time mu1 = spectralRadiusOfMapping(mapp, useKrylovKit=false)
    @show mu1

    println("Testing DiscreteMapping with useKrylovKit=true")
    @time mu2 = spectralRadiusOfMapping(mapp, useKrylovKit=true)
    @show mu2

    @assert isapprox(mu1, mu2, atol=1e-8)
    println("SUCCESS: Both solvers return the same spectral radius.")


end

verify_logic(10000)
verify_logic(100)
