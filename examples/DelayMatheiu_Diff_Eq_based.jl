# semidisc_dde_extractor.jl
# Compatible with the SemiDiscretizationMethod.jl package

using SemiDiscretizationMethod
using DifferentialEquations # Added for long-time DDE simulation
using ForwardDiff
using StaticArrays
using LinearAlgebra
using Plots
plotlyjs()
using LaTeXStrings

# =======================================================================
# 2. EXAMPLE: TESTING THE DELAY MATHIEU EQUATION
# =======================================================================
println("Configuring the system...")

# Parameters: δ, ε, b0, a1, T, τ
T_per = 2π 
τmax = T_per# 2π
params = (13.0, 0.2, -0.15, 0.1, T_per, τmax)

# Definition of the DDE "Black box" (standard rhs format)
# Equation: x'' + (δ + ε*cos(2π/T * t))x + a1*x' = b0*x(t-τ) + sin(4π/T * t)
function mathieu_dde_rhs(u, h, p, t)
    δ, ε, b0, a1, T, τ_val = p

    # Fetching the past
    hist = h(p, t - τ_val)
    tmod = mod(t, T)
    s = tmod < T / 5 ? 1 : 0
    du1 = u[2]
    du2 = sin(4π / T * t)
    du1 = u[2]
    du2 = -(δ + ε * cos(2π / T * t)) * u[1] - a1 * u[2] + b0 * hist[1] + s

    return @SVector [du1, du2]
end

# Automatic extraction of the system!
println("Automatically extracting system matrices using ForwardDiff...")
mathieu_lddep = extract_SDM_system(mathieu_dde_rhs, params, Val(2); t_test=0.0)

# Discretization method
method = SemiDiscretization(2, 0.05) # 2nd order discretization, Δt=0.05
n_steps_calc = Int((T_per + 100eps(T_per)) ÷ method.Δt)

# Calculate the mappings
@time mapping_LR = DiscreteMapping_LR(mathieu_lddep, method, τmax,
    n_steps=n_steps_calc, calculate_additive=true);

# Checking the spectral radius (Stability)
println("\nChecking the spectral radius (Stability):")
@show spectralRadiusOfMapping(mapping_LR);

# Calculating fix points (Stationary solution)
@time fp_LR = fixPointOfMapping(mapping_LR);# just a point series representing the full priodic solution



#Plotting one segment based on the Plot recipe for the PeriodicSolution structure
plot( sol_periodic, vars=1, label="SDM Periodic (order 1)", linewidth=2)

display(p4)
