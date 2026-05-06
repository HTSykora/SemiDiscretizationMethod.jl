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
T_per = 2π * 5.5
τmax = T_per * 0.5# 2π
params = (13.0, 0.2, -0.0015, 0.1, T_per, τmax)

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
println("Calculating mapping matrices...")
@time mapping = DiscreteMapping(mathieu_lddep, method, τmax,
    n_steps=n_steps_calc, calculate_additive=true);

@time mapping_LR = DiscreteMapping_LR(mathieu_lddep, method, τmax,
    n_steps=n_steps_calc, calculate_additive=true);

# Checking the spectral radius (Stability)
println("\nChecking the spectral radius (Stability):")
@show spectralRadiusOfMapping(mapping);
@show spectralRadiusOfMapping(mapping_LR);

# Calculating fix points (Stationary solution)
println("\nCalculating fix points (Stationary solution)...")
@time fp = fixPointOfMapping(mapping);# just a point series representing the (statespace) history of the starting point of the priodic solution
@time fp_LR = fixPointOfMapping(mapping_LR);# just a point series representing the full priodic solution

# =======================================================================
# 3. PLOTTING THE SDM RESULTS based on the raw datapoint of the periodic solution
# =======================================================================
println("Generating SDM plots...")

r_len = length(fp[1:2:end])
t_forward = (0:-1:-(r_len - 1)) .* method.Δt
x_fp_forward = fp[1:2:end]
dx_fp_forward = fp[2:2:end]

p1 = plot(t_forward, dx_fp_forward, xlabel="t", label=L"\dot{x}(t)", linewidth=4)
plot!(p1, t_forward, x_fp_forward, xlabel="t", label=L"x(t)", linewidth=4)

p_len = length(fp_LR[1:2:end])
t_LR_forward = (0:-1:-(p_len - 1)) .* method.Δt
x_LR_forward = fp_LR[1:2:end]
dx_LR_forward = fp_LR[2:2:end]

plot!(p1, t_LR_forward, dx_LR_forward, label=L"\dot{x}_{LR}(t)", linewidth=2)
plot!(p1, t_LR_forward, x_LR_forward,
    xlabel=L"t",
    title=L"Periodic solution \quad t \in [0,T]",
    guidefontsize=14,
    linewidth=2,
    label=L"x_{LR}(t)",
    legendfontsize=11,
    tickfont=font(10)
)

# Plotting the excitation for reference
tloc = 0.0:method.Δt:T_per
s = [t < T_per / 5 ? 1 : 0 for t in tloc]
plot!(p1, tloc .- T_per, s, label="Excitation s(t)")

display(p1)


# =======================================================================
# 4. LONG-TIME SIMULATION & COMPARISON WITH DIFFERENTIALEQUATIONS.JL
# compared to the interpolatable extracted periodic solution structure
# =======================================================================
println("\nRunning long-time simulation with DifferentialEquations.jl...")

# Define history function for the DDE solver (constant 0 for t < 0)
h_dde(p_arg, t_req; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : @SVector [0.0, 0.0]

# Time span: simulate for 200 periods to let transients die out completely
periods = 10
tspan = (0.0, periods * T_per)

# Initial condition (arbitrary non-zero state)
u0 = @SVector [1.0, 0.0]

# Define and solve the DDE problem
dde_prob = DDEProblem(mathieu_dde_rhs, u0, h_dde, tspan, params; constant_lags=[τmax])
sol = solve(dde_prob, reltol=1e-8, abstol=1e-8)
# Generate comparison plot
println("Generating comparison plot...")
p2 = plot(sol)

plot!(p2, t_LR_forward .+ tspan[end], dx_LR_forward, label=L"\dot{x}_{LR}(t)", linewidth=2)

display(p2)
println("Done! The extracted system works perfectly with the SemiDiscretizationMethod and matches the long-time DDE simulation.")


# --- Plotting the periodic solution vs transient ---
println("\nPlotting transient settling on periodic orbit...")
t_full = 0.0:0.02:tspan[end]
sol_periodic = get_periodic_solution(mapping_LR,2)# the dimension of the system must be provided
#alternativa function call:
#sol_periodic = get_periodic_solution(mapping_LR, mathieu_lddep, 1)# the last element is the order of interpolation, which can be inferred from the method as well:
long_periodic_sol = sol_periodic.(t_full) # automatcially uses interpolation to get values at any t, not just the discrete points, handling the periodicity as well

p4 = plot(sol, vars=(0, 1), label="DDE Transient (x)", alpha=0.6)
plot!(p4, t_full, getindex.(long_periodic_sol, 1),
    label="SDM Periodic Limit Cycle", linewidth=2, linestyle=:dash)
xlabel!(p4, "Time (t)")
ylabel!(p4, "x(t)")


#Plotting one segment based on the Plot recipe for the PeriodicSolution structure
plot!(p4, sol_periodic, vars=1, label="SDM Periodic (order 1)", linewidth=2)


display(p4)
