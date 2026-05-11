# sdm_periodic_sol.jl
# This file provides a DifferentialEquations.jl-like interface for periodic solutions 
# extracted from SemiDiscretizationMethod.jl mappings.

using LinearAlgebra
using RecipesBase
using StaticArrays


"""
    PeriodicSolution{d, tT, uT}

A structure representing a periodic solution extracted from a Semi-Discretization mapping.
It behaves like a function `sol(t)` and supports interpolation based on the method order used during discretization.

# Fields
- `t::Vector{tT}`: Time grid points within one period [0, Δt, ..., T].
- `u::Vector{uT}`: State vectors (typically `SVector`) at each time point.
- `T::tT`: The period length.
- `order::Int`: The interpolation order (0: Constant/ZOH, 1: Linear, etc.).

# Usage
The object is callable: `sol(t)` returns the interpolated state at time `t`. 
It automatically handles periodicity using `mod(t, T)`.
"""
struct PeriodicSolution{d, tT, uT}
    t::Vector{tT}      # Time points [0, Δt, ..., T]
    u::Vector{uT}      # State vectors at time points
    T::tT              # Period length
    order::Int         # Interpolation order
end

# --- Make the object callable like a function: sol(t) ---
function (sol::PeriodicSolution{d})(t_call) where d
    # Use modulo to map any time point into the principal period [0, T]
    t_loc = mod(t_call, sol.T)
    
    # Handle boundary conditions for safety
    if t_loc <= sol.t[1]
        return sol.u[1]
    elseif t_loc >= sol.t[end]
        return sol.u[end]
    end
    
    # Calculate index based on uniform grid (SDM default)
    dt = sol.t[2] - sol.t[1]
    idx = floor(Int, t_loc / dt) + 1
    idx = clamp(idx, 1, length(sol.t) - 1)
    
    if sol.order == 0
        return sol.u[idx]
    else
        # Higher order Lagrange interpolation matching SDM logic but as interpolation
        # To hit both u[idx] and u[idx+1], we use points u[idx+1], u[idx], u[idx-1], ...
        # and evaluate relative to the right endpoint t[idx+1].
        
        n = length(sol.u) - 1 # number of intervals
        res = zero(sol.u[1])
        
        # Local time relative to the RIGHT end of the interval
        θ = t_loc - sol.t[idx+1]
        # Effective step is negative because we go backwards from idx+1
        h_eff = -dt
        
        for k in 0:sol.order
            # Using the lagr_el0 from functions_utility.jl
            # lagr_el0(q, k, τerr, h, t)
            # Here point k=0 is u[idx+1] at θ=0, k=1 is u[idx] at θ=-dt, etc.
            L = lagr_el0(sol.order, k, 0.0, h_eff, θ)
            
            # Use mod1 for periodic indexing of u (1 to n)
            u_idx = mod1(idx + 1 - k, n)
            res += sol.u[u_idx] * L
        end
        return res
    end
end

"""
    get_periodic_solution(mappLR::DiscreteMapping_LR, d::Int, order::Int=1)

Core function to extract the periodic solution from a Semi-Discretization mapping.

# Arguments
- `mappLR`: The discrete mapping (must be calculated with `calculate_additive=true`).
- `d::Int`: The dimension of the state space.
- `order::Int`: The interpolation order to be used in the resulting object.

# Returns
- A `PeriodicSolution` object.
"""
function get_periodic_solution(mappLR, d::Int, order::Int=1)
    # Calculate the raw fixed point vector y* = (L - R) \ v
    y_raw = (mappLR.LmappingMX - mappLR.RmappingMX) \ Vector(mappLR.mappingVs[1])
    
    # In LR mapping, the number of steps is length(y_raw) / d
    p_steps = div(length(y_raw), d)
    T_period = mappLR.ts[end] - mappLR.ts[1]
    
    # We reconstruct the solution in forward time order: [0, Δt, ..., T]
    t_vals = collect(range(0, T_period, length=p_steps + 1))
    
    # Slice the raw vector into state vectors (backward order from mapping)
    u_blocks_backward = [y_raw[(1 + (i-1)*d):(i*d)] for i in 1:p_steps]
    u_blocks_forward = reverse(u_blocks_backward)
    
    # Complete the [0, T] range: u[1] = x(0) = x(T)
    # Use vcat to avoid flattening if the state is a vector
    u_vals = vcat([u_blocks_forward[end]], u_blocks_forward)
    
    return PeriodicSolution{d, typeof(T_period), eltype(u_vals)}(t_vals, u_vals, T_period, order)
end

# --- Convenience overloads ---

"""
    get_periodic_solution(mappLR, prob::LDDEProblem{d}, [method::DiscretizationMethod])

Convenience wrapper that automatically infers the dimension `d` and interpolation order.

# Example
```julia
sol = get_periodic_solution(mapping_LR, prob, method)
plot(sol, vars=1) # Smooth plot using the correct interpolation
```
"""
function get_periodic_solution(mappLR, prob::LDDEProblem{d}, method::DiscretizationMethod) where d
    return get_periodic_solution(mappLR, d, methodorder(method))
end

function get_periodic_solution(mappLR, prob::LDDEProblem{d}, order::Int=1) where d
    return get_periodic_solution(mappLR, d, order)
end

function get_periodic_solution(mappLR, d::Int, method::DiscretizationMethod)
    return get_periodic_solution(mappLR, d, methodorder(method))
end

# --- Plots.jl Recipe ---
@recipe function f(sol::PeriodicSolution; vars=1, n_plot=nothing)
    # Generate points for a smooth plot
    if n_plot === nothing
        n_plot = max(300, length(sol.t) * 10)
    end
    t_plot = range(0, sol.T, length=n_plot)
    
    # Determine which variables to plot
    if vars isa AbstractVector
        for v in vars
            @series begin
                label --> "x_$v(t)"
                t_plot, [sol(t)[v] for t in t_plot]
            end
        end
    else
        u_plot = [sol(t)[vars] for t in t_plot]
        xlabel --> "Time (t)"
        ylabel --> "x_$(vars)(t)"
        label  --> "Periodic Solution (SDM)"
        t_plot, u_plot
    end
end