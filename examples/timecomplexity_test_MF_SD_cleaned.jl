using SemiDiscretizationMethod
using StaticArrays
using Plots
using CSV
using DataFrames
using BenchmarkTools
using LaTeXStrings

"""
    createMathieuProblem(δ, ε, b0, a1; T=2π)

Creates a Mathieu LDDE problem for testing.
"""
function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    τ1 = t -> 2π
    BMx1 = DelayMX(τ1, @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, 1.0 * cos(2π / T * t * 4)])
    return LDDEProblem(AMx, BMx1, cVec)
end

"""
    power_fit(x::Vector{Float64}, y::Vector{Float64}, N_threshold)

Fits a power law y = c * x^p to the data for x > N_threshold.
"""
function power_fit(x::Vector{Float64}, y::Vector{Float64}, N_threshold)
    mask = (x .> N_threshold) .& (.!isnan.(y)) .& (y .> 0)
    if sum(mask) < 2
        return Float64[], Float64[]
    end
    xnew = x[mask]
    ynew = y[mask]

    logx = log.(xnew)
    logy = log.(ynew)

    X = hcat(ones(length(logx)), logx)
    β = X \ logy

    c = exp(β[1])
    p = β[2]

    fitted_y = c .* xnew .^ p
    println("Fitted power law: c * x^$p")
    return xnew, fitted_y
end

function run_complexity_test(; Nv=ceil.(10 .^ (1.0:0.1:4.0)), Twaitfor_SH=2.0)
    println("Starting complexity test...")

    # Preallocate result arrays
    n_v = length(Nv)

    # Short-Hand (Traditional) method times
    tmake_SH_PhiALL = fill(NaN, n_v)
    teig_SH = fill(NaN, n_v)
    fixSH = fill(NaN, n_v)

    # Multiplication-Free (LR) method times
    tmake_LR = fill(NaN, n_v)
    teig_LR = fill(NaN, n_v)
    fixLR = fill(NaN, n_v)

    # Standard deviations for ribbon plots
    tmake_SH_PhiALL_S = fill(0.0, n_v)
    teig_SH_S = fill(0.0, n_v)
    fixSH_S = fill(0.0, n_v)
    tmake_LR_S = fill(0.0, n_v)
    teig_LR_S = fill(0.0, n_v)
    fixLR_S = fill(0.0, n_v)

    domoreSH = true
    NEV = 1
    kpow = 1e-4 # tolerance for eigs

    # Set benchmark parameters
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1

    τmax = 2π
    T = 2π

    for (i, Ndisc) in enumerate(Nv)
        println("Testing Ndisc = $Ndisc ($(round(100*i/n_v, digits=1))%)")

        # Recreate problem and method for current resolution
        mathieu_lddep = createMathieuProblem(3.0, 3.0, -0.5, 0.2, T=T)
        method = SemiDiscretization(1, T / Ndisc)
        n_steps = Int((T + 100eps(T)) ÷ method.Δt)

        # 1. Multiplication-Free (LR) method
        # Mapping creation
        b_LR = @benchmark DiscreteMapping_LR($mathieu_lddep, $method, $τmax, n_steps=$n_steps, calculate_additive=true)
        t_LR = BenchmarkTools.median(b_LR)
        tmake_LR[i] = t_LR.time / 1e9
        tmake_LR_S[i] = BenchmarkTools.std(b_LR).time / 1e9

        mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=n_steps, calculate_additive=true)

        # Eigenvalue calculation
        b_eig_LR = @benchmark spectralRadiusOfMapping($mappingLR, nev=$NEV, tol=$kpow)
        t_eig_LR = BenchmarkTools.median(b_eig_LR)
        teig_LR[i] = t_eig_LR.time / 1e9
        teig_LR_S[i] = BenchmarkTools.std(b_eig_LR).time / 1e9

        # Fixed point calculation
        b_fix_LR = @benchmark fixPointOfMapping($mappingLR)
        t_fix_LR = BenchmarkTools.median(b_fix_LR)
        fixLR[i] = t_fix_LR.time / 1e9
        fixLR_S[i] = BenchmarkTools.std(b_fix_LR).time / 1e9

        # 2. Short-Hand (Traditional) method - only if it's not too slow
        if domoreSH
            # Mapping creation (includes matrix products)
            b_SH = @benchmark DiscreteMapping_1step($mathieu_lddep, $method, $τmax, n_steps=$n_steps, calculate_additive=true)
            t_SH = BenchmarkTools.median(b_SH)
            tmake_SH_PhiALL[i] = t_SH.time / 1e9
            tmake_SH_PhiALL_S[i] = BenchmarkTools.std(b_SH).time / 1e9

            mappingFull = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=n_steps, calculate_additive=true)

            # Eigenvalue calculation
            b_eig_SH = @benchmark spectralRadiusOfMapping($mappingFull, nev=$NEV, tol=$kpow)
            t_eig_SH = BenchmarkTools.median(b_eig_SH)
            teig_SH[i] = t_eig_SH.time / 1e9
            teig_SH_S[i] = BenchmarkTools.std(b_eig_SH).time / 1e9

            # Fixed point calculation
            b_fix_SH = @benchmark fixPointOfMapping($mappingFull)
            t_fix_SH = BenchmarkTools.median(b_fix_SH)
            fixSH[i] = t_fix_SH.time / 1e9
            fixSH_S[i] = BenchmarkTools.std(b_fix_SH).time / 1e9

            # Stop traditional method if mapping creation takes too long
            if tmake_SH_PhiALL[i] > Twaitfor_SH
                domoreSH = false
                println("Traditional method is too slow, skipping further resolutions for it.")
            end
        end
        p = plot_results((Nv=Nv, tmake_SH=tmake_SH_PhiALL, teig_SH=teig_SH, fixSH=fixSH,
        tmake_LR=tmake_LR, teig_LR=teig_LR, fixLR=fixLR,
        tmake_SH_S=tmake_SH_PhiALL_S, teig_SH_S=teig_SH_S, fixSH_S=fixSH_S,
        tmake_LR_S=tmake_LR_S, teig_LR_S=teig_LR_S, fixLR_S=fixLR_S))
        display(p)
    end

    return (Nv=Nv, tmake_SH=tmake_SH_PhiALL, teig_SH=teig_SH, fixSH=fixSH,
        tmake_LR=tmake_LR, teig_LR=teig_LR, fixLR=fixLR,
        tmake_SH_S=tmake_SH_PhiALL_S, teig_SH_S=teig_SH_S, fixSH_S=fixSH_S,
        tmake_LR_S=tmake_LR_S, teig_LR_S=teig_LR_S, fixLR_S=fixLR_S)
end

function plot_results(results)
    Nv = results.Nv

    # Scale ribbons for visibility
    ScaleSTD = 0.5

    p = plot(xaxis=:log10, yaxis=:log10, gridlinewidth=2, legend=:bottomright)

    # Plot data with ribbons
    plot!(p, Nv, results.tmake_SH, ribbon=results.tmake_SH_S .* ScaleSTD, label="Traditional: Mapping Creation (prodl)")
    plot!(p, Nv, results.teig_SH, ribbon=results.teig_SH_S .* ScaleSTD, label="Traditional: eigs(Φ)")
    plot!(p, Nv, results.fixSH, ribbon=results.fixSH_S .* ScaleSTD, label="Traditional: fixPoint")

    plot!(p, Nv, results.tmake_LR, ribbon=results.tmake_LR_S .* ScaleSTD, label="MFSD (LR): Mapping Creation", linewidth=2)
    plot!(p, Nv, results.teig_LR, ribbon=results.teig_LR_S .* ScaleSTD, label="MFSD (LR): eigs(ΦR, ΦL)", linewidth=2)
    plot!(p, Nv, results.fixLR, ribbon=results.fixLR_S .* ScaleSTD, label="MFSD (LR): fixPoint", linewidth=2)

    # Add power law fits (dashed lines)
    for (y, label) in [(results.tmake_SH, ""), (results.teig_SH, ""), (results.tmake_LR, ""), (results.teig_LR, "")]
        xf, yf = power_fit(Nv, y, 1e2)
        if !isempty(xf)
            plot!(p, xf, yf, label=label, linewidth=1, linecolor=:black, linestyle=:dash)
        end
    end

    xlabel!(p, L"number of steps, (N)")
    ylabel!(p, "CPU-time [s]")
    title!(p, "Time Complexity Analysis")

    return p
end

# Main execution
results = run_complexity_test(Nv=ceil.(10 .^ (1.0:0.25:5.0)), Twaitfor_SH=1.0)
pfinal = plot_results(results)
display(pfinal)

# # Save results to CSV
# df = DataFrame(results)
# CSV.write("TimeComplexity_Results.csv", df)
# println("Results saved to TimeComplexity_Results.csv")
