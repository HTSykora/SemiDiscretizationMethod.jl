# mu_convergence_test
###### Multiplication-Free Semi-Discretization Method for Time-Periodic Delayed Systems
###### The mapping is repsented by a left and a right matrix


using SemiDiscretizationMethod
using StaticArrays
using Plots

using CSV
using DataFrames
using BenchmarkTools

using LaTeXStrings
gr()

function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    #τ1=t->1+0.3*sin(t/T*2*pi*3)#2.75π # if function is needed, the use τ1 = t->foo(t)
    #τ1 = 0.5π # if function is needed, the use τ1 = t->foo(t])
    #τ1 =τ1 = t->0.5π-0.5π*cos(2π / T * t * 2) 
    τ1 = τ1 = t -> 1π - 0.5π * sin((t / T))
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    #τ2=1.0π # if function is needed, the use τ1 = t->foo(t)
    τ2 = τ1 = t -> 2.0π + 0.0π * (t / T)
    #TODO: nem ugyan az ha függvény vagy, ha konstans!?!?!?!?! 
    BMx2 = DelayMX(τ2, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, 1.0 * cos(2π / T * t * 4)])
    LDDEProblem(AMx, [BMx1, BMx2], cVec)
    # LDDEProblem(AMx, BMx1, cVec)
end


#----------------------------------------------------------
#-------------------- mu convergence test -----------------
#----------------------------------------------------------
Ndisc = 300_000
NSD_order = 3
τmax = 2π * 1.05 # the largest τ of the system
T = 2π #Principle period of the system (sin(t)=cos(t+T)) 
mathieu_lddep = createMathieuProblem(3.0, 2.9, -0.45, 0.2, T=T)
#mathieu_lddep = createMathieuProblem(3.0, 0, -0.45, 0.2, T=T)
method = SemiDiscretization(NSD_order, T / Ndisc) # 3rd order semi discretization with Δt=0.1


@time mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Ndisc, calculate_additive=true);#The discrete mapping of the system
@show @time μLR_ref = spectralRadiusOfMapping(mappingLR, nev=1, tol=1e-100)


Nv = ceil.(Int, 10 .^ (1.0:0.1:5))

NSD_orderv = 0:2
tt = zeros(size(NSD_orderv, 1), size(Nv, 1))
mui = zeros(size(NSD_orderv, 1), size(Nv, 1))


p1 = scatter(yaxis=:log10, xaxis=:log10, xlabel=L"CPU time [s]", ylabel=L"\mu_{\mathrm{error}}")
p2 = scatter(yaxis=:log10, xaxis=:log10, xlabel=L"number of steps, (N)", ylabel=L"\mu_{\mathrm{error}}")
p3 = scatter(yaxis=:log10, xaxis=:log10, xlabel=L"number of steps, (N)", ylabel=L"CPU time [s]")

@time for (iord, NSD_order) in enumerate(NSD_orderv)
    for kNdisc in vcat([1, 1], 1:length(Nv)) #the first is repated to get read of the first compliation time
        Ndisc = Nv[kNdisc]
        println([NSD_order, Ndisc])

        method = SemiDiscretization(NSD_order, T / Ndisc) # 3rd order semi discretization with Δt=0.1

        tt[iord, kNdisc] = @elapsed mui[iord, kNdisc] = spectralRadiusOfMapping(
            DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Ndisc, calculate_additive=true),
            nev=1, tol=1e-100)
    end

    err = abs.(mui[iord, :] .- μLR_ref)

    # Estimate convergence rate avoiding round-off plateau
    valid_idx = findall((err .> 1e-11) .& (err .< 1e-3))
    if length(valid_idx) >= 2
        x_fit = log10.(Nv[valid_idx])
        y_fit = log10.(err[valid_idx])
        A = hcat(x_fit, ones(length(x_fit)))
        slope, intercept = A \ y_fit
        rate = round(-slope, digits=2)
        label_p2 = "Order $NSD_order (rate: $rate)"
    else
        label_p2 = "Order $NSD_order"
    end

    scatter!(p1, tt[iord, :], err, labels="Order $NSD_order", m=:diamond)
    scatter!(p2, Nv, err, labels=label_p2, m=:utriangle)
    scatter!(p3, Nv, tt[iord, :], labels="Order $NSD_order", m=:xcross)

    display(plot(p1, p2, p3, layout=(1, 3), size=(1200, 400), bottom_margin=5Plots.mm, left_margin=5Plots.mm))
end

