# mu_convergence_test
######Product-Free full mapping
######The mapping is repsented by a left and a right matrix
5 + 5

###import Pkg
########
###Pkg.activate("")
######### Pkg.activate()# ez az eredeti, mindent tartalamzo
#########] dev SemiDiscretizationMethod

using Revise
using SemiDiscretizationMethod
using StaticArrays
using Plots

using CSV
using DataFrames
using BenchmarkTools

using LaTeXStrings
gr()
Threads.nthreads()

function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    #τ1=t->1+0.3*sin(t/T*2*pi*3)#2.75π # if function is needed, the use τ1 = t->foo(t)
    #τ1 = 0.5π # if function is needed, the use τ1 = t->foo(t])
    #τ1 =τ1 = t->0.5π-0.5π*cos(2π / T * t * 2) 
    τ1 = τ1 = t -> 1π - 0.5π * sin((t / T))
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    #τ2=1.0π # if function is needed, the use τ1 = t->foo(t)
    τ2 = τ1 = t -> 1.0π + 0.0π * (t / T)
    #TODO: nem ugyan az ha függvény vagy, ha konstans!?!?!?!?! 
    BMx2 = DelayMX(τ2, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, 1.0 * cos(2π / T * t * 4)])
    LDDEProblem(AMx, [BMx1, BMx2], cVec)
    # LDDEProblem(AMx, BMx1, cVec)
end


#----------------------------------------------------------
#-------------------- mu convergence test -----------------
#----------------------------------------------------------
Ndisc = Int(1e5)
NSD_order = 5
τmax = 2π # the largest τ of the system
T = 2π #Principle period of the system (sin(t)=cos(t+T)) 
mathieu_lddep = createMathieuProblem(3.0, 2.9, -0.45, 0.2, T=T)
method = SemiDiscretization(NSD_order, T / Ndisc) # 3rd order semi discretization with Δt=0.1


@time mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Ndisc, calculate_additive=true);#The discrete mapping of the system
@show @time μLR_ref = spectralRadiusOfMapping(mappingLR,nev=1, tol=1e-100)


Nv = ceil.(Int, 10 .^ (1.5:0.1:4.00))
NSD_orderv = 0:4
tt = zeros(size(NSD_orderv, 1), size(Nv, 1))
mui = zeros(size(NSD_orderv, 1), size(Nv, 1))


scatter()
for (iord, NSD_order) in enumerate(NSD_orderv)
a


    for kNdisc in vcat([1, 1], 1:length(Nv)) #the first is repated to get read of the first compliation time
        Ndisc = Nv[kNdisc]
        println([NSD_order, Ndisc])

        method = SemiDiscretization(NSD_order, T / Ndisc) # 3rd order semi discretization with Δt=0.1

        tt[iord, kNdisc] += @elapsed mui[iord, kNdisc] = spectralRadiusOfMapping(
            DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Ndisc, calculate_additive=true),
            nev=1, tol=1e-100)
    end
    #display(
    #    scatter!(tt[iord, :], (abs.(mui[iord, :] .- μLR_ref)),yaxis=:log10,xaxis=:log10, labels=NSD_order,
    #    xlabel=L"number of steps, (N)", ylabel=L"\mu error")
    #)
    display(
        scatter!(Nv, (abs.(mui[iord, :] .- μLR_ref)),yaxis=:log10,xaxis=:log10, labels=NSD_order,
        xlabel=L"CPU time [s]", ylabel=L"\mu _{error}")
    )
end