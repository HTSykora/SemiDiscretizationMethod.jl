######Product-Free full mapping
######The mapping is repsented by a left and a right matrix
## Time-complexity measurements for a publication

5+5
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
Threads.nthreads()


function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    #τ1=t->1+0.3*sin(t/T*2*pi*3)#2.75π # if function is needed, the use τ1 = t->foo(t)
    #τ1 = 0.5π # if function is needed, the use τ1 = t->foo(t])
    #τ1 =τ1 = t->0.5π-0.5π*cos(2π / T * t * 2) 
    #τ1 = τ1 = t -> 2π - 0.5π * (t / T)
    τ1 = τ1 = t -> 2π
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    #τ2=1.0π # if function is needed, the use τ1 = t->foo(t)
    #τ2 = τ1 = t -> 1.0π + 1.0π * (t / T)
    ##TODO: nem ugyan az ha függvény vagy, ha konstans!?!?!?!?! 
    #BMx2 = DelayMX(τ2, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, 1.0 * cos(2π / T * t * 4)])
    #LDDEProblem(AMx, [BMx1, BMx2], cVec)
    LDDEProblem(AMx, BMx1, cVec)
end

Ndisc = 200

τmax = 2π # the largest τ of the system
T = 20π #Principle period of the system (sin(t)=cos(t+T)) 
mathieu_lddep = createMathieuProblem(3.0, 3.0, -0.5, 0.2, T=T) # LDDE problem for Hayes equation
method = SemiDiscretization(5, T / Ndisc) # 3rd order semi discretization with Δt=0.1
Nsteps = Int((T + 100eps(T)) ÷ method.Δt)

#mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)#The discrete mapping of the system
#μLR = spectralRadiusOfMapping(mappingLR)

@time mapping = DiscreteMapping(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
@show @time μ = spectralRadiusOfMapping(mapping)
@time mapping = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
@show @time μ = spectralRadiusOfMapping(mapping)
@time mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
@show @time μLR = spectralRadiusOfMapping(mappingLR)



BenchmarkTools.DEFAULT_PARAMETERS.samples = 100.0
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3.0
@show t = @benchmark spectralRadiusOfMapping($mappingLR)
BenchmarkTools.median(t).time / 1e9


@benchmark abs(eigs(mappingLR.RmappingMX, mappingLR.LmappingMX)[1][1])
@benchmark abs(eigen(collect(mappingLR.RmappingMX),collect(mappingLR.LmappingMX),sortby=abs).values[end])

#tiem of eig test for different precision-------
tt = []
muerror = []
kpowv = -15:0.25:0
for kpow in kpowv
    tloc = @elapsed μLR = spectralRadiusOfMapping(mapping, nev=1, tol=10.0^kpow)
    push!(tt, tloc)
    push!(muerror, μ - μLR)
    println(μ - μLR)
end
scatter(tt, log.(abs.(muerror)))
#scatter!(tt,log.(abs.(muerror)))
#scatter(kpowv,log.(abs.(muerror)))
#scatter!(tt,kpowv)

#--------- profiling the code ---------------

using Profile
function foo(n)
    for k in 1:n
        map = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)#The discrete mapping of the system
        μLR_SP = spectralRadiusOfMapping(map, nev=1, tol=1e-5)
    end
end
foo(4)

@profview foo(100)







## --------------- time test -------------------------
## -----------------------------------------------------

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000.0
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2.0


Nv = ceil.(10 .^ (1.0:0.01:5.61)) #5.61  #100:100:3000
Nv = ceil.(10 .^ (1.0:0.05:6))
Twaitfor_SH = 10.0;


BenchmarkTools.DEFAULT_PARAMETERS.samples = 50.0
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1
Nv = ceil.(10 .^ (1.0:0.025:3.00))
Twaitfor_SH = 1.0;

kpow = 1e-4
NEV = 1
#kpow=1e-20  #has stron effect in T<tau
#NEV=20

tmake_SH_PRi = zeros(Float64, length(Nv));
tmake_SH_Ci = zeros(Float64, length(Nv));
tmake_SH_PhiALL = zeros(Float64, length(Nv));
teig_SH = zeros(Float64, length(Nv));
tfixP_SH = zeros(Float64, length(Nv));
μSH = zeros(Float64, length(Nv));
fixSH = zeros(Float64, length(Nv));

#deviation
tmake_SH_PRi_S = zeros(Float64, length(Nv));
tmake_SH_Ci_S = zeros(Float64, length(Nv));
tmake_SH_PhiALL_S = zeros(Float64, length(Nv));
teig_SH_S = zeros(Float64, length(Nv));
tfixP_SH_S = zeros(Float64, length(Nv));
fixSH_S = zeros(Float64, length(Nv));


domoreSH = true;
domoreSH_PR = true;
#@show domoreSH=false;

tmake_LR = zeros(Float64, length(Nv));
teig_LR = zeros(Float64, length(Nv));
tfixP_LR = zeros(Float64, length(Nv));
μLR = zeros(Float64, length(Nv));
fixLR= zeros(Float64, length(Nv));
#deviation
tmake_LR_S = zeros(Float64, length(Nv));
teig_LR_S = zeros(Float64, length(Nv));
tfixP_LR_S = zeros(Float64, length(Nv));
fixLR_S= zeros(Float64, length(Nv));


f1 = (x) -> log(x) ./ log(10)
f2 = (x) -> log(x) ./ log(10)

#kNdisc=10
#@profview  
for kNdisc in vcat([1, 1], 1:length(Nv)) #the first is repated to get read of the first compliation time
    if kNdisc == 1
        tmake_SH_PRi[kNdisc] = 0.0
        tmake_SH_Ci[kNdisc] = 0.0
        tmake_SH_PhiALL[kNdisc] = 0.0
        teig_SH[kNdisc] = 0.0
        tfixP_SH[kNdisc] = 0.0
        μSH[kNdisc] = 0.0
        fixSH[kNdisc] = 0.0

        domoreSH = true
        domoreSH_PR = true
        tmake_LR[kNdisc] = 0.0
        teig_LR[kNdisc] = 0.0
        tfixP_LR[kNdisc] = 0.0
        μLR[kNdisc] = 0.0
        fixLR[kNdisc] = 0.0
    end
    Ndisc = Nv[kNdisc]
    println([kNdisc / length(Nv), Ndisc])

    #@show Naver = maximum([ceil(4 - (log(Ndisc) / log(10))) * 2, 1])
    @show Naver = 1
    τmax = 2π # the largest τ of the system
    T = 20π #Principle period of the system (sin(t)=cos(t+T)) 
    mathieu_lddep = createMathieuProblem(3.0, 3.0, -0.5, 0.2, T=T) # LDDE problem for Hayes equation

    method = SemiDiscretization(1, T / Ndisc) # 3rd order semi discretization with Δt=0.1

    Nsteps = Int((T + 100eps(T)) ÷ method.Δt)
    for _ = 1:Naver


        #tmake_SH_PRi[kNdisc] += @elapsed resultPR = SemiDiscretizationMethod.calculateResults(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)
        t = @benchmark SemiDiscretizationMethod.calculateResults($mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)

        tmake_SH_PRi[kNdisc] += BenchmarkTools.median(t).time / 1e9
        tmake_SH_PRi_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
        resultPR = SemiDiscretizationMethod.calculateResults(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)

        if domoreSH_PR
            #tmake_SH_Ci[kNdisc] += @elapsed mmpp2 = SemiDiscretizationMethod.DiscreteMappingSteps(resultPR)
            t = @benchmark SemiDiscretizationMethod.DiscreteMappingSteps($resultPR)
            @show t
            tmake_SH_Ci[kNdisc] += BenchmarkTools.median(t).time / 1e9
            tmake_SH_Ci_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
            mmpp2 = SemiDiscretizationMethod.DiscreteMappingSteps(resultPR)


            domoreSH_PR = tmake_SH_Ci[kNdisc] < Twaitfor_SH * Naver
        else
            tmake_SH_Ci[kNdisc] = NaN
        end
        if domoreSH
            println("domoreSH")
            #tmake_SH_PhiALL[kNdisc] += @elapsed mappingFull = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system
            t = @benchmark DiscreteMapping_1step($mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system
            @show t
            tmake_SH_PhiALL[kNdisc] += BenchmarkTools.median(t).time / 1e9
            tmake_SH_PhiALL_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
            mappingFull = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system


            #teig_SH[kNdisc] += @elapsed μSH[kNdisc] = spectralRadiusOfMapping(mappingFull,nev=NEV,tol=10.0^kpow)
            t = @benchmark spectralRadiusOfMapping($mappingFull, nev=NEV, tol=10.0^kpow)
            @show t
            teig_SH[kNdisc] += BenchmarkTools.median(t).time / 1e9
            teig_SH_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
            μSH[kNdisc] = spectralRadiusOfMapping(mappingFull, nev=NEV, tol=10.0^kpow)

            # fixSH[kNdisc]  += @elapsed μSH[kNdisc] = fixPointOfMapping(mappingFull)
            t = @benchmark fixPointOfMapping($mappingFull)
            @show t
            fixSH[kNdisc] += BenchmarkTools.median(t).time / 1e9
            fixSH_S[kNdisc] += BenchmarkTools.std(t).time / 1e9


            domoreSH = tmake_SH_PhiALL[kNdisc] < Twaitfor_SH * Naver
        else
            tmake_SH_PhiALL[kNdisc] = NaN
            teig_SH[kNdisc] = NaN
            fixSH[kNdisc] = NaN
        end

        #tmake_LR[kNdisc] += @elapsed mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)#The discrete mapping of the system
        t = @benchmark DiscreteMapping_LR($mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)#The discrete mapping of the system
        @show t
        tmake_LR[kNdisc] += BenchmarkTools.median(t).time / 1e9
        tmake_LR_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
        mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)#The discrete mapping of the system



        #teig_LR[kNdisc] += @elapsed μLR[kNdisc] = spectralRadiusOfMapping(mappingLR,nev=NEV,tol=10.0^kpow)
        t = @benchmark spectralRadiusOfMapping($mappingLR, nev=NEV, tol=10.0^kpow)
        @show t
        teig_LR[kNdisc] += BenchmarkTools.median(t).time / 1e9
        teig_LR_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
        μLR[kNdisc] = spectralRadiusOfMapping(mappingLR, nev=NEV, tol=10.0^kpow)

        # tfixP_LR[kNdisc] += @elapsed xPLR = fixPointOfMapping(mappingLR)
        t = @benchmark fixPointOfMapping($mappingLR)
        @show t
        fixLR[kNdisc] += BenchmarkTools.median(t).time / 1e9
        fixLR_S[kNdisc] += BenchmarkTools.std(t).time / 1e9

    end

    tmake_SH_PRi[kNdisc] /= Naver
    tmake_SH_Ci[kNdisc] /= Naver
    tmake_SH_PhiALL[kNdisc] /= Naver
    teig_SH[kNdisc] /= Naver
    fixSH[kNdisc] /= Naver

    tmake_LR[kNdisc] /= Naver
    teig_LR[kNdisc] /= Naver
    fixLR[kNdisc] /= Naver


    tmake_SH_PRi_S[kNdisc] /= Naver
    tmake_SH_Ci_S[kNdisc] /= Naver
    tmake_SH_PhiALL_S[kNdisc] /= Naver
    teig_SH_S[kNdisc] /= Naver
    fixSH_S[kNdisc] /= Naver

    tmake_LR_S[kNdisc] /= Naver
    teig_LR_S[kNdisc] /= Naver
    fixLR_S[kNdisc] /= Naver


    #if T=τmax
    #@show norm(xP-xPLR)
    #plot(xP[1:2:end])
    #plot!(xPLR[1:2:end])

    # @show norm(μSH[kNdisc]-μLR[kNdisc])

    #if (kNdisc>5000  ||  mod(kNdisc,10)==0)
    df = DataFrame(
        Nv_data=Nv,
        tmake_SH_data=tmake_SH_PhiALL,
        teig_SH_data=teig_SH,
        tfixP_SH_data=tfixP_SH,
        tmake_LR_data=tmake_LR,
        teig_LR_data=teig_LR,
        tfixP_LR_data=tfixP_LR,
        Nv_data_log=f2.(Nv),
        tmake_SH_data_log=f1.(tmake_SH_PhiALL),
        teig_SH_data_log=f1.(teig_SH),
        tfixP_SH_data_log=f1.(tfixP_SH),
        tmake_LR_data_log=f1.(tmake_LR),
        teig_LR_data_log=f1.(teig_LR),
        tfixP_LR_data_log=f1.(tfixP_LR))

    CSV.write("CPU_Timev__T_" * string(T) * "_tau_" * string(τmax) * ".csv", df)
    println("data saved")
    #end


    

    tmake_SH_PRi_S[isnan.(tmake_SH_PRi_S)] .= 0.0
    tmake_SH_Ci_S[isnan.(tmake_SH_Ci_S)] .= 0.0
    tmake_SH_PhiALL_S[isnan.(tmake_SH_PhiALL_S)] .= 0.0
    teig_SH_S[isnan.(teig_SH_S)] .= 0.0
    fixSH_S[isnan.(fixSH_S)] .= 0.0

    tmake_LR_S[isnan.(tmake_LR_S)] .= 0.0
    teig_LR_S[isnan.(teig_LR_S)] .= 0.0
    fixLR_S[isnan.(fixLR_S)] .= 0.0
    

    #f1 = (x) -> log(x) ./ log(10)
    #f2 = (x) -> log(x) ./ log(10)

f1 = (x) -> x == 0.0 ? NaN : x
f2 = (x) -> x == 0.0 ? NaN : x

ScaleSTD=0.5

#fplot=scatter
#fplot! = scatter!
fplot = plot
fplot! = plot!

#fplot(f2.(Nv), f1.(tmake_SH_PRi .+tmake_SH_PRi_S), labels="tmake_SH_PRi")
#fplot!(f2.(Nv), f1.(tmake_SH_PRi ), labels="tmake_SH_PRi")
#fplot!(f2.(Nv), f1.(tmake_SH_PRi .-tmake_SH_PRi_S), labels="tmake_SH_PRi")
fplot(f2.(Nv), f1.(tmake_SH_PRi), ribbon=tmake_SH_PRi_S .* ScaleSTD, labels="tmake_SH_PRi")
fplot!(f2.(Nv), f1.(tmake_SH_Ci), ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
#fplot!(f2.(Nv), f1.(tmake_SH_Ci),ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
fplot!(f2.(Nv), f1.(tmake_SH_PhiALL), ribbon=tmake_SH_PhiALL_S .* ScaleSTD, labels="tmake_SH_Phi_1step...prod(...)")
fplot!(f2.(Nv), f1.(teig_SH), ribbon=teig_SH_S .* ScaleSTD, labels="teig_SH: eigs(PHI) only")
#fplot!(f2.(Nv), f1.(tfixP_SH),ribbon=tfixP_SH_S .* ScaleSTD, labels="teig_SH: fix Point")
fplot!(f2.(Nv), f1.(tmake_LR), ribbon=tmake_LR_S .* ScaleSTD, labels="tmake_LR", linewidth=2)
fplot!(f2.(Nv), f1.(teig_LR), ribbon=teig_LR_S .* ScaleSTD, labels="teig_LR: eigs(ΦR,ΦL) only", linewidth=2)

fplot!(f2.(Nv), f1.(fixSH), ribbon=fixSH_S .* ScaleSTD, labels="fixSH", linewidth=1)
fplot!(f2.(Nv), f1.(fixLR), ribbon=fixLR_S .* ScaleSTD, labels="fixLR", linewidth=2)

#fplot!(xticks =collect( 10 .^(1:0.25:5))) 
#fplot!(yticks =collect( 10.0.^ (-5:1:5))) 
fplot!(yaxis=(:log10, [0.00001, :auto]))
fplot!(xaxis=:log10)
fplot!(gridlinewidth=2)
display(
    fplot!(labels="teig_LR: fix Point",
        xlabel=L"number of steps, (N)", ylabel=L"CPU-time", legend=:bottomright,
        xticks=10.0 .^ (1:6), yticks=10.0 .^ (-5:3))
)
    
    savefig( "myplot3.svg") 
end


tmake_SH_PRi_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(tmake_SH_PRi, tmake_SH_PRi_S)]
tmake_SH_Ci_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(tmake_SH_Ci, tmake_SH_Ci_S)]
tmake_SH_PhiALL_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(tmake_SH_PhiALL, tmake_SH_PhiALL_S)]
teig_SH_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(teig_SH, teig_SH_S)]
fixSH_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(fixSH, fixSH_S)]

tmake_LR_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(tmake_LR, tmake_LR_S)]
teig_LR_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(teig_LR, teig_LR_S)]
fixLR_S .= [-maximum([-a*0.99, -b]) for (a, b) in zip(fixLR, fixLR_S)]



#tmake_SH_PRi_S[isnan.(tmake_SH_PRi_S)] .= 0.0
#tmake_SH_Ci_S[isnan.(tmake_SH_Ci_S)] .= 0.0
#tmake_SH_PhiALL_S[isnan.(tmake_SH_PhiALL_S)] .= 0.0
#teig_SH_S[isnan.(teig_SH_S)] .= 0.0
#tfixP_SH_S[isnan.(tfixP_SH_S)] .= 0.0
#tmake_LR_S[isnan.(tmake_LR_S)] .= 0.0
#teig_LR_S[isnan.(teig_LR_S)] .= 0.0
#tfixP_LR_S[isnan.(tfixP_LR_S)] .= 0.0


f1 = (x) -> x == 0.0 ? NaN : x
f2 = (x) -> x == 0.0 ? NaN : x

ScaleSTD=0.5

#fplot=scatter
#fplot! = scatter!
fplot = plot
fplot! = plot!

#fplot(f2.(Nv), f1.(tmake_SH_PRi .+tmake_SH_PRi_S), labels="tmake_SH_PRi")
#fplot!(f2.(Nv), f1.(tmake_SH_PRi ), labels="tmake_SH_PRi")
#fplot!(f2.(Nv), f1.(tmake_SH_PRi .-tmake_SH_PRi_S), labels="tmake_SH_PRi")
fplot(f2.(Nv), f1.(tmake_SH_PRi), ribbon=tmake_SH_PRi_S .* ScaleSTD, labels="tmake_SH_PRi")
fplot!(f2.(Nv), f1.(tmake_SH_Ci), ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
#fplot!(f2.(Nv), f1.(tmake_SH_Ci),ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
fplot!(f2.(Nv), f1.(tmake_SH_PhiALL), ribbon=tmake_SH_PhiALL_S .* ScaleSTD, labels="tmake_SH_Phi_1step...prod(...)")
fplot!(f2.(Nv), f1.(teig_SH), ribbon=teig_SH_S .* ScaleSTD, labels="teig_SH: eigs(PHI) only")
#fplot!(f2.(Nv), f1.(tfixP_SH),ribbon=tfixP_SH_S .* ScaleSTD, labels="teig_SH: fix Point")
fplot!(f2.(Nv), f1.(tmake_LR), ribbon=tmake_LR_S .* ScaleSTD, labels="tmake_LR", linewidth=2)
fplot!(f2.(Nv), f1.(teig_LR), ribbon=teig_LR_S .* ScaleSTD, labels="teig_LR: eigs(ΦR,ΦL) only", linewidth=2)


fplot!(f2.(Nv), f1.(fixSH), ribbon=fixSH_S .* ScaleSTD, labels="fixSH", linewidth=1)
fplot!(f2.(Nv), f1.(fixLR), ribbon=fixLR_S .* ScaleSTD, labels="fixLR", linewidth=2)



#fplot!(xticks =collect( 10 .^(1:0.25:5))) 
#fplot!(yticks =collect( 10.0.^ (-5:1:5))) 
fplot!(yaxis=(:log10, [0.00001, :auto]))
fplot!(xaxis=:log10)
fplot!(gridlinewidth=2)



function power_fit(x::Vector{Float64}, y::Vector{Float64},N)
    ynew=deepcopy(y[x .>N])
    xnew=deepcopy(x[x .> N])
    xnew=(xnew[.!isnan.(ynew)])
    ynew=(ynew[.!isnan.(ynew)])
    logx = log.(xnew)
    logy = log.(ynew)
    
    # Linear fit to the transformed problem
    X = hcat(ones(length(logx)), logx)
    β = X \ logy  # Solves for log(c) and p
    
    c = exp(β[1])  # Convert back from log(c) to c
    p = β[2]
    
    # Calculate fitted y values using the original x values
    fitted_y = c .* xnew .^ p
    println([c,p])
    return xnew,fitted_y #c, p, 
end

fplot!(power_fit(Nv, tmake_SH_PRi,1e2)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, tmake_SH_Ci,1e2)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, tmake_SH_PhiALL,1e2)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, teig_SH,700)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)

fplot!(power_fit(Nv, tmake_LR,1e2)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, teig_LR,1e2)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)

fplot!(power_fit(Nv, fixSH,1e2)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, fixLR,1e2)...,  labels="", linewidth=1, linecolor=:black, linestyle=:dash)


display(
    fplot!(labels="teig_LR: fix Point",
        xlabel=L"number of steps, (N)", ylabel=L"CPU-time", legend=:bottomright,
        xticks=10.0 .^ (1:6), yticks=10.0 .^ (-5:3))
)


savefig( "myplot_final32STD05_T20pi_2__eigs_test.svg") 