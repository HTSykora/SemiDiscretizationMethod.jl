######Product-Free full mapping
######The mapping is repsented by a left and a right matrix
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
    τ1 =τ1 = t->1π-0.5π*(  t /  T) 
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    #τ2=1.0π # if function is needed, the use τ1 = t->foo(t)
    τ2 =τ1 = t->2.0π-1.0π*(  t /  T) 
    #TODO: nem ugyan az ha függvény vagy, ha konstans!?!?!?!?! 
    BMx2 = DelayMX(τ2,t->@SMatrix [0. 0.; b0 0.]);
    cVec = Additive(t -> @SVector [0.0, 1.0 * cos(2π / T * t * 4)])
    LDDEProblem(AMx,[BMx1,BMx2],cVec)
   # LDDEProblem(AMx, BMx1, cVec)
end

Ndisc=150

τmax = 2π # the largest τ of the system
T = 0.1π #Principle period of the system (sin(t)=cos(t+T)) 
mathieu_lddep = createMathieuProblem(3.0, 2.9, -0.45, 0.2, T=T) # LDDE problem for Hayes equation
method = SemiDiscretization(5 ,T  / Ndisc) # 3rd order semi discretization with Δt=0.1
Nsteps = Int((T + 100eps(T)) ÷ method.Δt)

#mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)#The discrete mapping of the system
#μLR = spectralRadiusOfMapping(mappingLR)

@time mapping = DiscreteMapping(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
@show @time μ = spectralRadiusOfMapping(mapping)
@time mapping = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
@show @time μ = spectralRadiusOfMapping(mapping)

@time mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
@show @time μLR = spectralRadiusOfMapping(mappingLR_SP)

using Profile
function foo(n)
for k in  1:n
    map = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
    μLR_SP = spectralRadiusOfMapping(map)
end
end
foo(4)

@profview foo(10)

# --------------- time test -------------------------
# -----------------------------------------------------



Nv = ceil.(10 .^ (1.0:0.01:5.61)) #5.61  #100:100:3000
Nv = ceil.(10 .^ (1.0:0.02:6))
Twaitfor_SH = 10.0;

#Nv = ceil.(10 .^ (1.0:0.05:5.00))
#Twaitfor_SH = 1.0;

tmake_SH_PRi = zeros(Float64, length(Nv));
tmake_SH_Ci = zeros(Float64, length(Nv));
#tmake_SH_Ci = zeros(Float64, length(Nv));# approx zero time
tmake_SH_PhiALL = zeros(Float64, length(Nv));
teig_SH = zeros(Float64, length(Nv));
tfixP_SH = zeros(Float64, length(Nv));
μSH = zeros(Float64, length(Nv));


domoreSH = true;
domoreSH_PR = true;
#@show domoreSH=false;
tmake_LR = zeros(Float64, length(Nv));
teig_LR = zeros(Float64, length(Nv));
tfixP_LR = zeros(Float64, length(Nv));
μLR = zeros(Float64, length(Nv));


f1 = (x) -> log(x) ./ log(10)
f2 = (x) -> log(x) ./ log(10)

#kNdisc=10
#@profview  
for kNdisc in vcat([1,1],1:length(Nv)) #the first is repated to get read of the first compliation time
    if kNdisc==1
        tmake_SH_PRi[kNdisc] = 0.0;
        tmake_SH_Ci[kNdisc] = 0.0;
        tmake_SH_PhiALL[kNdisc] = 0.0;
        teig_SH[kNdisc] = 0.0;
        tfixP_SH[kNdisc] = 0.0;
        μSH[kNdisc] = 0.0;
        domoreSH = true;
        domoreSH_PR = true;
        tmake_LR[kNdisc] = 0.0;
        teig_LR[kNdisc] = 0.0;
        tfixP_LR[kNdisc] = 0.0;
        μLR[kNdisc] = 0.0;
    end
    Ndisc = Nv[kNdisc]
    println([kNdisc / length(Nv), Ndisc])

    #@show Naver = maximum([ceil(4 - (log(Ndisc) / log(10))) * 2, 1])
    @show Naver = 10
    τmax = 2π # the largest τ of the system
    T = 2π #Principle period of the system (sin(t)=cos(t+T)) 
    mathieu_lddep = createMathieuProblem(3.0, 2.9, -0.45, 0.2, T=T) # LDDE problem for Hayes equation

    method = SemiDiscretization(0, 2 * pi / Ndisc) # 3rd order semi discretization with Δt=0.1

    Nsteps = Int((T + 100eps(T)) ÷ method.Δt)
    for _ = 1:Naver
 

        tmake_SH_PRi[kNdisc] += @elapsed resultPR = SemiDiscretizationMethod.calculateResults(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)
        if domoreSH_PR
            tmake_SH_Ci[kNdisc] += @elapsed mmpp2 = SemiDiscretizationMethod.DiscreteMappingSteps(resultPR)
            #tmake_SH_Ci[kNdisc] += @elapsed mmpp3 = SemiDiscretizationMethod.DiscreteMapping(mmpp2...)
            domoreSH_PR = tmake_SH_Ci[kNdisc] < Twaitfor_SH * Naver
        end
        if domoreSH
            println("domoreSH")
            #mappingFull = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system
            ##t = @benchmark DiscreteMapping_1step($mathieu_lddep, $method, $τmax, n_steps=$Nsteps, calculate_additive=true) #The discrete mapping of the system
            ##tmake_SH_Phi[kNdisc] = BenchmarkTools.median(t).time / 1e9
            ##μSH[kNdisc] = spectralRadiusOfMapping(mappingFull)
            ##t = @benchmark  spectralRadiusOfMapping($mappingFull)
            ##teig_SH[kNdisc] = BenchmarkTools.median(t).time / 1e9
            ##t = @benchmark  fixPointOfMapping($mappingFull);#xP =
            ##tfixP_SH[kNdisc] = BenchmarkTools.median(t).time / 1e9


            tmake_SH_PhiALL[kNdisc] += @elapsed mappingFull = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system
            teig_SH[kNdisc] += @elapsed μSH[kNdisc] = spectralRadiusOfMapping(mappingFull)
            #  tfixP_SH[kNdisc] += @elapsed xP = fixPointOfMapping(mappingFull)
            domoreSH = tmake_SH_PhiALL[kNdisc] < Twaitfor_SH * Naver
        end

        ##mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)#The discrete mapping of the system
        ##t = @benchmark DiscreteMapping_LR($mathieu_lddep, $method, $τmax, n_steps=$Nsteps, calculate_additive=true)#The discrete mapping of the system
        ##tmake_LR[kNdisc] = BenchmarkTools.median(t).time / 1e9
        ##μLR[kNdisc] = spectralRadiusOfMapping(mappingLR)
        ##t = @benchmark spectralRadiusOfMapping($mappingLR)
        ##teig_LR[kNdisc] = BenchmarkTools.median(t).time / 1e9
        ##t = @benchmark fixPointOfMapping($mappingLR);#xPLR = 
        ##tfixP_LR[kNdisc] = BenchmarkTools.median(t).time / 1e9

        tmake_LR[kNdisc] += @elapsed mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)#The discrete mapping of the system
        teig_LR[kNdisc] += @elapsed μLR[kNdisc] = spectralRadiusOfMapping(mappingLR)
        # tfixP_LR[kNdisc] += @elapsed xPLR = fixPointOfMapping(mappingLR)

    end

    tmake_SH_PRi[kNdisc] /= Naver
    tmake_SH_Ci[kNdisc] /= Naver
    #tmake_SH_Ci[kNdisc] /= Naver
    tmake_SH_PhiALL[kNdisc] /= Naver
    #tmake_SH_PhiALL[kNdisc] -=  tmake_SH_Ci[kNdisc]+tmake_SH_PRi[kNdisc]
    teig_SH[kNdisc] /= Naver

    tfixP_SH[kNdisc] /= Naver
    tmake_LR[kNdisc] /= Naver
    teig_LR[kNdisc] /= Naver
    tfixP_LR[kNdisc] /= Naver

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

#end
f1 = (x) -> log(x) ./ log(10)
f2 = (x) -> log(x) ./ log(10)
#f1 = (x) -> x
#f2 = (x) -> x
fplot=scatter
fplot! = scatter!
#fplot=plot
#fplot!=plot!
fplot(f2.(Nv), f1.(tmake_SH_PRi), labels="tmake_SH_PRi")
fplot!(f2.(Nv), f1.(tmake_SH_Ci), labels="tmake_SH_Ci")
#fplot!(f2.(Nv), f1.(tmake_SH_Ci), labels="tmake_SH_Ci")
fplot!(f2.(Nv), f1.(tmake_SH_PhiALL), labels="tmake_SH_PhiALL")
fplot!(f2.(Nv), f1.(teig_SH), labels="teig_SH: eigs(prod(...)) only")
#fplot!(f2.(Nv), f1.(tfixP_SH), labels="teig_SH: fix Point")
fplot!(f2.(Nv), f1.(tmake_LR), labels="tmake_LR", linewidth=2)
fplot!(f2.(Nv), f1.(teig_LR), labels="teig_LR: eigs(ΦR,ΦL) only", linewidth=2)
#fplot!(f2.(Nv), f1.(tfixP_LR))

#fplot!(xticks =collect( 10 .^(1:0.25:5))) 
#fplot!(yticks =collect( 10.0.^ (-5:1:5))) 
#fplot!(yaxis=:log10) 
#fplot!(xaxis=:log10)
fplot!(gridlinewidth=2)
display(
    fplot!(labels="teig_LR: fix Point",
        xlabel=L"number of steps, log_{10}(N)", ylabel=L"log_{10}(time)", legend=:bottomright)
)


end
5 + 5
#####@show norm(μSH - μLR)
#####
#####
#####plot(f2.(Nv), f1.(μSH))
#####plot!(f2.(Nv), f1.(μLR))
#####plot(f2.(Nv), (μSH .- μLR))
#####
#####
###### --------------------------------------------------------
#####
#####
#####Ndisc = 600
#####method = SemiDiscretization(1, 2 * pi / Ndisc) # 3rd order semi discretization with Δt=0.1
#####@elapsed mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)
#####@elapsed μLR = spectralRadiusOfMappingLR(mappingLR)
#####@elapsed xPLR = fixPointOfMappingLR(mappingLR)
#####
#####mappingLR.LmappingMX
#####mappingLR.RmappingMX
#####
######using TimerOutputs
######tmr = TimerOutput()
######DiscreteMapping_LR(mathieu_lddep,method,τmax,n_steps=Int((T+100eps(T))÷method.Δt),calculate_additive=true)
######show(tmr)
#####
######Windows: cmd
######set JULIA_NUM_THREADS=4
#####using Profile, PProf
#####
#####DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)
#####@profile DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)
#####pprof(; webport=58699)
#####
#####
#####using BenchmarkTools
######println("Drop")
######println("No-drop")
#####@benchmark mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)
######seconds=10
######@btime mappingLR=DiscreteMapping_LR(mathieu_lddep,method,τmax,n_steps=Int((T+100eps(T))÷method.Δt),calculate_additive=true)
#####
#####
#####
#####@time mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true);