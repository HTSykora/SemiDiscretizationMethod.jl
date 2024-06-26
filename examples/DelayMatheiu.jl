


using SemiDiscretizationMethod

function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    τ1 = t -> 2π # if function is needed, the use τ1 = t->foo(t)
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, sin(4π / T * t)])
    LDDEProblem(AMx, [BMx1], cVec)
end;

τmax = 2π # the largest τ of the system
T = 2π #Principle period of the system (sin(t)=sin(t+P)) 
mathieu_lddep = createMathieuProblem(3.0, 0.2, -0.15, 0.1, T=T); # LDDE problem for Hayes equation
method = SemiDiscretization(0, 0.05) # 3rd order semi discretization with Δt=0.1
# if T = τmax, then n_steps is automatically calculated
@time mapping = DiscreteMapping(mathieu_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true); #The discrete mapping of the system
@time mapping_LR = DiscreteMapping_LR(mathieu_lddep, method, τmax,
    n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true); #The discrete mapping of the system

# spectral radius ρ of the mapping matrix (ρ>1 unstable, ρ<1 stable)
@time @show spectralRadiusOfMapping(mapping);
@time @show spectralRadiusOfMapping(mapping_LR);

# stationary solution of the hayes equation (equilibrium position)
@time fp = fixPointOfMapping(mapping);
@time fp_LR = fixPointOfMapping(mapping_LR);



using Plots
gr();
using LaTeXStrings

r = length(fp[1:2:end])
plot((1:r) .* method.Δt, fp[2:2:end], xlabel="-s", label=L"\dot{x}(t-s)")
plot!((1:r) .* method.Δt, fp[1:2:end], xlabel="-s", label=L"x(t-s)")

p = length(fp_LR[1:2:end])
plot!((0:p-1) .* method.Δt, fp_LR[2:2:end], xlabel="-s", label=L"\dot{x}_{LR}(t-s)")
plot!((0:p-1) .* method.Δt, fp_LR[1:2:end], xlabel=L"-s", title=L"t \in [nP,(n+1)P],\quad n \to \infty", guidefontsize=14, linewidth=3, label=L"x_{LR}(t-s)", legendfontsize=11, tickfont=font(10))

plot!(0.0:method.Δt:T, sin.(4pi * (0.0:method.Δt:T) ./ T), label=L"sin(4\pi  t /T)")


#--------- Stability map ----------
using MDBM

a1 = 0.1;
ε = 1;
τmax = 2π;
T = 1π;
method = SemiDiscretization(2, T / 40);

foo(δ, b0) = log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ, ε, b0, a1, T=T), method, τmax,
        n_steps=Int((T + 100eps(T)) ÷ method.Δt)), nev=1, tol=1e-3)); # No additive term calculated

axis = [Axis(-1:0.5:5.0, :δ),
    Axis(-2:0.5:1.5, :b0)]

iteration = 4;
@time stab_border_points = getinterpolatedsolution(solve!(MDBM_Problem(foo, axis), iteration));

scatter(stab_border_points...,
    label="", title="Stability border of the delay Mathieu equation", xlabel=L"\delta", ylabel=L"b_0",
    guidefontsize=14, tickfont=font(10), markersize=2, markerstrokewidth=0)



#--------- Stability map of the Mathieu equation (no delay)----------
using MDBM

a1 = 0.01;
τmax = 2π;
T = 2π;
b0 = 0.0
method = SemiDiscretization(2, T / 40);

foo(δ, ε) = log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ, ε, b0, a1, T=T), method, τmax,
        n_steps=Int((T + 100eps(T)) ÷ method.Δt)), nev=1, tol=1e-4)); # No additive term calculated

axis = [Axis(-2:1:5.0, :δ),
    Axis(-0.01:1:5, :ε)]

iteration = 5;
@time stab_border_points = getinterpolatedsolution(solve!(MDBM_Problem(foo, axis), iteration));

scatter(stab_border_points...,
    label="", title="Stability border of the delay Mathieu equation", xlabel=L"\delta", ylabel=L"\epsilon",
    guidefontsize=14, tickfont=font(10), markersize=2, markerstrokewidth=0)


#--------- 3D stability map ----------
using MDBM
plotly()
a1 = 0.01;
τmax = 2π;
T = 2π;
method = SemiDiscretization(2, T / 40);

foo(δ, b0, ε) = log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ, ε, b0, a1, T=T), method, τmax,
        n_steps=Int((T + 100eps(T)) ÷ method.Δt)), nev=1, tol=1e-4)); # No additive term calculated

axis = [Axis(-2:0.5:5.0, :δ),
    Axis(-2:0.5:1.5, :b0),
    Axis(-0.01:0.5:5, :ε)]

iteration = 2;
@time stab_border_points = getinterpolatedsolution(solve!(MDBM_Problem(foo, axis), iteration));

scatter(stab_border_points...,
    label="", title="Stability border of the delay Mathieu equation", xlabel=L"\delta", ylabel=L"b_0", zlabel=L"\epsilon",
    guidefontsize=14, tickfont=font(10), markersize=1, markerstrokewidth=0)

