# SemiDiscretizationMethod.jl

Julia package to investigate the behaviour of linear delay differential equations based on the book 
[Semi-Discretization for Time-Delay Systems (by Insperger and Stepan)](http://link.springer.com/10.1007/978-1-4614-0335-7).

This package provides a tool to approximate the stability properties and stationary behaviour of linear periodic delay systems of the form:

$$ \dot{\mathbf{x}}(t) = \mathbf{A}(t) \mathbf{x}(t) + \sum_{j=1}^g \mathbf{B}_j(t) \mathbf{x}(t-\tau_j(t))+\mathbf{c}(t)$$
by transforming the underlying differential equation into the mapping:
$$\mathbf{y}_{n+1} = \mathbf{F}_n\mathbf{y}_n+\mathbf{f}_n,$$

where $n$ is the discrete time ($t_n = n \Delta t$), $\mathbf{F}_n$ is the mapping matrix constructed using $\mathbf{A}(t)$, $\mathbf{B}(t)$ and $\tau_j(t)$, while the vector $\mathbf{y}_n$ is the discretized state space vector:

$$ \mathbf{y}_n = \left(\mathbf{x}(t_n)^\top, \mathbf{x}(t_{n-1})^\top,\ldots,\mathbf{x}(t_{n-r})\right)^\top\!.$$

Each coefficient matrices of delay differential equations are periodic, with a principle period of $T$, namely:
$\mathbf A(t)=\mathbf A(t+T),\; \mathbf B_j(t)=\mathbf B_j(t+T),\; \tau_j(t)=\tau_j(t+T)$) and $\mathbf{c}(t)=\mathbf{c}(t+T)$
Furthermore, the integer $r$ is chosen in a way, that $r\Delta t\geq \max_{t \in \left[0,T\right],j=1\ldots g}\tau_j(t)$
 (the discretized "history function" contains all possible delayed values).  

With the use of the discrete (one-step) mappings the Monodromy matrix (or One-Period mapping) can be determied:
$$\mathbf{\Phi}_n =\prod_{i=0}^{p-1}\mathbf{F}_{n+i}$$

The stability of the original system can be investigated (approximately), by the spectral radius $\rho$ of $\mathbf{\Phi}_n$:

$$\rho\left(\mathbf{\Phi}_n\right): \quad
    \begin{matrix}
    <1 & \Rightarrow & \text{the mapping is stable}\\
    >1 & \Rightarrow & \text{the mapping is unstable}
    \end{matrix}
    $$

Furthermore, the stationary solution can be determined by the periodic fix point (stationary orbit) of the mapping.

# Citing

If you use this package as part of your research, teaching, or other activities, we would be grateful if you could cite the book it is based on (BibTeX entry):
```
@book{Insperger2011,
address = {New York, NY},
author = {Insperger, Tam{\'{a}}s and St{\'{e}}p{\'{a}}n, G{\'{a}}bor},
doi = {10.1007/978-1-4614-0335-7},
isbn = {978-1-4614-0334-0},
publisher = {Springer New York},
series = {Applied Mathematical Sciences},
title = {{Semi-Discretization for Time-Delay Systems}},
url = {http://link.springer.com/10.1007/978-1-4614-0335-7},
volume = {178},
year = {2011}
}
```

If you are interested in the behaviour of your linear delay model in the presence of Gaussian white noise, please consider the [StochasticSemiDiscretizationMethod.jl](https://github.com/HTSykora/StochasticSemiDiscretizationMethod.jl) package.

# Updated: Multiplication Free Semi-Discretiztaion Method:
It is found that the One-Period Mapping can be describen by two matrices:
$$\mathbf{\Phi _L}_n \mathbf{y}_{n+1} = \mathbf{\Phi _R}_n\mathbf{y}_n+\mathbf{v}_n,$$

which can be created based on the one-step matrices. The quadratic time complexity is reduced to linear  with this modification.

# Usage with examples
## Installation
```julia
julia> ] add SemiDiscretizationMethod
```

## Hayes equations
$$\dot{x}(t) = a \,x(t) + b \,x(t-1) + 1.$$

Here 

$$ \mathbf{A}(t) \equiv \begin{bmatrix} a \end{bmatrix},
\quad \mathbf{B}_1(t) \equiv \begin{bmatrix}b\end{bmatrix},
\quad \tau_1(t) \equiv 1 , 
\quad \text{and} \quad \mathbf{c}(t) \equiv \begin{bmatrix} 1 \end{bmatrix}.$$

(Page 13 of the book referenced above)

```julia
using SemiDiscretizationMethod
```

```julia
function createHayesProblem(a,b)
    AMx =  ProportionalMX(a*ones(1,1));
    τ1=1. 
    BMx1 = DelayMX(τ1,b*ones(1,1));
    cVec = Additive(ones(1))
    LDDEProblem(AMx,[BMx1],cVec)
end
```

```julia
hayes_lddep=createHayesProblem(-1.,-1.); # LDDE problem for Hayes equation
method=SemiDiscretization(1,0.1) # 1st order semi discretization with Δt=0.1
τmax=1. # the largest τ of the system
mapping=DiscreteMapping_LR(hayes_lddep,method,τmax,n_steps=1,calculate_additive=true); #The discrete mapping of the system
```

```julia
@show spectralRadiusOfMapping(mapping); # spectral radius ρ of the mapping matrix (ρ>1 unstable, ρ<1 stable)
@show fixPointOfMapping(mapping); # stationary solution of the hayes equation (equilibrium position)

# spectralRadiusOfMapping(mapping) = 0.941189374166563
# fixPointOfMapping(mapping) = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
```
### Stability borders of the Hayes Equation
```julia
using MDBM

using Plots
gr();
using LaTeXStrings
```
```julia
method=SemiDiscretization(4,0.1);
τmax=1.

foo(a,b) = log(spectralRadiusOfMapping(DiscreteMapping_LR(createHayesProblem(a,b),method,τmax,
    n_steps=1))); # No additive term calculated

axis=[Axis(-15.0:15.,:a),
    Axis(-15.0:15.,:b)]

iteration=3;
stab_border_points=getinterpolatedsolution(solve!(MDBM_Problem(foo,axis),iteration));

scatter(stab_border_points...,xlim=(-15.,15.),ylim=(-15.,15.),
    label="",title="Stability border of the Hayes equation",xlabel=L"a",ylabel=L"b",
    guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)
```
![](./assets/HayesStability.png)
## Delay Mathieu equation

$$\ddot{x}(t) + a_1 \dot{x}(t)+(\delta + \varepsilon \cos(t))x(t)=\\ b_0 x(t-2\pi) + \sin(2t)$$

<!--Here 

$$ \mathbf{x}(t) = \begin{bmatrix} x(t) \\ \dot{x}(t) \end{bmatrix} , \quad
\mathbf{A}(t) = \begin{bmatrix} 0 & 1 \\ -\delta - \varepsilon \cos(t) & -a_1 \end{bmatrix},
\quad \mathbf{B}_1(t) = \begin{bmatrix}0 & 0 \\ b_0 & 0\end{bmatrix},
\quad \tau_1(t) \equiv 2\pi, 
\quad \text{and} \quad \mathbf{c}(t) = \begin{bmatrix} 0 \\ \sin(2t) \end{bmatrix}.$$ -->

(Page 77 of the book referenced above)

```julia
function createMathieuProblem(δ, ε, b0, a1; T=2π)
    AMx = ProportionalMX(t -> @SMatrix [0.0 1.0; -δ-ε*cos(2π / T * t) -a1])
    τ1 = t -> 2π # if function is needed, the use τ1 = t->foo(t)
    BMx1 = DelayMX(τ1, t -> @SMatrix [0.0 0.0; b0 0.0])
    cVec = Additive(t -> @SVector [0.0, sin(4π / T * t)])
    LDDEProblem(AMx, [BMx1], cVec)
end;
```
```julia
τmax=2π # the largest τ of the system
T=2π #Principle period of the system (sin(t)=sin(t+T)) 
mathieu_lddep=createMathieuProblem(3.,2.,-0.15,0.1,T=T); # LDDE problem for Hayes equation
method=SemiDiscretization(1,0.01) # 1st order semi discretization with Δt=0.01
# if P = τmax, then n_steps is automatically calculated
mapping=DiscreteMapping_LR(mathieu_lddep,method,τmax,
    n_steps=Int((T+100eps(T))÷method.Δt),calculate_additive=true); #The discrete mapping of the system

@show spectralRadiusOfMapping(mapping); # spectral radius ρ of the mapping matrix (ρ>1 unstable, ρ<1 stable)
# spectralRadiusOfMapping(mapping) = 0.5131596340374617
fp=fixPointOfMapping(mapping); # stationary solution of the hayes equation (equilibrium position)

# fixPointOfMapping plotted
using Plots
gr();
using LaTeXStrings
n_steps=Int((T+100eps(T))÷method.Δt)
t_plot=0.0:method.Δt:((size(fp,1)/2-1)*method.Δt)
plot(t_plot,fp[1:2:end],
    xlabel=L"-s",title=L"t \in [nT,(n+1)T],\quad n \to \infty",guidefontsize=14,linewidth=3,
    label=L"x(t-s)",legendfontsize=11,tickfont = font(10))
plot!(t_plot,fp[2:2:end],
    xlabel=L"-s",linewidth=3,
    label=L"\dot{x}(t-s)")
plot!(t_plot,sin.(2 .* t_plot),linewidth=3,label=L"\sin(2t)")
```
![](./assets/MathieuStationary.png)
### Stability Chart of the delayed Mathieu equation
```julia
using MDBM
a1=0.1;
ε=1;
τmax=2π;
T=1π;
method=SemiDiscretization(2,T/40);

foo(δ,b0) = log(spectralRadiusOfMapping(DiscreteMapping_LR(createMathieuProblem(δ,ε,b0,a1,T=T),method,τmax,
    n_steps=Int((T+100eps(T))÷method.Δt)))); # No additive term calculated
using MDBM
axis=[Axis(-1:0.2:5.,:δ),
    Axis(-2:0.2:1.5,:b0)]
    
iteration=3;
stab_border_points=getinterpolatedsolution(solve!(MDBM_Problem(foo,axis),iteration));

scatter(stab_border_points...,xlim=(-1.,5),ylim=(-2.,1.5),
    label="",title="Stability border of the delay Mathieu equation",xlabel=L"\delta",ylabel=L"b_0",
    guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)
```
![](./assets/MathieuStability.png)
