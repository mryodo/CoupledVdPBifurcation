[![DOI](https://zenodo.org/badge/459579858.svg)](https://zenodo.org/badge/latestdoi/459579858)


# Coupled VdP Bifurcation [Julia code]
Julia code for the examination of dissipatively coupled van der Pol equations in the form

![equation](https://user-images.githubusercontent.com/6823593/154071584-f720fe84-0247-4b2a-aeb5-52cd45ddce6a.png)

where ![equation-2](https://user-images.githubusercontent.com/6823593/154071905-39ce2b80-7e41-4129-b425-4ee314170cf2.png) and ![equation-3](https://user-images.githubusercontent.com/6823593/154071976-93436c32-83d4-4ea7-9548-cf33f8586192.png) are given parameters.

**For the more detailed discussion and theoretical analysis of the system, please, consult the paper [TBD][TBD].**

The code provides the following functionality:

1. High tolerance solution of the VdP system through `Tsit5()` routine with interpolation to the unfiorm mesh
2. Extraction of the naive limit cycle out of the solution
3. Perturbation of the LC with the given distance to the LC
4. Numerical definition of the vertical displacement and the phase different of the system
5. Numerical definition of the time to the limit cycle
6. Plotting routine of the paper [TBD][TBD]


### Necessary header and dependencies

We advice the following list of `Julia` modules for the successful run of the code:

```julia
using BenchmarkTools
#math
using DifferentialEquations
using Interpolations
using Peaks
using Polynomials
using Statistics
#plots
using Plots
pgfplotsx()
theme(:mute)
using ColorSchemes
cols=ColorSchemes.Spectral_11;
Plots.scalefontsizes(1.5);

using Printf
using LaTeXStrings
```


## `VdPModule.jl` functions

| function  | Inputs | Outputs| Notes |
| ------------- | ------------- | ------------- | ------------- |
| `f(du, u, p, t)`  | `du` – LHS <br> `u`, `t` – variables <br> `p`  – collection of the parameters, `tuple` | RHS of the system  | Template RHS for `DifferentialEqualtions` module  |
| `vpdSolve(problem::ODEProblem, interp::Bool, mult::Int64)`  | `ODEProblem` input (see below) <br> `interp` – keep `true` for now, check if the mesh should be uniform <br> `mult` – multiplier of number of points for interpolation  | `t, x, dx, y, dx` arrays  | Wrapper on Solver + interpolation to the uniform mesh  |
| `getLimCycleNaive(t, x, dx, y, dy)` | solution of `vpdSolve` | `γ_x, γ_dx, γ_y, γ_dy` – LC | Naive LC exteractor, returns last cycle |
| `getDD(x::Float64, dx::Float64, y::Float64, dy::Float64, p::Tuple{Float64, Float64})` | `p` – `tuple` of parameters | `ddx`, `ddy` | Extractor of the second derivatives in the given point |
| `getNewDot(eps::Float64, x::Float64, dx::Float64, y::Float64, dy::Float64, p::Tuple{Float64, Float64})` | `eps` – the distanceto the LC <br> `x, dx, y, dy` – point of the LC being perturbed <br> `p` – `tuple` of parameters | `u0` – initial point for the solver | Extractor of the initial perturbation exactly `eps`-away from the LC |
| `dsquare(x,y)` | 2 2-d points|  | returns euclidian distance between to points |
| `minimumdistance(x,y)` | `x` – given point, 2x1 <br> `y` – array of points, n x 2  | minimal distance from `x` to `y` and the index of the smallest distance | |
| `getDists(x, dx, y, dy, γ, radius::Int)` | `x, dx, y, dy` – output of `vdpSolve`<br> `γ` – LC <br> radius – technical; how far away from the previous point we search the minimal distance for the new point   | `dists_x`, `dists_y` – arrays of distances from the solution to the corresponding LC | |

### Minimal run of the solver
```julia
T=2*π;
n=100;
tspan=(0.0, n*T);
Δω=0.2; μ=0.6; p=(Δω, μ);

u0=[3; 0; 3; 0];
problem=ODEProblem(f, u0, tspan, p);
t, x, dx, y, dy=vpdSolve(problem, true, 10);
```

## `PhiC_calc.jl` functions
| function  | Inputs | Outputs| Notes |
| ------------- | ------------- | ------------- | ------------- |
| `getPhiC(mu::Float64, Tstar::Float64)` | `mu` – coupling <br> `Tstar` – moment of calculation (consider the paper) | `phi` – phase difference <br> `C` – vertical displacement | Returns the vertical displacement and phase difference directly from definition |

### Minimal run of the Phi/C extraction

```julia
T=2*π;
mu=1.5;
getPhiC(mu, 100*T)
```

**Additionally, `PhiC_calc.jl` performs run of `getPhiC` through several couplings and plots the `C\phi` diagram from the paper.**

## `time2LC_calc..jl` functions

| function  | Inputs | Outputs| Notes |
| ------------- | ------------- | ------------- | ------------- |
| `getDistsSmooth(dists_x, t, pNum))` | `dists_x` – distances to `x/dx`-LC, output from `getDists` <br> `t` – time, output from `vpdSolve`, `pNum` – number of dots per LC, `size(γ_x, 1)` | `dists_x_sm` – smoothed distances <br> `t` – correponding (truncated) time | Smoothes by averaging the distances to the LC in the period window |
| `coef(x, y)` |`x` and `y` – 1-d array | `k` – optimal slope by pseudo-inverse |  `Julia`-optimal linreg slope  |

### Minimal run to get distance to the LC
```julia
T=2*π;
n=100;
tspan=(0.0, n*T);
Δω=0.2; μ=1.5; p=(Δω, μ);

#solve the system to get the LC
u0=[3; 0; 3; 0];
problem=ODEProblem(f, u0, tspan, p);
mult=10;
t, x, dx, y, dy=vpdSolve(problem, true, mult);
γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
γ=(γ_x, γ_dx, γ_y, γ_dy);

# get the perturbed point
D_0=0.1;
indx=pNum ÷ 3;
u0=getNewDot(D_0, γ_x[indx], γ_dx[indx], γ_y[indx], γ_dy[indx], p);

#solve the system from the perturbed point
n2=10;
mult2=mult;
tspan2=(0.0, n2*T);
problem2=ODEProblem(f, u0, tspan2, p);
t2, x, dx, y, dy=vpdSolve(problem2, true, mult2);
γ_x2, γ_dx2, γ_y2, γ_dy2=getLimCycleNaive(t2, x, dx, y, dy);
pNum2=size(γ_x2, 1);

#get the distances
till=8;
radius=20;
dists_x, dists_y=getDists(x[1:Int(round(till*size(x, 1)/n2))],
    dx[1:Int(round(till*size(x, 1)/n2))],
    y[1:Int(round(till*size(x, 1)/n2))],
    dy[1:Int(round(till*size(x, 1)/n2))],
    γ, radius
);

#smooth the distances
dists_x3, t32=getDistsSmooth(dists_x, t2[1:Int(round(till*size(x, 1)/n2))], pNum2);
dists_x3 /= D_0; 
dists_x3, t32=getDistsSmooth(dists_x3, t32, pNum2);

#find the moment when the distance is small enough
thr= 10*1e-3; 
indx=findfirst(x->x<thr, dists_x3); 
time2LC=indx/pNum2;
```
**Additionally, `time2LC_calc.jl` performs run of extracting the distance to the LC through several couplings and plots the `time to the LC` diagram from the paper for different phases of the initially pertrubed point.**

## Illustrative plotting examples 

### Examples of the Solutions

```julia
T=2*π;
n=100;
tspan=(0.0, n*T);
Δω=0.2; μ=0.25; p=(Δω, μ);

u0=[3; 0; 1; 0];
problem=ODEProblem(f, u0, tspan, p);
t, x, dx, y, dy=vpdSolve(problem, true, 10);
γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
γ=(γ_x, γ_dx, γ_y, γ_dy);
pNum=size(γ_x, 1);

u=0.5*(x+y); v=0.5*(x-y);
du=0.5*(dx+dy); dv=0.5*(dx-dy);

pgfplotsx()

plot(layout=grid(2, 2, heights=[0.5 ,0.5, 0.5, 0.5], widths=[0.75, 0.25, 0.75, 0.25]))
plot!(t[end-5*pNum:end], x[end-5*pNum:end], color=cols[11], sp=1, lw=3, labels=L"x/\dot{x}")
plot!(t[end-5*pNum:end], y[end-5*pNum:end], color=cols[2], sp=1, lw=3, labels=L"y/\dot{y}",)
title!(L"\textrm{S}\textrm{olutions\; in \;time}", titlefontsize=14, sp=1)
ylabel!(L"\textrm{oscillators \; }\{x(t), y(t)\}", yguidefontsize=12, sp=1, )
plot!(lest_margin=8Plots.mm, sp=1)

plot!(t[end-5*pNum:end], u[end-5*pNum:end], color=cols[9], sp=3, lw=3, labels=L"u/\dot{u}")
plot!(t[end-5*pNum:end], v[end-5*pNum:end], color=cols[4], sp=3, lw=3, labels=L"v/\dot{v}")
ylabel!(L"\textrm{half-sum/diff \; }\{u(t), v(t)\}", yguidefontsize=12, sp=3, )
xlabel!(L"\textrm{time,\;} t", xguidefontsize=14, sp=3,)

plot!(x[1:10*pNum], dx[1:10*pNum], line=(:dot), sp=2, color=cols[11], alpha=0.75, lw=3, labels="")
plot!(y[1:10*pNum], dy[1:10*pNum], line=(:dot), sp=2, color=cols[2], alpha=0.75, lw=3,  labels="")
ylabel!(L"\textrm{derivatives}",yguidefontsize=12, sp=2, )
title!(L"\textrm{P}\textrm{hase\; portraits}", titlefontsize=14, sp=2)

plot!(u[1:10*pNum], du[1:10*pNum], line=(:dot), sp=4, color=cols[9], alpha=0.75, lw=3, labels="")
plot!(v[1:10*pNum], dv[1:10*pNum], line=(:dot), sp=4, color=cols[4], alpha=0.75, lw=3, labels="", legend_position=:bottomright, legendmarkerstroke=1)
ylabel!(L"\textrm{derivatives}",yguidefontsize=12, sp=4, )
xlabel!(L"\textrm{functions}", xguidefontsize=14, sp=4, )
plot!(size=(1000, 500), bottom_margin=4Plots.mm, left_margin=4Plots.mm)
```
![figure1](https://user-images.githubusercontent.com/6823593/154083776-9789dbb4-cdde-4708-9759-fbe407007c15.png)

### Amplitudes and half-sum/half-difference example

```julia
amps=zeros(4);
mus=[0.25, 0.5, 1, 3];
refcols=[cols[1]; cols[11]; cols[9]; cols[10] ];
plot(layout=@layout [
    a{0.83w, 1.0h} [b{0.9h}; _]])
for i in 1:size(mus,1)
    global mus, amps
    μ=mus[i];
    t=2*π;
    n=100;
    tspan=(0.0, n*T);
    Δω=0.2; p=(Δω, μ);

    u0=[3; 0; 1; 0];
    problem=ODEProblem(f, u0, tspan, p);
    t, x, dx, y, dy=vpdSolve(problem, true, 10);
    γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
    γ=(γ_x, γ_dx, γ_y, γ_dy);
    pNum=size(γ_x, 1);
    
    u=0.5*(x+y); v=0.5*(x-y);
    du=0.5*(dx+dy); dv=0.5*(dx-dy);
    amps[i]=maximum(v[end-pNum:end]);
    plot!(t[1:10*pNum], v[1:10*pNum], lw=3, color=refcols[i], sp=1, labels="")
    scatter!([mus[i]], [amps[i]], xscale=:log10, yscale=:log10, sp=2,
        labels=L"\mu=%$(μ)", 
        legend_position=:top,
        marker=(:circle, 5),
        xtickfont = font(10),
        ytickfont = font(10),
        color=refcols[i]
    )
end
ylims!(10^(-1.2), 10^(1.9), sp=2)
xlabel!(L"\mathrm{time, \; } t", sp=1)
ylabel!(L"\mathrm{half-difference,\;} v(t)", sp=1)

yticks!([0.1, 1], sp=2)
xlabel!(L"\mathrm{coupling, \; } \mu", sp=2, xguidefontsize=11)
ylabel!(L"\mathrm{amplitude}", sp=2, yguidefontvalign =:bottom, yguidefontsize=11)
plot!(size=(900, 400), bottom_margin=5Plots.mm, left_margin=4Plots.mm)
```

![figure2](https://user-images.githubusercontent.com/6823593/154084090-2bbd0822-b2af-48d0-9e16-63a9c0751f76.png)

