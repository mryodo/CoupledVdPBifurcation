# CoupledVdPBifurcation
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


## `VdPMoudule.jl` functions

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
