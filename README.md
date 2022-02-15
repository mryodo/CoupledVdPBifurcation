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
| `f(du, u, p, t)`  | `du` -- LHS
`u`, `t` -- variables
`p`  -- collection of the parameters, `tuple` | RHS of the system  | Template RHS for `DifferentialEqualtions` module  |
| `vpdSolve(problem::ODEProblem, interp::Bool, mult::Int64)`  | Content Cell  | Content Cell  | Content Cell  |
| `getLimCycleNaive(t, x, dx, y, dy)` | | | |
| `getDD(x::Float64, dx::Float64, y::Float64, dy::Float64, p::Tuple{Float64, Float64})` | | | |
| `getNewDot(eps::Float64, x::Float64, dx::Float64, y::Float64, dy::Float64, p::Tuple{Float64, Float64})` | | | |
| `dsquare(x,y)` | | | |
| `minimumdistance(x,y)` | | | |
| `getDists(x, dx, y, dy, γ, radius::Int)` | | | |
