# CoupledVdPBifurcation
Julia code for the examination of dissipatively coupled van der Pol equations in the form

![equation](https://user-images.githubusercontent.com/6823593/154071584-f720fe84-0247-4b2a-aeb5-52cd45ddce6a.png)



The code provides the following functionality:
1. 


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
