# CoupledVdPBifurcation
Julia code for the examination of dissipatively coupled van der Pol equations in the form
$$
\begin{cases}
  \ddot{x}-(1-x^2)\dot x +(1-\Delta \omega)x - \mu  ( \dot x -\dot y)=0 \\
    \ddot{x}-(1-y^2)\dot y +(1+\Delta \omega)y - \mu  ( \dot y -\dot x)=0
\end{cases}
$$

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
