using LinearAlgebra, Arpack
using BenchmarkTools
#math
using DifferentialEquations
using Interpolations
using Peaks
using Polynomials

#plots
using Plots
pgfplotsx()
theme(:mute)
using ColorSchemes
cols=ColorSchemes.Spectral_11;
Plots.scalefontsizes(1.25);
using LaTeXStrings


include("vpdModule.jl")


function getPhiC(mu::Float64, Tstar::Float64)
    u0=[3; 0; 3; 0];
    dw=0.2; #mu=0.25; inherited from input
    p=(dw, mu);
    T=2*π;
    n=200;
    tspan=(0.0, n*T);  
    problem=ODEProblem(f, u0, tspan, p);
    t, x, dx, y, dy=vpdSolve(problem, true, 10);
    t_ind=findfirst(x->x==true, t.>Tstar);

    h=t[2]-t[1];
    p_approx=Int(round(T/h));
    step=2;
    phi=Int(round(p_approx/4));
    shifts=[maximum(x[t_ind+phi:end]-y[t_ind:end-phi]) for phi in 1:step:p_approx-1];
    C, phi_ind=findmin(shifts);
    phi=(1:step:p_approx-1)[phi_ind]*h
    phi=phi % (2*π);
    return phi, C
end

T=2*π;
getPhiC(0.25, 50*T)
getPhiC(1.5, 50*T)
getPhiC(6.0, 50*T)

mus=[range(0.22, 1.75, 100); range(1.8, 6, 20)];
Tstar=100*T;

Cs=zeros(size(mus, 1));
phis=zeros(size(mus, 1));

using Printf
for i in 1:size(mus,1)
    global n, tspan, ks, mus
    @time phis[i], Cs[i]=getPhiC(mus[i], Tstar);
    @printf "coupling: %f   |||  phase dif:  %f   :::  vert shift: %f\n " mus[i] phis[i] Cs[i]
end



#HERE LIES THE PLOTTER, ONLY AESTHETICS

ind=findmax(Cs[3:end])[2];
plot_font = "Computer Modern";

ψ=range(π/12+π/2, 0, 1000);
rad=0.25

plot()
scatter!(Cs[4:end], phis[4:end], xscale=:log10, yscale=:log10, color=:firebrick, labels="")
plot!(Cs[4:end], phis[4:end], xscale=:log10, yscale=:log10, alpha=0.2, lw=4, color=:firebrick, linestyle=:solid, labels="", )
scatter!([Cs[ind-1]], [phis[ind-1]], color=:black, markersize=7, labels="")
annotate!(Cs[ind-1]*1.1, 1.1*phis[ind-1], 
            text(L"$\mathbf{C(\mu^*)} \approx 0.083$", :left, 11
            )
            )
annotate!(Cs[ind-1]*1.1, 0.90*phis[ind-1], 
            text(L"$\mathbf{\varphi(\mu^*)}\approx 0.28$", :left, 11
            )
            )

annotate!(1.6*Cs[3], 0.8*phis[3], text(L"$\mathbf{\mu \approx \Delta\omega}$", :right, 10))
annotate!(0.9*Cs[end], 1.1*phis[end], text(L"$\mathbf{\mu \to \infty}$"*("\n")*L"$(C(\mu), \varphi(\mu) \to (0, 0)$", :right, 10))

plot!(10 .^ (rad*cos.(ψ) .- 1.75), 10 .^ (rad*sin.(ψ) .- 0.6), color=:black, lw=2, labels="")
scatter!( [10 .^ (rad*cos.(ψ[end]) .- 1.75)], [10 .^ (rad*sin.(ψ[end]) .- 0.6)], color=:black, markershape=:dtriangle, labels="")
annotate!(10 .^ (rad*cos.(ψ[1]) .- 1.75), 1.05*10 .^ (rad*sin.(ψ[1]) .- 0.6), text(L"\mathbf{\mu}\textrm{\; increases}", :left, :bottom, 10))

xlims!((0.0075, 0.2))
xlabel!(L"\textrm{V}\textrm{ertical \; displacement, }C(\mu)")
ylabel!(L"\textrm{P}\textrm{hase \; difference, }\varphi(\mu)")
plot!(size=(500, 500))


savefig("figure6.png")
savefig("figure6.svg")

savefig("Cphi.png")
savefig("Cphi.svg")

