using LinearAlgebra, Arpack
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

include("vpdModule.jl")


function getDistsSmooth(dists_x, t, pNum)
    dists_x_sm=[mean((dists_x[maximum([1, i - pNum ÷ 2]) : minimum([i+ pNum ÷ 2, size(dists_x, 1)])])) for i in 1:size(dists_x,1)];
   
    t_sm=t;
    return dists_x_sm, t_sm
end

function coef(x, y)
    n=size(x, 1);
    xy=sum(x .* y); sx=sum(x); sy=sum(y);
    sx2=sum(x .* x);
    noma=sy*sx2-sx*xy;
    nomb=n*xy-sum(x)*sum(y);
    denom=n*sx2-sum(x)*sum(x);
    return [noma/denom; nomb/denom]
    
end


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


#SEVERAL COUPLINGS


mus=[range(0.25, 0.40, 16); range(0.41, 0.80, 80); range(0.81, 1.1, 15); range(sqrt(1.15), sqrt(10), 30) .^ 2];
num=7;
share002=zeros(size(mus, 1), num);


for i in 1:size(mus, 1)
    global mus, share002, share005, share01, share001
    T=2*π;
    n=100;
    tspan=(0.0, n*T);
    Δω=0.2; μ=mus[i]; p=(Δω, μ);
    @time begin
        u0=[3; 0; 3; 0];
        problem=ODEProblem(f, u0, tspan, p);
        mult=10;
        t, x, dx, y, dy=vpdSolve(problem, true, mult);
        γ_x, γ_dx, γ_y, γ_dy=getLimCycleNaive(t, x, dx, y, dy);
        γ=(γ_x, γ_dx, γ_y, γ_dy);

        D_0=0.5; thr=0.01;

        pNum=size(γ_x, 1);
        num=7;
        for j in 2:2:num-1
            ϕ=π*j/num;
            k=tan(ϕ);
            indx=findmin(abs.(γ_dx ./ γ_x .- k) )[2];
            (2* indx > pNum+12) ? (indx-= Int(round(pNum/2))) : 0;
        
            u0=getNewDot(D_0, γ_x[indx], γ_dx[indx], γ_y[indx], γ_dy[indx], p);
    
            n2=10;
            mult2=mult;

            tspan2=(0.0, n2*T);
            problem2=ODEProblem(f, u0, tspan2, p);
            t2, x, dx, y, dy=vpdSolve(problem2, true, mult2);
        
            γ_x2, γ_dx2, γ_y2, γ_dy2=getLimCycleNaive(t2, x, dx, y, dy);
            pNum2=size(γ_x2, 1);
            
            till=8;
        
            radius=20;
            dists_x, dists_y=getDists(x[1:Int(round(till*size(x, 1)/n2))],
                dx[1:Int(round(till*size(x, 1)/n2))],
                y[1:Int(round(till*size(x, 1)/n2))],
                dy[1:Int(round(till*size(x, 1)/n2))],
                γ, radius
            );
        
            dists_x3, t32=getDistsSmooth(dists_x, t2[1:Int(round(till*size(x, 1)/n2))], pNum2);
            dists_x3 /= D_0; 
        
            dists_x3, t32=getDistsSmooth(dists_x3, t32, pNum2);
            thr= 10*1e-3; indx=findfirst(x->x<thr, dists_x3); share002[i, j+1]=indx/pNum2;
        end
    end
    @printf "coupling: %f (%d) \n" round(mus[i]; digits=2) i;
end

# HERE START THE PLOTTING; ONLY AESTHETICS

share002[:, [3, 7]].+=10
pgfplotsx()
plot()
plot!(mus, share002[:, [3, 5, 7]], lw=3, alpha=1.0, labels=[L"\phi=\frac{2\pi}{7}"   L"\phi=\frac{4\pi}{7}"  L"\phi=\frac{6\pi}{7}" ], color=[cols[4] cols[10] cols[1]], legend_position=:topleft, legend_columns=-1, legendfont=font(16))
xlabel!(L"\textrm{C}\textrm{oupling,} \, \mu")
ylabel!(L"\textrm{T}\textrm{ime to  the  LC,} \, \tau_{LC}/(2\pi)")
share002[:, [3, 7]].-=10;
plot!(mus[15:80], 
        share002[15:80, [3, 5, 7]], 
        lw=3, alpha=1.0,
        labels="",
        marker=(:circle, 3),
        inset=(1, bbox(0.5, 0.05, 0.5, 0.55, :bottom, :left)),
        subplot=2,
        xtickfont = font(10),
        ytickfont = font(10),
        framestyle = :box, color=[cols[4] cols[10] cols[1]]
)
plot!(size=(900, 350), left_margin = 8Plots.mm, bottom_margin=8Plots.mm)
ylims!(0.9, 6, sp=1)


savefig("figure3.png")
savefig("figure3.svg")

savefig("t2LC_0407_big_3phases+inset.png")
savefig("t2LC_0407_big_3phases+inset.svg")




