using LinearAlgebra
using Polynomials
using Plots; gr()

# Bottom topography function
function bottom(x)
  b=zeros(size(x)); N=length(x);
  for j=1:N
    if x[j] > 4.5 && x[j] < 5.5
      b[j] = cos(pi*x[j]); 
    end
  end
  return b
end

# Setting up the conditions we will need to implement
x_range=[0,10]; T=10; M=5;

# Building the matrix E
Phi = [Polynomial([1]),
       Polynomial([0, 1]),
       Polynomial([-1/2, 0, 3/2]),
       Polynomial([0, -3/2, 0, 5/2]),
       Polynomial([3/8, 0, -30/8, 0, 35/8])];
for k=1:M
  Phi[k] = Phi[k] / sqrt( integrate(0.5*Phi[k]*Phi[k], -1, 1) );
end

include("build_E.jl")
E = build_E(Phi);

MAX_M = 4;
cell_sizes = [100, 200, 400];
#mean_err = zeros(length(cell_sizes), MAX_M);
#sd_err = zeros(length(cell_sizes), MAX_M);
#
#include("burgers_wbsg.jl")
#for m=1:MAX_M, j=1:length(cell_sizes)
#  println("Testing with (M, Nx) = ($m, $(cell_sizes[j]))");
#  En = E[1:m, 1:m, 1:m];
#  cell_count = cell_sizes[j]; dx = 10/cell_count; 
#
#  x=0.5*dx:dx:10-0.5*dx;
#
#  # Building u0_h
#  u0_h = zeros(m,cell_count+1);
#  u0_h[1,1] = 2;
#
#  # Buildng b_h
#  coefficients=zeros(m,1);
#  pz = Polynomial([2,1]);
#  for k=1:m
#    coefficients[k] = integrate(0.5*pz*Phi[k], -1, 1);
#  end
#  b = bottom(x);
#  b_h = hcat(zeros(m,1), coefficients*b');
#
#  steady_mean = 2 .- 2*b;
#  steady_sd   = abs.(b) / sqrt(3);
#
#  # Well Balanced Scheme
#  @time u_wb = burgers_wbsg(u0_h, b_h, En, T, dx, 0.0025/8);
#  mean_wb = u_wb[1,2:end]; 
#  sd_wb = norm.( eachcol(u_wb[2:end, 2:end]) );
#
#  mean_err[j,m] = dx*norm(steady_mean .- mean_wb, 1);
#  sd_err[j,m] = dx*norm(steady_sd .- sd_wb, 1);
#end

p1 = plot(yaxis=:log);
p1 = plot!(cell_sizes, mean_err[:,end], c=:orange, lw=3, label="")
p1 = plot!(cell_sizes, sd_err[:,end], c=:green, ls=:dash, lw=3, label="")

p2 = plot(yaxis=:log);
p2 = plot!(0:MAX_M-1, mean_err[end,:], c=:orange, lw=3, label="Mean Errors")
p2 = plot!(0:MAX_M-1, sd_err[end,:], c=:green, ls=:dash, lw=3, label="SD Errors")

plot(p1, p2, layout=(1,2), size=(800,400));
savefig("Figures/convergence.pdf")
