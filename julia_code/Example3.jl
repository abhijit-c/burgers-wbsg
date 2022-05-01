using LinearAlgebra
using QuadGK
using Polynomials
using Plots

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
# Use numerical quadrature to estimate mean
function u_mean(x)
  b = bottom(x); N = length(b);
  mean = zeros(N);
  for i=1:N
    b_i = b[i];
    f(z) = 0.5*(3*(2-(2+z)*b_i))^(1/3);
    mean[i], err = quadgk(f,-1,1);
  end
  return mean
end
# Use numerical quadrature to estimate sd
function u_sd(x)
  b = bottom(x); mean = u_mean(x); N = length(mean);
  sd = zeros(N);
  for i=1:N
    b_i = b[i];
    f(z) = 0.5*(3*(2-(2+z)*b_i))^(1/1.5);
    sd[i], err = quadgk(f,-1,1);
    sd[i] = sd[i] - mean[i]^2;
  end
  return sqrt.(sd);
end

# Setting up the conditions we will need to implement
cell_count = 100; dx = 10/cell_count; x_range=[0,10]; T=10; M=5;
x=0.5*dx:dx:10-0.5*dx;

# Building the matrix E
Phi = [Polynomial([1]),
       Polynomial([0, 1]),
       Polynomial([-1/2, 0, 3/2]),
       Polynomial([0, -3/2, 0, 5/2]),
       Polynomial([3/8, 0, -30/8, 0, 35/8])];
for k=1:M
  Phi[k] = Phi[k] / sqrt( integrate(0.5*Phi[k]*Phi[k], -1, 1) );
end

include("build_E4.jl")
E = build_E4(Phi);

# Building u0_h

# Our initial condition is u0 = 0 so coefficient functions should be
# uniformly 0 with the exception of the boundary condition

# First column corresponds to cell 0 which reflects our boundary condition
# of u = 2. This boundary condition is constant in the random variable so
# E[2 Phi_k] = 2 E[Phi_k] = 2 E[Phi_1 Phi_k] = 2 dirac_{1k}, so only the
# first entry of the first column should be nonzero
u0_h = zeros(M,cell_count+1);
u0_h[1,1] = 6^(1/3);

# Buildng b_h
#
# We now construct the bottom topography function. The bottom topography
# can be written as (2+z)*bottom(x) where bottom is the function defined
# below. 
coefficients=zeros(M,1);
pz = Polynomial([2,1]);
for k=1:M
  coefficients[k] = integrate(0.5*pz*Phi[k], -1, 1);
end
b = bottom(x);
b_h = hcat(zeros(M,1), coefficients*b');

steady_mean = u_mean(x);
steady_sd = u_sd(x);
                                            
# Well Balanced Scheme

include("u4_wbsg.jl")
u_wb = u4_wbsg(u0_h, b_h, E, T, dx, 0.0025/8);
mean_wb = u_wb[1,2:end]; 
sd_wb = norm.( eachcol(u_wb[2:end, 2:end]) );

include("u4_nwbsg.jl")
u_nwb = u4_nwbsg(u0_h, b_h, E, T, dx, 0.0025/8);
mean_nwb = u_nwb[1,2:end]; 
sd_nwb = norm.( eachcol(u_nwb[2:end, 2:end]) );

plot(x, steady_mean, c=:black, ls=:dash, lw=1.5, label="Steady State");
scatter!(x, mean_wb, shape=:hexagon, ms=8, label="WBSG");
scatter!(x, mean_nwb, shape=:utriangle, ms=5, label="NWBWG");
savefig("Figures/mean_u4_sgwb.pdf")

plot(x, steady_sd, c=:black, ls=:dash, lw=1.5, label="Steady State");
scatter!(x, sd_wb, shape=:hexagon, ms=8, label="WBSG");
scatter!(x, sd_nwb, shape=:utriangle, ms=5, label="NWBWG");
savefig("Figures/sd_u4_sgwb.pdf")
