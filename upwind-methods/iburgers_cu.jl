using LinearAlgebra;
"""
  iburgers_cu(xx, u0, T, dx, dt)

Solve the inviscid Burgers equation
\\[ u_t + (u^2/2)_x = 0 (*) \\]
with some initial condition u0(x) up to time T on spatial domain [a,b]. We
employ a conservative second order central upwind scheme along with a 3rd order
SSP Runge-Kutta method to numerically approximate the solution.
...
# Arguments
x_range Tuple{Float64,Float64}: Spatial domain to solve PDE over.
u0 Function: Initial condition u0(xl,xr) returns the average of u0 over [xl,xr].
T Float64: Final time to solve PDE to.
nc Int64: Number of cells to use in the spatial discretization.
...
"""
function iburgers_cu(x_range, u0, T, nc)

  # Spatial Discretization Setup
  a = x_range[1]; b = x_range[end]; 
  dx = (b-a) / nc;
  xx = [a + (j-1/2)*dx for j=0:nc+1]; #xx[k] is the midpoint of cells[k]
  cells = [a + (j-1)*dx for j=0:nc+2]; #Cells[1] & Cells[end] are ghost cells

  curr_t = 0.0;
  u = u0.(cells[1:end-1], cells[2:end]);
  u1 = zeros(length(u));
  while curr_t < T #2nd order SSP scheme.
    rhs, dt = burgers_reconstruct(dx, u, true);
    u1[2:end-1] = u[2:end-1] .+ dt*rhs;
    u1[1] = u[1]; u1[end] = u[end];

    rhs, _ = burgers_reconstruct(dx, u1, false);
    u[2:end-1] = 0.5*(u[2:end-1] .+ u1[2:end-1]) .+ 0.5*dt*rhs;
    # No need to update u ghost cells, they haven't been touched.
    curr_t = curr_t + dt;
  end
  return xx[2:end-1], u[2:end-1];
end

function burgers_reconstruct(dx, u, compute_dt)
  # Calculate slopes via flux limiter
  theta = 1.5;
  t_forward = theta*(u[2:end-1] .- u[1:end-2]) / dx;
  central = (u[3:end] .- u[1:end-2]) / (2*dx)
  t_backward = theta*(u[3:end] .- u[2:end-1]) / dx;
  sigma = minmod.(t_forward, central, t_backward);

  # Build U_{j+1/2}^{\pm} for j=1,...,n and get U_{1/2}^-, U_{n+1/2}^+ from b.c.
  Um = zeros(length(u)-1); Up = zeros(length(u)-1);
  Um[1] = u[1]; # TODO: Verify!
  Um[2:end] = u[2:end-1] .+ (dx/2) * sigma;
  Up[1:end-1] = u[2:end-1] .- (dx/2) * sigma;
  Up[end] = u[end]; # TODO: Verify!

  # Compute local speeds a_{j+1/2}^\pm
  zero_blk = zeros(length(u)-1);
  am = 2*min.(Up, Um, zero_blk);
  ap = 2*max.(Up, Um, zero_blk);

  # Compute numerical fluxes F_{j+1/2}
  f(x) = x^2 / 2; #Burgers scalar conservation law
  denom = ap .- am;
  F = zeros(length(denom));
  for k=1:length(denom)
    if denom[k] > 1e-8
      F[k] = ((ap[k]*f(Um[k]) - am[k]*f(Up[k])) + (ap[k]*am[k])*(Up[k]-Um[k]));
      F[k] /= denom[k];
    else #Fallback to central difference if denom too small
      F[k] = (f(Um[k]) - f(Up[k])) / (2*dx);
    end
  end

  # Compute RHS of ODE system 
  rhs = -(F[2:end] .- F[1:end-1]) / dx;

  # if compute_dt estimate the time step
  dt = -1.0;
  if compute_dt
    amax = maximum(max.(ap, -am));
    dt = dx / (2*amax);
  end
  return rhs, dt
end

function minmod(a1, a2, a3)
  if a1 > 0 && a2 > 0 && a3 > 0
    return min(a1, a2, a3);
  elseif a1 < 0 && a2 < 0 && a3 < 0
    return max(a1, a2, a3);
  else
    return 0
  end
end
