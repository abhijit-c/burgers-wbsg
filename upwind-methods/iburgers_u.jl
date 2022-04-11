"""
  iburgers_u(xx, u0, T, dx, dt)

Solve the inviscid Burgers equation
\\[ u_t + (u^2/2)_x = 0 (*) \\]
with some initial condition u0(x) up to time T on spatial domain [a,b]. We
employ a conservative upwind scheme to numerically approximate the solution.
...
# Arguments
xx N-element Vector{Float64}: Grid to solve PDE over.
u0 Function: Initial condition u0(xl,xr) returns the average of u0 over [xl,xr].
T Float64: Final time to solve PDE to.
dx Float64: Width of spatial discretization.
dt Float64: Width of temporal discretization.
...
"""
function iburgers_u(xx, u0, T, dx, dt)
  f(u) = u^2 / 2; df(u) = u; #Burgers conservation law & derivative
  function Flux(u,v) #Numerical flux function
    F = zeros(length(u));
    for k=1:length(F)
      F[k] = ( (f(u[k])-f(v[k]))/(u[k]-v[k])>=0 ) ? f(u[k]) : f(v[k]);
    end
    return F;
  end

  a = xx[1]; b = xx[end]; 
  xx = [a-dx; xx; b+dx]; #Adjust xx to include ghost points.
  cell_centers = xx[1:end-1] .+ dx/2;
  u = u0.(xx[1:end-1], xx[2:end]);
  t = 0.0;
  while t < T
    diff_Flux = Flux(u[2:end-1],u[3:end]) .- Flux(u[1:end-2],u[2:end-1]);
    u[2:end-1] = u[2:end-1] .- (dt/dx) .* diff_Flux;
    t = min(t + dt, T);
  end
  return cell_centers[2:end-1], u[2:end-1];
end
