function [u_h] = iburgers_wbsg(x_range, u0_h, ul, E, T, Nc)
% burgers_wbsg.m solves the Burgers equation 
% \[
%   u_t(x,t,z) + f_x(u(x,t,z)) = -b'(x,z)u(x,t,z), (*)
%   u(x,0) = u0(x),
%   u(0,t) = ul(t).
% \] 
% in space (x) and time (t) given some source term $-b'(x,z)u(t,x)$ where $z$ is
% some user specified random variable and $f(u) = u^2/2$. The random input is
% expanded with the generalized polynomial chaos (gPC) setting:
% \begin{align*}
%   u(x,t,z) &= u_N(x,t,z) = \sum_{m=1}^M \hat{u}_m(x,t)\Phi_m(z) \\
%   b(x,z)   &= b_N(x,z)   = \sum_{m=1}^M \hat{b}_m(x,t)\Phi_m(z) \\
% \end{align*}
% where $(\Phi_m)_{m=1}^M$ is a collection of orthogonal polynomials chosen with
% respect to the distribution zn obeys.
% ------------------------------------------------------------------------------
% Parameter List:
% x_range : [x_range(1),x_range(2)] denotes the spatial domain to solve (*)
%           over.
% u0_h : Anonymous function describing the initial condition of the model
%        represented in the basis of the orthogonal polynomials, i.e. 
%        \[
%          u0(x) = \sum_{m=1}^M u0_h(m) \Phi_m(x)
%        \]
%        Must be defined over x_range.
% ul : Anonymous function describing the (Dirichlet) boundary condition at
%      x_range(1).
% E : M x M x M matrix where $E(k,m,n) = Expectation[\Phi_k \Phi_m \Phi_n]$,
%     where $(\Phi_m)_{m=1}^M$ is the orthogonal polynomial basis descibed
%     earlier, and the expectation is with respect to the law z obeys.
% T : Final time to solve (*) up to.
% Nc : Number of cells to divide x_range into.
% Returned Variables:
% u_hat: Coefficients of u in its generalized polynomial chaos expansion in the
%        provided set of orthogonal polynomials
% ------------------------------------------------------------------------------
% Authors: Andrew Shedlock & Abhijit Chowdhary, @ North Carolina State.
% Date  : 2022/04/14
  f = @(u) u^2 / 2;
end

function A = build_A(E, u_h)
  M = length(E);
  for m=1:M
    for n=1:M
      A(m,n) = u_h' * E[:,m,n];
    end
  end
end

function B = build_B(E, b_h)
  M = length(E);
  for m=1:M
    for n=1:M
      b(m,n) = b_h' * E[:,m,n];
    end
  end
end
