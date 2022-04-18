function u_h = burgers_non_wbsg(u0_h, b_h, E, T, dx)
% burgers_wbsg.m solves the Burgers equation 
% \[
%   u_t(x,t,z) + f_x(u(x,t,z)) = -b'(x,z)u(x,t,z), (*)
%   u(x,0) = u0(x),
%   u(0,t) = 2
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
% u0_h : M x N matrix where the jth column corresponds to the cell average of
%        the initial condition over the jth spatial cell expanded out in the
%        orthogonal polynomial basis corresponding to z, i.e.:
%        \[
%          avg(u0(x_j)) \approx \sum_{m=1}^M u0_h(m,j) polyval(Phi{m},x_j)
%        \]
% b_h : M x N matrix where the jth column corresponds to the cell average of the
%       initial condition over the jth spatial cell expanded out in the
%       orthogonal polynomial basis corresponding to z, i.e.:
%       \[
%         avg(b(x_j)) \approx \sum_{m=1}^M b_h(m,j) polyval(Phi{m},x_j)
%       \]
% E : M x M x M matrix where $E(k,m,n) = Expectation[\Phi_k \Phi_m \Phi_n]$,
%     where $(\Phi_m)_{m=1}^M$ is the orthogonal polynomial basis descibed
%     earlier, and the expectation is with respect to the law z obeys.
% T : Final time to solve (*) up to.
% dx : Parameter controlling width of spatial discretization.
% Returned Variables:
% u_hat: Coefficients of u in its generalized polynomial chaos expansion in the
%        provided set of orthogonal polynomials
% ------------------------------------------------------------------------------
% Authors: Andrew Shedlock & Abhijit Chowdhary, @ North Carolina State.
% Date  : 2022/04/14
  
  
  [M, N] = size(u0_h);
  B = zeros(M,M,N);       %B(:,:,j) is B_j defined by (3.15)
  for j=1:N
    B(:, :, j) = build_Bj(E, b_h(:, j));
  end

  u_h = u0_h;
  curr_t = 0.0;
  temp = zeros( size(u_h) );
  while curr_t < T
      % Updating the cell of each node
    for j=2:N
      Aj = build_Aj(E, u_h(:, j)); Aj_1 = build_Aj(E, u_h(:, j-1));

      %dt = dx^2 %TODO: Figure out how to choose this.
      % from the paper, they took dt = 0.0025/8
      dt = 0.0025/8;
      
      rhs = -( Aj*u_h(:,j) - Aj_1*u_h(:,j-1) ) / (2*dx);
      rhs = rhs - ( B(:,:,j)-B(:,:,j-1) )*( u_h(:,j) ) / (dx);
      % Updating the j column
      temp(:,j) = u_h(:,j) + dt * rhs;
    end
    u_h(:,2:end) = temp(:,2:end);   % first column which corresponds to 
                                    % cell 0 should remain constant
    curr_t = curr_t + dt;
  end
end

function Aj = build_Aj(E, u_hj)
  M = length(E);
  Aj = zeros(M, M);
  for m=1:M
    for n=1:M
      Aj(m,n) = u_hj' * E(:,m,n);
    end
  end
end

function Bj = build_Bj(E, b_hj)
  M = length(E);
  Bj = zeros(M, M);
  for m=1:M
    for n=1:M
      Bj(m,n) = b_hj' * E(:,m,n);
    end
  end
end

