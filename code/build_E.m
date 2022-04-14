function E = build_E(x_range, Phi, mu)
% build_E.m constructs the 3D matrix where:
% \[
%   E(k,m,n) 
%   = Expectation[ Phi_k Phi_m Phi_n ] 
%   = \int_a^b \Phi_k(x) \Phi_m(x) \Phi_n(x) mu(x) ~dx
% \]
% ------------------------------------------------------------------------------
% Parameter List:
% x_range: Range to integrate over, x_range = [a,b].
% Phi: Cell array with length(Phi) = M. Contains the polynomial coefficients of
%      the given orthogonal polynomials, so that Phi_n(x) = polyval(Phi{n},x).
% mu: Anonymous function such that $x \mapsto mu(x)$.
% Returned Variables:
% E: M x M x M matrix as described above.
% ------------------------------------------------------------------------------
% Authors: Andrew Shedlock & Abhijit Chowdhary, @ North Carolina State.
% Date  : 2022/04/14
  a = x_range(1); b = x_range(2);
  M = length(Phi)
  E = zeros(M,M,M);
  for k=1:M %Triple loop, but M is small so it's fine.
    for m=1:M
      for n=1:M
        f = @(x) polyval(Phi{k},x)*polyval(Phi{m},x)*polyval(Phi{n},x)*mu(x);
        E(k,m,n) = integrate(f, a, b);
      end
    end
  end
end
