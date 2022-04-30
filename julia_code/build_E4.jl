using Polynomials
function build_E4(Phi)
# build_E4.m constructs the 5D matrix where:
# \[
#   E(k,m,n) 
#   = Expectation[ Phi_l Phi_m Phi_n Phi_p Phi_q] 
#   = 0.5*\int_-1^1 \Phi_l(x) \Phi_m(x) \Phi_n(x) \Phi_p(x) \Phi_q(x) dx
# \]
# ------------------------------------------------------------------------------
# Parameter List:
# Phi: Cell array with length(Phi) = M. Contains the polynomial coefficients of
#      the given orthogonal polynomials, so that Phi_n(x) = polyval(Phi{n},x)
# Returned Variables:
# E: M x M x M x M x M matrix as described above.
# ------------------------------------------------------------------------------
# Authors: Andrew Shedlock & Abhijit Chowdhary, @ North Carolina State.
# Date  : 2022/04/14

  M = length(Phi);
  E = zeros(M,M,M,M,M);
  for l=1:M, m=1:M, n=1:M, p=1:M, q=1:M
    E[l,m,n,p,q] = integrate(0.5*Phi[l]*Phi[m]*Phi[n]*Phi[p]*Phi[q], -1, 1);
  end
  return E
end
