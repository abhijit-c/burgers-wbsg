using Polynomials
function build_E(Phi)
# build_E.jl constructs the 3D matrix where:
# \[
#   E(k,m,n) 
#   = Expectation[ Phi_k Phi_m Phi_n ] 
#   = \int_a^b \Phi_k(x) \Phi_m(x) \Phi_n(x) mu(x) ~dx
# \]
# Here, it's assumed that we're in a Unif(-1,1) setting, hence mu(x) = 0.5 and
# Phi_k are the orthonormal Legendre polynomials.
# ------------------------------------------------------------------------------
# Parameter List:
# Phi: Array of Polynomial objects where Phi[k] = Phi_k as descirbed above.
# Returned Variables:
# E: M x M x M matrix as described above.
# ------------------------------------------------------------------------------
# Authors: Andrew Shedlock & Abhijit Chowdhary, @ North Carolina State.
# Date  : 2022/04/14

  M = length(Phi);
  E = zeros(M,M,M);
  for k=1:M, m=1:M, n=1:M
    E[k,m,n] = integrate(0.5*Phi[k]*Phi[m]*Phi[n], -1, 1)
  end
  return E
end
