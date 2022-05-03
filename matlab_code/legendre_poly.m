% legendre_poly.m by David Terr, Raytheon, 5-10-04
% Given nonnegative integer n, compute the 
% Legendre polynomial P_n. Return the result as a vector whose mth
% element is the coefficient of x^(n+1-m).
% polyval(legendre_poly(n),x) evaluates P_n(x).
function pk = legendre_poly(n)
  if n==0 
    pk = 1;
  elseif n==1
    pk = [sqrt(3) 0];
  else
    pkm2 = zeros(1,n+1);
    pkm2(n+1) = 1;
    pkm1 = zeros(1,n+1);
    pkm1(n) = 1;
    for k=2:n
      pk = zeros(1,n+1);
      for e=n-k+1:2:n
        pk(e) = (2*k-1)*pkm1(e+1) + (1-k)*pkm2(e);
      end
      pk(n+1) = pk(n+1) + (1-k)*pkm2(n+1);
      pk = pk/k;
      if k<n
        pkm2 = pkm1;
        pkm1 = pk;
      end
    end
    % Normalizing the eigenvectors
    coeff=conv(pk,pk);
    q=polyint(coeff);
    a=0.5*diff(polyval(q,[-1 1])); a=sqrt(a);
    pk=pk/a;
  end
