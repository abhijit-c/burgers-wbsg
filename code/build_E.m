function E = build_E(Phi)
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

  M = length(Phi);
  E = zeros(M,M,M);
  for k=1:M %Triple loop, but M is small so it's fine.
    for m=1:M
      for n=1:M
%           
%         f = @(x) polyval(Phi{k},x)*polyval(Phi{m},x)*polyval(Phi{n},x)*mu(x);
%         E(k,m,n) = integrate(f, a, b);
            
        % Alternative approach using polynomial properties and matlab
        % functions
        u=Phi{k}; v=Phi{m}; w=Phi{n};
        coeff = conv(conv(u,v),w); 
        q=polyint(coeff);
        E(k,m,n) = 0.5*diff(polyval(q,[-1 1]));
        
        % Idea for this implementation came from
        % https://www.mathworks.com/help/matlab/ref/polyint.html
        % https://www.mathworks.com/help/matlab/ref/conv.html
      end
    end
  end
end
