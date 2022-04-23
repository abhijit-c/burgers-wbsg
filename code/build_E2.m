function E = build_E2(Phi)
% build_E2.m constructs the 5D matrix where:
% \[
%   E(k,m,n) 
%   = Expectation[ Phi_l Phi_m Phi_n Phi_p Phi_q] 
%   = 0.5*\int_-1^1 \Phi_l(x) \Phi_m(x) \Phi_n(x) \Phi_p(x) \Phi_q(x) dx
% \]
% ------------------------------------------------------------------------------
% Parameter List:
% Phi: Cell array with length(Phi) = M. Contains the polynomial coefficients of
%      the given orthogonal polynomials, so that Phi_n(x) = polyval(Phi{n},x)
% Returned Variables:
% E: M x M x M x M x M matrix as described above.
% ------------------------------------------------------------------------------
% Authors: Andrew Shedlock & Abhijit Chowdhary, @ North Carolina State.
% Date  : 2022/04/14

  M = length(Phi);
  E = zeros(M,M,M,M,M);
  for l=1:M %Quintuple loop, but M is small so it's fine.
    for m=1:M
      for n=1:M
        for p=1:M
           for q=1:M
            % Alternative approach using polynomial properties and matlab
            % functions
            x=Phi{l}; y=Phi{m}; z=Phi{n}; u=Phi{p}; v=Phi{q};
            coeff = conv(conv(conv(x,y),conv(u,v)),z); % polynomial multiplication
            poly=0.5*polyint(coeff);
            E(l,m,n,p,q) = diff(polyval(poly,[-1 1]));
        
            % Idea for this implementation came from
            % https://www.mathworks.com/help/matlab/ref/polyint.html
            % https://www.mathworks.com/help/matlab/ref/conv.html
           end
        end
      end
    end
  end
end

