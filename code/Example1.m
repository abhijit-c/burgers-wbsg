%% Setting up the conditions we will need to implement
cell_count = 200; dx = 10/cell_count; x_range=[0,10]; T=10;

% building the matrix E
Phi=cell(5,1);
Phi{1}=legendre_poly(0);
Phi{2}=legendre_poly(1);
Phi{3}=legendre_poly(2);
Phi{4}=legendre_poly(3);
Phi{5}=legendre_poly(4);
E = build_E(Phi); M=length(Phi);


% Our initial condition is u0 = 0 so coefficient functions should be
% uniformly 0 with the exception of the boundary condition
u0_h = zeros(M,cell_count+1);
% First column corresponds to cell 0 which reflects our boundary condition
% of u = 2. This boundary condition is constant in the random variable so
% E[2 Phi_k] = 2 E[Phi_k] = 2 E[Phi_1 Phi_k] = 2 dirac_{1k}, so only the
% first entry of the first column should be nonzero
u0_h(1,1) = 2;



% We now construct the bottom topography function. The bottom topography
% can be written as (2+z)*bottom(x) where bottom is the function defined
% below. 
coefficients=zeros(M,1);
for k=1:M
    coeff = conv([1 2], Phi{k});    % polynomial multiplication with z+2
    q=0.5*polyint(coeff);           % inner product
    coefficients(k)=diff(polyval(q,[-1 1]));
end

b_h = [zeros(M,1) coefficients*bottom(x)];  % We need to add a column for 
                                            % the left boundary condition 
                                            % which happens to be just zeros

u_h = burgers_wbsg(x_range, u0_h, b_h, E, T, dx);



%% Functions Used

%Has to be called one time for each instance of random variable
function b = bottom(x)
    b=zeros(size(x)); N=length(x);
    for j=1:N
       if x(j) > 4.5 && x(j) < 5.5
          b(j) = cos(pi*x(j)); 
       end
    end
end