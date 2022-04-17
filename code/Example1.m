%% Initial condition, boundary condition and bottom function 
left = 2;
u0 = @(x) 0;
b = @(x,z) bottom(x,z);


%% Setting up the conditions we will need to implement
dx=0.05; x=0:dx:10; x_range=[0,10];

Phi=cell(3,1);
Phi{1}=legendre_poly(0);
Phi{2}=legendre_poly(1);
Phi{3}=legendre_poly(2);
Phi{4}=legendre_poly(3);
E = build_E(Phi)

%u_h = burgers_wbsg(x_range, u0_h, b_h, left, E, T, dx);























%% Functions Used

%Has to be called one time for each instance of random variable
function b = bottom(x,z)
    b=zeros(size(x)); N=length(x);
    for j=1:N
       if x(j) > 4.5 && x(j) < 5.5
          b(j) = (2+z)*cos(pi*x(j)); 
       end
    end
end