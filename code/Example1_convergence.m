clc; clear; clf


%% Showing convergence for Example 1
N=1:5; cells=[100, 200, 400]; x_range=[0,10]; T=10;
mean_err = zeros(length(cells),length(N));
sd_err = zeros(length(cells),length(N));

% Building E tensor separately and will call subtensor as needed

Phi=cell(N(end),1);
for j=1:N(end)
   Phi{j}=legendre_poly(j-1);
end
E = build_E(Phi);


for n=1:N(end)
    En = E(1:n,1:n,1:n);
    for j=1:length(cells)
        cell_count=cells(j); dx = 10/cell_count;
        x=0.5*dx:dx:10-0.5*dx;
        
        u0_h = zeros(n,cell_count+1);
        u0_h(1,1) = 2;
       
       
        coefficients=zeros(n,1);
        for k=1:n
            coeff = conv([1 2], Phi{k});    
            q=0.5*polyint(coeff);           
            coefficients(k)=diff(polyval(q,[-1 1]));
        end
        b = bottom(x);
        b_h = [zeros(n,1) coefficients*b];  

        steady_mean = 2-2*b;
        steady_sd   = abs(b)/sqrt(3);
        
        u_wb = burgers_wbsg(u0_h, b_h, En, T, dx);       
        mean_wb = u_wb(1,2:end); 
        sd_wb = u_wb.^2; sd_wb = sd_wb(2:end,2:end); sd_wb=sqrt(sum(sd_wb));
        
        % Error analysis
        mean_err(j,n) = err(steady_mean,mean_wb);
        sd_err(j,n) = err(steady_sd,sd_wb);
   end
end
%% Functions used

% Bottom topography function
function b = bottom(x)
    b=zeros(size(x)); N=length(x);
    for j=1:N
       if x(j) > 4.5 && x(j) < 5.5
          b(j) = cos(pi*x(j)); 
       end
    end
end

% Computes the L1 error
function e = err(actual,approx)
    cell_count=length(actual);
    dx = 10/cell_count;
    e=dx*norm(actual-approx);
end