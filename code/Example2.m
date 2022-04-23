clear; clc; clf;
%% Setting up the conditions we will need to implement
cell_count = 100; dx = 10/cell_count; x_range=[0,10]; T=15; M=5;
x=0.5*dx:dx:10-0.5*dx; z=0.5*dx/64:dx/64:10-0.5*dx/64;

steady_mean = 2-2*bottom(z);
steady_sd   = abs(bottom(z))/sqrt(3);
%% Building the matrix E
Phi=cell(M,1);
for j=1:M
   Phi{j}=legendre_poly(j-1);
end
E = build_E(Phi);

%% Building u0_h

% Our initial condition is u0 = 0 so coefficient functions should be
% uniformly 0 with the exception of the boundary condition
u0_h = zeros(M,cell_count+1);
% First column corresponds to cell 0 which reflects our boundary condition
% of u = 2. This boundary condition is constant in the random variable so
% E[2 Phi_k] = 2 E[Phi_k] = 2 E[Phi_1 Phi_k] = 2 dirac_{1k}, so only the
% first entry of the first column should be nonzero
u0_h(1,1) = 2;


%% Building b_h

% We now construct the bottom topography function. The bottom topography
% can be written as (2+z)*bottom(x) where bottom is the function defined
% below. 
coefficients=zeros(M,1);
for k=1:M
    coeff = conv([1 2], Phi{k});    % polynomial multiplication with z+2
    q=0.5*polyint(coeff);           % inner product
    coefficients(k)=diff(polyval(q,[-1 1]));
end
b = bottom(x);
b_h = [zeros(M,1) coefficients*b];  % We need to add a column for 
                                    % the left boundary condition 
                                    % which happens to be just zeros
                                            
%% Well Balanced Scheme

u_wb = burgers_wbsg(u0_h, b_h, E, T, dx);       % WB      eq 3.14
mean_wb = u_wb(1,2:end); 
sd_wb = u_wb.^2; sd_wb = sd_wb(2:end,2:end); sd_wb=sqrt(sum(sd_wb));

%% WB Plotting
figure(1)
plot(z,steady_mean,'k--','LineWidth', 1.5); hold on
plot(x,mean_wb,'o','Color',[0.6350 0.0780 0.1840],'LineWidth', 1.5);
ylim([1.75,2.25])
legend({'Steady State','Stochastic WB'},"Location","NorthEast",'FontSize',16);
ax = gca; ax.FontSize = 18; 

figure(2)
plot(z,steady_sd,'k--','LineWidth', 1.5); hold on
plot(x,sd_wb,'o','Color',[0.6350 0.0780 0.1840],'LineWidth', 1.5);
ylim([0,0.07])
legend({'Steady State','Stochastic WB'},"Location","NorthEast",'FontSize',16);
ax = gca; ax.FontSize = 18; 

%% Non-Well Balanced Scheme

u_non = burgers_non_wbsg(u0_h, b_h, E, T, dx);  % Non WB  eq 3.21
mean_non = u_non(1,2:end); sd_non = u_non.^2; 
sd_non = sd_non(2:end,2:end); sd_non=sqrt(sum(sd_non));

%% Non-WB Plotting

figure(3)
plot(z,steady_mean,'k--','LineWidth', 1.5); hold on
plot(x,mean_non,'^','Color',[0.9290 0.6940 0.1250],'LineWidth', 1.5);
ylim([1.75,2.25])
legend({'Steady State','gPC Non-WB'},"Location","NorthEast",'FontSize',16);
ax = gca; ax.FontSize = 18; 

figure(4)
plot(z,steady_sd,'k--','LineWidth', 1.5); hold on
plot(x,sd_non,'^','Color',[0.9290 0.6940 0.1250],'LineWidth', 1.5);
ylim([0,0.07])
legend({'Steady State','gPC Non-WB'},"Location","NorthEast",'FontSize',16);
ax = gca; ax.FontSize = 18; 


%% Functions Used

% Bottom topography function
function b = bottom(x)
    b=zeros(size(x)); N=length(x);
    for j=1:N
       if x(j) > 5 && x(j) < 6
          b(j) = 0.1*cos(pi*x(j)); 
       end
    end
end
