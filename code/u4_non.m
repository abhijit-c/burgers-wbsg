function u_h = u4_non(u0_h, b_h, E, T, dx)
  [M, N] = size(u0_h);
  B = zeros(M,M,N);       %B(:,:,j) is B_j defined by (3.15)
  for j=1:N
    B(:, :, j) = build_Bj(E, b_h(:, j));
  end

  u_h = u0_h;
  curr_t = 0.0;
  temp = zeros( size(u_h) );
  while curr_t < T
      % Updating the cell of each node
    Sj_prev =build_Sj(E, u_h(:, 1));
    for j=2:N
      Sj = build_Sj(E, u_h(:, j));
      dt = 0.0025/8;
      
      rhs = -( Sj*u_h(:,j) - Sj_prev*u_h(:,j-1) ) / (4*dx);
      rhs = rhs - ( B(:,:,j)-B(:,:,j-1) )*u_h(:,j)/dx;
      % Updating the j column
      temp(:,j) = u_h(:,j) + dt * rhs;
      Sj_prev = Sj;
    end
    u_h(:,2:end) = temp(:,2:end);   % first column which corresponds to 
                                    % cell 0 should remain constant
    curr_t = curr_t + dt;
  end
end

function Sj = build_Sj(E, u_hj)
  M = length(E);
  Sj = zeros(M, M);
  for p=1:M
    for q=p:M
        s=0;
        for l=1:M
           for m=1:M
              for n=1:M
                s = s + u_hj(l)*u_hj(m)*u_hj(n)*E(l,m,n,p,q);
              end
           end
        end
        Sj(p,q) = s;
        Sj(q,p) = s;
    end
  end
end

function Bj = build_Bj(E, b_hj)
  M = length(E);
  Bj = zeros(M, M);
  for m=1:M
    for n=1:M
      Bj(m,n) = b_hj' * E(:,m,n,1,1);
    end
  end
end




