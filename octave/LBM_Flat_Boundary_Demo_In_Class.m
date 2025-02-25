clear all;
%% Define parameters
% Domain
N_x=50;
N_y=50;
dx=1;
dy=1;
% LBM related
Ksi=[0 1 0 -1  0 1 -1  -1  1;...
     0 0 1  0 -1 1  1  -1 -1];
w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
c_s=1/sqrt(3);
Tau=0.58655;
%% Initialization
Rho_in=1;
U_x=0.1*c_s;
Rho=ones(N_y,N_x)*Rho_in;
u=zeros(N_y,N_x);
v=zeros(N_y,N_x);
f_eq=zeros(N_y,N_x,9);
for j=1:N_y
    for i=1:N_x
        for k=1:9
            f_eq(j,i,k)=w(k)*Rho(j,i)*(1+...
                                       Ksi(:,k)'*[u(j,i);v(j,i)]/c_s^2+...
                                       (Ksi(:,k)'*[u(j,i);v(j,i)])^2/(2*c_s^4)-...
                                       [u(j,i),v(j,i)]*[u(j,i);v(j,i)]/(2*c_s^2));
            f(j,i,k)=Rho_in*w(k);
        end
    end
end
f_new=f;
T=100;
%% Solving
for t=1:T
    disp(t)
    % Streaming
    for j=1:N_y
        for i=1:N_x
            if j==1 % Top surface
                if i==1 % Top-Left corner
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,3)=f(j+1,i,3);
                    f_new(j,i,4)=f(j,i+1,4);
                    f_new(j,i,7)=f(j+1,i+1,7);

                    Rho_ji = (Rho(j+1,i) + Rho(j,i+1) + Rho(j+1,i+1)) / 3.0;
                    f_new(j,i,2)=f_new(j,i,4);
                    f_new(j,i,5)=f_new(j,i,3);
                    f_new(j,i,6)=0.5*(Rho_ji - f_new(j,i,1) - 2*(f_new(j,i,3)+f_new(j,i,4)+f_new(j,i,7)));
                    f_new(j,i,8)=f_new(j,i,6);
                    f_new(j,i,9)=f_new(j,i,7);
                elseif i==N_x % Top-Right corner
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,2)=f(j,i-1,2);
                    f_new(j,i,3)=f(j+1,i,3);
                    f_new(j,i,6)=f(j+1,i-1,6);

                    Rho_ji = (Rho(j+1,i) + Rho(j,i-1) + Rho(j+1,i-1)) / 3.0;
                    f_new(j,i,4)=f_new(j,i,2);
                    f_new(j,i,5)=f_new(j,i,3);
                    f_new(j,i,7)=0.5*(Rho_ji - f_new(j,i,1) - 2*(f_new(j,i,2)+f_new(j,i,3)+f_new(j,i,6)));
                    f_new(j,i,8)=f_new(j,i,6);
                    f_new(j,i,9)=f_new(j,i,7);
                else % All other nodes on the top surface
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,2)=f(j,i-1,2);
                    f_new(j,i,3)=f(j+1,i,3);
                    f_new(j,i,4)=f(j,i+1,4);
                    f_new(j,i,6)=f(j+1,i-1,6);
                    f_new(j,i,7)=f(j+1,i+1,7);

                    % f_new(j,i,5)=f_new(j,i,3);
                    % f_new(j,i,8)=0.5*(f_new(j,i,1)+f_new(j,i,3)+2*(f_new(j,i,3)+f_new(j,i,4)+f_new(j,i,7))-Rho_in*(U_x-1));
                    % f_new(j,i,9)=0.5*(Rho_in*(U_x+1)+f_new(j,i,7)-f_new(j,i,1) -2*(f_new(j,i,2)+f_new(j,i,3)+f_new(j,i,6)));
                    f_new(j,i,5)=f_new(j,i,3);
                    f_new(j,i,8)=f_new(j,i,6) + 0.5*(-1*Rho(j,i)*U_x + f_new(j,i,2) - f_new(j,i,4));
                    f_new(j,i,9)=f_new(j,i,7) + 0.5*(Rho(j,i)*U_x + f_new(j,i,4) - f_new(j,i,2));
                end
            elseif j==N_y % Bottom surface
                if i==1 % Bottom-Left corner
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,4)=f(j,i+1,4);
                    f_new(j,i,5)=f(j-1,i,5);
                    f_new(j,i,8)=f(j-1,i+1,8);

                    Rho_ji = (Rho(j-1,i) + Rho(j,i+1) + Rho(j-1,i+1)) / 3.0;
                    f_new(j,i,2)=f_new(j,i,4);
                    f_new(j,i,3)=f_new(j,i,5);
                    f_new(j,i,6)=f_new(j,i,8);
                    f_new(j,i,7)=0.5*(Rho_ji - f_new(j,i,1) - 2*(f_new(j,i,4)+f_new(j,i,5)+f_new(j,i,8)));
                    f_new(j,i,9)=f_new(j,i,7);
                elseif i==N_x % Bottom-Right corner
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,2)=f(j,i-1,2);
                    f_new(j,i,5)=f(j-1,i,5);
                    f_new(j,i,9)=f(j-1,i-1,9);

                    Rho_ji = (Rho(j-1,i) + Rho(j,i-1) + Rho(j-1,i-1)) / 3.0;
                    f_new(j,i,3)=f_new(j,i,5);
                    f_new(j,i,4)=f_new(j,i,2);
                    f_new(j,i,6)=0.5*(Rho_in - f_new(j,i,1) - 2*(f_new(j,i,2)+f_new(j,i,5)+f_new(j,i,9)));
                    f_new(j,i,7)=f_new(j,i,9);
                    f_new(j,i,8)=f_new(j,i,6);
                else % All other nodes on the bottom surface
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,2)=f(j,i-1,2);
                    f_new(j,i,4)=f(j,i+1,4);
                    f_new(j,i,5)=f(j-1,i,5);
                    f_new(j,i,8)=f(j-1,i+1,8);
                    f_new(j,i,9)=f(j-1,i-1,9);

                    f_new(j,i,3)=f_new(j,i,5);
                    f_new(j,i,6)=(f_new(j,i,4)-f_new(j,i,2))/2+f_new(j,i,8);
                    f_new(j,i,7)=(f_new(j,i,2)-f_new(j,i,4))/2+f_new(j,i,9);
                end
            elseif i==1 % Left surface
                f_new(j,i,1)=f(j,i,1);
                f_new(j,i,3)=f(j+1,i,3);
                f_new(j,i,4)=f(j,i+1,4);
                f_new(j,i,5)=f(j-1,i,5);
                f_new(j,i,7)=f(j+1,i+1,7);
                f_new(j,i,8)=f(j-1,i+1,8);

                f_new(j,i,2)=f_new(j,i,4);
                f_new(j,i,6)=0.5*(f_new(j,i,5) - f_new(j,i,3) + 2*f_new(j,i,8));
                f_new(j,i,9)=0.5*(f_new(j,i,3) - f_new(j,i,5) + 2*f_new(j,i,7));
            elseif i==N_x % Right surface
                f_new(j,i,1)=f(j,i,1);
                f_new(j,i,2)=f(j,i-1,2);
                f_new(j,i,3)=f(j+1,i,3);
                f_new(j,i,5)=f(j-1,i,5);
                f_new(j,i,6)=f(j+1,i-1,6);
                f_new(j,i,9)=f(j-1,i-1,9);

                f_new(j,i,4)=f_new(j,i,2);
                f_new(j,i,7)=0.5*(f_new(j,i,5) - f_new(j,i,3) + 2*f_new(j,i,9));
                f_new(j,i,8)=0.5*(f_new(j,i,3) - f_new(j,i,3) + 2*f_new(j,i,6));
            else % All interior nodes
                f_new(j,i,1)=f(j,i,1);
                f_new(j,i,2)=f(j,i-1,2);
                f_new(j,i,3)=f(j+1,i,3);
                f_new(j,i,4)=f(j,i+1,4);
                f_new(j,i,5)=f(j-1,i,5);
                f_new(j,i,6)=f(j+1,i-1,6);
                f_new(j,i,7)=f(j+1,i+1,7);
                f_new(j,i,8)=f(j-1,i+1,8);
                f_new(j,i,9)=f(j-1,i-1,9);
            end
        end
    end
    % Collsion
    % Computing moments
    for j=1:N_y
        for i=1:N_x
            Rho(j,i)=sum(f_new(j,i,:));
            u(j,i)=(f_new(j,i,2)+f_new(j,i,6)+f_new(j,i,9)-f_new(j,i,4)-f_new(j,i,7)-f_new(j,i,8))/Rho(j,i);
            v(j,i)=(f_new(j,i,3)+f_new(j,i,6)+f_new(j,i,7)-f_new(j,i,5)-f_new(j,i,8)-f_new(j,i,9))/Rho(j,i);
        end
    end
    % Computing f_eq
    for j=1:N_y
        for i=1:N_x
            f_eq(j,i,:)=eqm_d2q9(Rho(j,i),[u(j,i);v(j,i)]);
        end
    end
    % Collision & Update
    f=f_new-1/Tau*(f_new-f_eq);
end


% Visualization

% Load experimental data
% Define y* values
y_star = [0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1];

% Define u* values
u_star = [0, -0.03717, -0.04192, -0.04775, -0.06434, -0.1015, -0.15662, -0.2109, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1];


% Define x* values
x_star = [0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5, ...
         0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1];
% Define the v* values
v_star = [0, 0.09233, 0.10091, 0.1089, 0.12317, 0.16077, 0.17507, 0.17527, ...
         0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0];


% Sample simulation data
u_min = min(u);
v_min = min(v);
u_range = range(u);
v_range = range(v);

u_norm_sim = zeros(N_y, N_x);
v_norm_sim = zeros(N_y, N_x);

for j=1:N_y
  for i=1:N_x

    u_norm_sim(j,i) = (u(j,i)-u_min)/u_range;
    v_norm_sim(j,i) = (v(j,i) - v_min)/v_range;

  endfor
endfor

%u_norm_sim = (u-u_min)/u_range
%v_norm_sim = (v-v_min)/v_range

disp('u_norm_sim')
disp(size(u_norm_sim(25,:)))
disp(u_norm_sim(:, 25))
%disp(u_norm_sim)

disp('u')
disp(size(u))



x_star_sim = linspace(0, 1, N_x);
y_star_sim = linspace(0, 1, N_y);


% Graphs the comparison
figure
title("v* vs x*")
ylabel("v*")
xlabel("x*")
scatter(x_star, v_star, mkr=".")
plot(x_star_sim, v_norm_sim(25, :))

figure
title("u* vs y*")
ylabel("y*")
xlabel("x*")
scatter(u_star,y_star,mkr=".")
plot(u_norm_sim(25, :), y_star_sim)



figure
quiver(flipud(u),flipud(v),10)
axis equal tight

figure
contourf(Rho,30)
axis equal tight
