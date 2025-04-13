clear all;

%% Define parameters
% Domain
N_x=51;
N_y=51;
dx=1;
dy=1;
% LBM related
Ksi=[0 1 0 -1  0 1 -1  -1  1;...
     0 0 1  0 -1 1  1  -1 -1];
w=[0 1/6 1/6 1/6 1/6 1/12 1/12 1/12 1/12];
c_s=1;
Tau=1;
%% Initialization
Rho_in=1;
Rho=ones(N_y,N_x)*Rho_in;
T=ones(N_y,N_x);
R=8.314; % Gas constant
k_t = 1;
h = 1;
T_inf = 0.8;
T_H=1;
T_L=0.1;
g_eq=zeros(N_y,N_x,9);
for j=1:N_y
    for i=1:N_x
        for k=1:9
            g_eq(j,i,k)=w(k)*Rho(j,i)*R*T(j,i);
        end
    end
end
g=g_eq;
g_new=g;

Timer=10000;
%% Solving
for t=1:Timer
    disp(t)
    % Streaming
    for j=1:N_y
        for i=1:N_x
            if j==1 % Top surface
                if i==1 % Top-Left corner
                    g_new(j,i,1)=g(j,i,1);
                    g_new(j,i,3)=g(j+1,i,3);
                    g_new(j,i,4)=g(j,i+1,4);
                    g_new(j,i,7)=g(j+1,i+1,7);

                    g_new(j,i,2)=g_new(j,i,4);
                    g_new(j,i,5)=g_new(j,i,3);
                    T_w=T(j+1,i);
                    % T_w=T_H; % second option
                    g_new(j,i,9)=g_new(j,i,7);
                    g_new(j,i,6)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,7)-g_new(j,i,9))/2;
                    g_new(j,i,8)=g_new(j,i,6);
                elseif i==N_x % Top-Right corner
                    g_new(j,i,1)=g(j,i,1);
                    g_new(j,i,2)=g(j,i-1,2);
                    g_new(j,i,3)=g(j+1,i,3);
                    g_new(j,i,6)=g(j+1,i-1,6);

                    g_new(j,i,4)=g_new(j,i,2);
                    g_new(j,i,5)=g_new(j,i,3);
                    g_new(j,i,8)=g_new(j,i,6);
                    T_w=T(j+1,i);
                    % T_w=T_L % Second option
                    g_new(j,i,7)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,6)-g_new(j,i,8))/2;
                    g_new(j,i,9)=g_new(j,i,7);
                else % All other nodes on the top surface
                    g_new(j,i,1)=g(j,i,1);
                    g_new(j,i,2)=g(j,i-1,2);
                    g_new(j,i,3)=g(j+1,i,3);
                    g_new(j,i,4)=g(j,i+1,4);
                    g_new(j,i,6)=g(j+1,i-1,6);
                    g_new(j,i,7)=g(j+1,i+1,7);

                    g_new(j,i,5)=g_new(j,i,3);
                    T_w=T_L;
                    g_new(j,i,8)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,6)-g_new(j,i,7))/2;
                    g_new(j,i,9)=g_new(j,i,8);
                end
            elseif j==N_y % Bottom surface
                if i==1 % Bottom-Left corner
                    g_new(j,i,1)=g(j,i,1);
                    g_new(j,i,4)=g(j,i+1,4);
                    g_new(j,i,5)=g(j-1,i,5);
                    g_new(j,i,8)=g(j-1,i+1,8);

                    g_new(j,i,2)=g_new(j,i,4);
                    g_new(j,i,3)=g_new(j,i,5);
                    g_new(j,i,6)=g_new(j,i,8);
                    T_w=T(j-1,1);
                    % T_w=T_H; % Second option
                    g_new(j,i,7)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,6)-g_new(j,i,8))/2;
                    g_new(j,i,9)=g_new(j,i,7);
                elseif i==N_x % Bottom-Right corner
                    g_new(j,i,1)=g(j,i,1);
                    g_new(j,i,2)=g(j,i-1,2);
                    g_new(j,i,5)=g(j-1,i,5);
                    g_new(j,i,9)=g(j-1,i-1,9);

                    g_new(j,i,3)=g_new(j,i,5);
                    g_new(j,i,4)=g_new(j,i,2); 
                    g_new(j,i,7)=g_new(j,i,9);
                    T_w=T(j-1,i);
                    % T_w=T_L; % Second option
                    g_new(j,i,6)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,7)-g_new(j,i,9))/2;
                    g_new(j,i,8)=g_new(j,i,6);
                else % All other nodes on the bottom surface
                    g_new(j,i,1)=g(j,i,1);
                    g_new(j,i,2)=g(j,i-1,2);
                    g_new(j,i,4)=g(j,i+1,4);
                    g_new(j,i,5)=g(j-1,i,5);
                    g_new(j,i,8)=g(j-1,i+1,8);
                    g_new(j,i,9)=g(j-1,i-1,9);

                    g_new(j,i,3)=g_new(j,i,5);
                    T_w = (k_t/dy * T(j-1,i) + h * T_inf) / (k_t/dy + h);
                    g_new(j,i,6)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,8)-g_new(j,i,9))/2;
                    g_new(j,i,7)=g_new(j,i,6);
                end
            elseif i==1 % Left surface
                g_new(j,i,1)=g(j,i,1);
                g_new(j,i,3)=g(j+1,i,3);
                g_new(j,i,4)=g(j,i+1,4);
                g_new(j,i,5)=g(j-1,i,5);
                g_new(j,i,7)=g(j+1,i+1,7);
                g_new(j,i,8)=g(j-1,i+1,8);

                g_new(j,i,2)=g_new(j,i,4);
                T_w=T_H;
                g_new(j,i,6)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,7)-g_new(j,i,8))/2;
                g_new(j,i,9)=g_new(j,i,6);
            elseif i==N_x % Right surface
                g_new(j,i,1)=g(j,i,1);
                g_new(j,i,2)=g(j,i-1,2);
                g_new(j,i,3)=g(j+1,i,3);
                g_new(j,i,5)=g(j-1,i,5);
                g_new(j,i,6)=g(j+1,i-1,6);
                g_new(j,i,9)=g(j-1,i-1,9);

                g_new(j,i,4)=g_new(j,i,2);
                T_w=T(j, i-1);
                g_new(j,i,7)=(Rho(j,i)*R*T_w-g_new(j,i,1)-g_new(j,i,2)-g_new(j,i,3)-g_new(j,i,4)-g_new(j,i,5)-g_new(j,i,6)-g_new(j,i,9))/2;
                g_new(j,i,8)=g_new(j,i,7);
            else % All interior nodes
                g_new(j,i,1)=g(j,i,1);
                g_new(j,i,2)=g(j,i-1,2);
                g_new(j,i,3)=g(j+1,i,3);
                g_new(j,i,4)=g(j,i+1,4);
                g_new(j,i,5)=g(j-1,i,5);
                g_new(j,i,6)=g(j+1,i-1,6);
                g_new(j,i,7)=g(j+1,i+1,7);
                g_new(j,i,8)=g(j-1,i+1,8);
                g_new(j,i,9)=g(j-1,i-1,9);
            end
        end
    end
    % Collsion
    % Computing moments
    for j=1:N_y
        for i=1:N_x
            T(j,i)=sum(g_new(j,i,:))/R/Rho(j,i);
        end
    end
    % Computing f_eq
    for j=1:N_y
        for i=1:N_x
            g_eq(j,i,:)=w'*Rho(j,i)*R*T(j,i);
        end
    end
    % Collision & Update
    g=g_new-1/Tau*(g_new-g_eq);
end

%% Visualization
% Analytical Solution
x_benchmark=0:0.01:1;
for i=1:length(x_benchmark)
    T_benchmark(i)=-1*(x_benchmark(i)-1);
end
figure;
plot(x_benchmark,T_benchmark,"red")
% Simultion result
for i=1:N_x
    T_sim(i)=T((N_y-1)/2+1,i);
end
hold on
plot((0:1:N_x-1)/(N_x-1),(T_sim-T_L)/(T_H-T_L),'blue');

figure
contourf(flipud(T),30)
axis equal tight

