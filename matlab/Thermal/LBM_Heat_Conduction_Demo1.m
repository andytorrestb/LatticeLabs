clc;
clear all;

%% Define parameters

% Cylinder Geometry
R=11; % Cylinder radius will be used to determine the domain size
D=2*R;

% Domain Size and Resolution
L=5*R;
N_x=L;
N_y=L;
dx=1;
dy=1;

% Cylinder center coordinates
x_circ=(N_x-1)/2;
y_circ=(N_y-1)/2;

% Establish arrays to store coordinate data for graphing.
for i=1:N_x
    x(i)=dx*(i-1);
end

for j=1:N_y
    y(j)=dy*(j-1);
end


% LBM related
Ksi=[0 1 0 -1  0 1 -1  -1  1;...
     0 0 1  0 -1 1  1  -1 -1];
w=[0 1/6 1/6 1/6 1/6 1/12 1/12 1/12 1/12];
c_s=1;
Tau=1;
%% Initialization

% Heat conduction parameters
R=8.314; % Gas constant
k_t = 0.125;
h = 1;
T_inf = 0.8;
T_H=0.33;
T_L=0.5;

%% Zoning
%% Domain_ID=0 --- Solid Domain
%% Domain_ID=1 --- Fluid Domain
Domain_ID=zeros(N_y,N_x);
for j=1:N_y
    for i=1:N_x
        if test_circle(i-1,j-1,R,x_circ,y_circ)
            Domain_ID(j,i)=0;
        else
            Domain_ID(j,i)=1;
        end
    end
end
contourf(flipud(Domain_ID),30)
title("Domain_ID")
axis equal tight

%% Zone_ID=0 --- Dead zone
%% Zone_ID=1 --- Boundary/Solid nodes
%% Zone_ID=2 --- Fluid Nodes
%% Zone_ID=3 --- Regular nodes in the fluid domain (including the nodes on the outer boundaries)
Zone_ID=zeros(N_y,N_x);
for j=1:N_y
    for i=1:N_x
        if Domain_ID(j,i)==0 % Solid Domain
            if Domain_ID(j,i+1)==1 || Domain_ID(j,i-1)==1 || Domain_ID(j-1,i)==1 || Domain_ID(j+1,i)==1 || Domain_ID(j-1,i+1)==1 || Domain_ID(j-1,i-1)==1 || Domain_ID(j+1,i+1)==1 || Domain_ID(j+1,i-1)==1
                Zone_ID(j,i)=1;
            else
                Zone_ID(j,i)=0;
            end
        else % Fluid Domain
            if j==1 || j==N_y || i==1 || i==N_x
                Zone_ID(j,i)=3;
            else
                if Domain_ID(j,i+1)==0 || Domain_ID(j,i-1)==0 || Domain_ID(j-1,i)==0 || Domain_ID(j+1,i)==0 || Domain_ID(j-1,i+1)==0 || Domain_ID(j-1,i-1)==0 || Domain_ID(j+1,i+1)==0 || Domain_ID(j+1,i-1)==0
                    Zone_ID(j,i)=2;
                else
                    Zone_ID(j,i)=3;
                end
            end
        end
    end
end
contourf(flipud(Zone_ID),30)
title("Zone_ID")
axis equal tight

% Energy fields
Rho_in=1;
Rho=ones(N_y,N_x)*Rho_in;
T=ones(N_y,N_x);

% Initialization of equalibrium PDF.
g_eq=zeros(N_y,N_x,9);
for j=1:N_y
    for i=1:N_x
        if Domain_ID(j,i) == 1
            for k=1:9
                g_eq(j,i,k)=w(k)*Rho(j,i)*R*T(j,i);
            end
        end
    end
end
g=g_eq;
g_new=g;

%% Main simulation loop.
Timer=10000;
for t=1:Timer
    disp(t)
    % Streaming
    for j=1:N_y
        for i=1:N_x
            if Zone_ID(j,i) == 0 % Dead Zone
                % Do Nothing
            elseif Zone_ID(j,i) == 1 % Boundary Nodes
                % Do Nothing
            elseif Zone_ID(j,i) == 2 % Fluid Nodes
            else 
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

