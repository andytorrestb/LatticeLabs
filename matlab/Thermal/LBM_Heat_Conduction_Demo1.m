clc;
clear all;

%% Define parameters

% Cylinder Geometry
% Cylinder radius will be used to determine the domain size (this must be
% an odd number to ensure easy sampling at the middle plane). 
radius=11;
D=2*radius;

% Domain Size and Resolution
L=5*radius;
N_x=L;
N_y=L;
dx=1;
dy=1;

% Cylinder center coordinates
x_circ=(N_x-1)/2;
y_circ=(N_y-1)/2;

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
T_1=0.5;
T_2=0.33;
T_3 = 1;

%% Zoning
%% Domain_ID=0 --- Solid Domain
%% Domain_ID=1 --- Fluid Domain
Domain_ID=zeros(N_y,N_x);
for j=1:N_y
    for i=1:N_x
        if test_circle(i-1,j-1,radius,x_circ,y_circ)
            Domain_ID(j,i)=0;
        else
            Domain_ID(j,i)=1;
        end
    end
end
% contourf(flipud(Domain_ID),30)
% title("Domain_ID")
% axis equal tight

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
% contourf(flipud(Zone_ID),30)
% title("Zone_ID")
% axis equal tight

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
                        T_w=T_1;
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
                    T_w=T_2;
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

for j=1:N_y
    for i=1:N_x
        % Manually change 0 values to 1
        if T(j,i) == 0
            T(j,i) = 1;
        end
    end
end
%% Visualization
figure;
contourf(flipud(T),30)
colormap(hot)
axis equal tight


% Establish arrays to store coordinate data for graphing.
for i=1:N_x
    x(i)=dx*(i-1);
end

for j=1:N_y
    y(j)=dy*(j-1);
end

% Normalize array data
x_norm = x / L;
y_norm = y / L;


% Sample and normalize temperature data.
T_x = T(x_circ, :);
T_y = T(:, y_circ);

T_x = (T_x - T_2) / (T_3 - T_2);
T_y = (T_y - T_2) / (T_3 - T_2);


% load benchmark data into bm
bm = load("Project3_Benchmark Data.mat")

% Calculate L2 error for each graph.
T_sim_interp_x = interp1(x_norm, T_x, bm.x_benchmark, 'linear');
L2_x = sqrt(sum((T_sim_interp_x - bm.T_benchmark_hori).^2));

T_sim_interp_y = interp1(y_norm, T_y, bm.y_benchmark, 'linear');
L2_y = sqrt(sum((T_sim_interp_y - bm.T_benchmark_vert).^2));

% Scatter plot
figure;
scatter(bm.x_benchmark, bm.T_benchmark_hori,'filled');  
hold on;

% Plot the smooth curve
plot(x_norm, T_x, 'r-', 'LineWidth', 2);

% Labels and title
xlabel('normlalized x-coordinates (x*)');
ylabel('Normalized Temperature (T*)');
title('Normalized Temperature Along the Horizontal Mid Plane (T* vs x*)');
legend('Data Points', 'Smooth Curve');
grid on;
% Add L2 norm as annotation
text(0.05, 0.95, sprintf('L_2 Norm: %.4f', L2_x), ...
    'Units', 'normalized', 'FontSize', 12, 'Color', 'k');

figure;
% Scatter plot
scatter(bm.y_benchmark, bm.T_benchmark_vert,'filled');  
hold on;

% Plot the smooth curve
plot(y_norm, T_y, 'r-', 'LineWidth', 2);

% Labels and title
xlabel('normlalized y-coordinates (y*)');
ylabel('Normalized Temperature (T*)');
title('Normalized Temperature Along the Vertical Mid Plane (T* vs y*)');
legend('Data Points', 'Smooth Curve');
grid on;

% Add L2 norm as annotation
text(0.05, 0.95, sprintf('L_2 Norm: %.4f', L2_y), ...
    'Units', 'normalized', 'FontSize', 12, 'Color', 'k');