clc;
clear all;

%% Define parameters
% Domain
R=10;
D=2*R;
x_circ=5*D;
y_circ=4.5*D;
N_x=40*D;
N_y=9*D;
dx=1;
dy=1;
for i=1:N_x
    x(i)=dx*(i-1);
end

for j=1:N_y
    y(j)=dy*(j-1);
end
% LBM related
Ksi=[0 1 0 -1  0 1 -1  -1  1;...
     0 0 1  0 -1 1  1  -1 -1];
w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
c_s=1/sqrt(3);
U_in=0.1*c_s;
Re=20;
Tau=0.5+(U_in*D)/(Re*c_s*c_s);

C = [0 1 0 -1  0 1 -1 -1  1;
     0 0 1  0 -1 1  1 -1 -1];
Cx = reshape(C(1, :), [1, 1, 9]);
Cy = reshape(C(2, :), [1, 1, 9]);
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
%% Initialization
Rho_in=1;
Rho=ones(N_y,N_x)*Rho_in;
u=zeros(N_y,N_x);
v=zeros(N_y,N_x);
f_eq=zeros(N_y,N_x,9);
for i=1:N_x
    for j=1:N_y
        if Domain_ID(j,i)==1
            for k=1:9
                f_eq(j,i,k)=w(k)*Rho_in*(1+Ksi(:,k)'*[0.01;0]/c_s^2+(Ksi(:,k)'*[0.01;0])^2/c_s^4/2-[0.01,0]*[0.01;0]/c_s^2/2);
            end
        else
            for k=1:9
                f_eq(j,i,k)=w(k)*Rho_in*(1+Ksi(:,k)'*[0;0]/c_s^2+(Ksi(:,k)'*[0;0])^2/c_s^4/2-[0,0]*[0;0]/c_s^2/2);
            end
        end
    end
end
f=f_eq;
f_new=f;

T=20;
u_t = zeros(N_y, N_x, T);
v_t = zeros(N_y, N_x, T);
%% Solving
for t=1:T
    disp(t)
    % Streaming
    for j=1:N_y
        for i=1:N_x
            if Zone_ID(j,i)==0 % Dead zone
                % Do Nothing;
            elseif Zone_ID(j,i)==1 % Boundary/Solid nodes
                % Do Nothing;
            elseif Zone_ID(j,i)==2 % Fluid node
                % Implemet the curved boundary condition
                %% Direction 1
                f_new(j,i,1)=f(j,i,1);

                %% Direction 2
                if Zone_ID(j,i+1)==1
                    x1=x(i+1);
                    y1=y(j);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j, i-1);
                    v_ff = v(j, i-1);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node
                    % velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,2) - f_eq(j,i,2);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,4)=f_eq(j,i,2) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,4)=f(j,i+1,4);
                end

                %% Direction 3
                if Zone_ID(j-1,i)==1
                    x1=x(i);
                    y1=y(j-1);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j+1, i);
                    v_ff = v(j+1, i);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node
                    % velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,3) - f_eq(j,i,3);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,5)=f_eq(j,i,3) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,5)=f(j-1,i,5);
                end

                %% Direction 4
                if Zone_ID(j,i-1)==1
                    x1=x(i-1);
                    y1=y(j);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j, i+1);
                    v_ff = v(j, i+1);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node
                    % velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,4) - f_eq(j,i,4);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,2)=f_eq(j,i,4) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,2)=f(j,i-1,2);
                end

                %% Direction 5
                if Zone_ID(j+1,i)==1
                    x1=x(i);
                    y1=y(j+1);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j-1, i);
                    v_ff = v(j-1, i);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node
                    % velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,5) - f_eq(j,i,5);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,3)=f_eq(j,i,5) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,3)=f(j+1,i,3);
                end

                %% Direction 6
                if Zone_ID(j-1,i+1)==1
                    x1=x(i+1);
                    y1=y(j-1);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j+1, i-1);
                    v_ff = v(j+1, i-1);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node
                    % velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,6) - f_eq(j,i,6);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,8)=f_eq(j,i,6) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,8)=f(j-1,i+1,8);
                end

                %% Direction 7
                if Zone_ID(j-1,i-1)==1
                    x1=x(i-1);
                    y1=y(j-1);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j+1, i+1);
                    v_ff = v(j+1, i+1);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,7) - f_eq(j,i,7);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,9)=f_eq(j,i,7) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,9)=f(j-1,i-1,9);
                end

                %% Direction 8
                if Zone_ID(j+1,i-1)==1
                    x1=x(i-1);
                    y1=y(j+1);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j-1, i+1);
                    v_ff = v(j-1, i+1);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node
                    % velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,8) - f_eq(j,i,8);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,6)=f_eq(j,i,8) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,6)=f(j+1,i-1,6);
                end

                %% Direction 9
                if Zone_ID(j+1,i+1)==1
                    x1=x(i+1);
                    y1=y(j+1);
                    x2=x(i);
                    y2=y(j);
                    C_w=find_the_wall_point(x1,y1,x2,y2,R,x_circ,y_circ);
                    delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                    % Grab velocity data from the fluid node.
                    u_f = u(j,i);
                    v_f = v(j,i);

                    % Grab data from one node past the fluid node. 
                    u_ff = u(j-1, i-1);
                    v_ff = v(j-1, i-1);

                    % Wall node is assumed to not move.
                    u_w = 0;
                    v_w = 0;
                    
                    % Calculate the components of the boundary node
                    % velocity.
                    u_b1 = ((delta - 1)*u_f+u_w)/delta;
                    u_b2 = ((delta - 1)*u_ff+2*u_w)/(1+delta);
    
                    v_b1 = ((delta - 1)*v_f+v_w)/delta;
                    v_b2 = ((delta - 1)*v_ff+2*v_w)/(1+delta);

                    % Calculate the boundary node velocity using GZS scheme.
                    f_non_eq = f_new(j,i,9) - f_eq(j,i,9);
                    if delta>=0.75
                        u_ff=u_b1;
                        v_ff=v_b1;
                    else
                        u_ff=delta*u_b1+(1-delta)*u_b2;
                        v_ff=delta*v_b1+(1-delta)*v_b2;
                        f_non_eq = delta*f_non_eq+(1-delta)*f_non_eq;
                    end
                    f_new(j,i,7)=f_eq(j,i,9) + (1 - (1/Tau))*f_non_eq;
                else
                    f_new(j,i,7)=f(j+1,i+1,7);
                end

            else % Zone_ID=3
                if j==1 % Top surface
                    if i==1 % Top-Left corner
                        f_new(j,i,1)=f(j,i,1);
                        f_new(j,i,3)=f(j+1,i,3);
                        f_new(j,i,4)=f(j,i+1,4);
                        f_new(j,i,7)=f(j+1,i+1,7);

                        f_new(j,i,2)=f_new(j,i,4);
                        f_new(j,i,5)=f_new(j,i,3);
                        f_new(j,i,6)=Rho_in/2-f_new(j,i,1)/2-f_new(j,i,3)-f_new(j,i,4)-f_new(j,i,7);
                        f_new(j,i,8)=f_new(j,i,6);
                        f_new(j,i,9)=f_new(j,i,7);
                    elseif i==N_x % Top-Right corner
                        f_new(j,i,1)=f(j,i,1);
                        f_new(j,i,2)=f(j,i-1,2);
                        f_new(j,i,3)=f(j+1,i,3);
                        f_new(j,i,6)=f(j+1,i-1,6);

                        f_new(j,i,4)=f_new(j,i,2);
                        f_new(j,i,5)=f_new(j,i,3);
                        f_new(j,i,7)=f_new(j,i-1,7);
                        f_new(j,i,8)=f_new(j,i,6);
                        f_new(j,i,9)=f_new(j,i,7);
                    else % All other nodes on the top surface
                        f_new(j,i,1)=f(j,i,1);
                        f_new(j,i,2)=f(j,i-1,2);
                        f_new(j,i,3)=f(j+1,i,3);
                        f_new(j,i,4)=f(j,i+1,4);
                        f_new(j,i,6)=f(j+1,i-1,6);
                        f_new(j,i,7)=f(j+1,i+1,7);

                        % PDF entering the top boundary must be equal to
                        % the PDF leaving the bottom boundary.
                        f_new(j,i,9)=f(N_y,i,9);
                        f_new(j,i,5)=f(N_y,i,5);
                        f_new(j,i,8)=f(N_y,i,8);                    
                    end
                elseif j==N_y % Bottom surface
                    if i==1 % Bottom-Left corner
                        f_new(j,i,1)=f(j,i,1);
                        f_new(j,i,4)=f(j,i+1,4);
                        f_new(j,i,5)=f(j-1,i,5);
                        f_new(j,i,8)=f(j-1,i+1,8);

                        f_new(j,i,2)=f_new(j,i,4);
                        f_new(j,i,3)=f_new(j,i,5);
                        f_new(j,i,6)=f_new(j,i,8);
                        f_new(j,i,7)=Rho_in/2-f_new(j,i,1)/2-f_new(j,i,4)-f_new(j,i,5)-f_new(j,i,8);
                        f_new(j,i,9)=f_new(j,i,7);
                    elseif i==N_x % Bottom-Right corner
                        f_new(j,i,1)=f(j,i,1);
                        f_new(j,i,2)=f(j,i-1,2);
                        f_new(j,i,5)=f(j-1,i,5);
                        f_new(j,i,9)=f(j-1,i-1,9);

                        f_new(j,i,3)=f_new(j,i,5);
                        f_new(j,i,4)=f_new(j,i,2);
                        f_new(j,i,6)=f_new(j,i-1,6);
                        f_new(j,i,7)=f_new(j,i,9);
                        f_new(j,i,8)=f_new(j,i,6);
                    else % All other nodes on the bottom surface
                        f_new(j,i,1)=f(j,i,1);
                        f_new(j,i,2)=f(j,i-1,2);
                        f_new(j,i,4)=f(j,i+1,4);
                        f_new(j,i,5)=f(j-1,i,5);
                        f_new(j,i,8)=f(j-1,i+1,8);
                        f_new(j,i,9)=f(j-1,i-1,9);

                        % PDF enetering the bottom surace must be equal to
                        % the PDF leaving the top surface. 
                        f_new(j,i,6)=f(1,i-1,6);
                        f_new(j,i,3)=f(1,i,3);
                        f_new(j,i,7)=f(1,i+1,7);
                    end
                elseif i==1 % Left surface
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,3)=f(j+1,i,3);
                    f_new(j,i,4)=f(j,i+1,4);
                    f_new(j,i,5)=f(j-1,i,5);
                    f_new(j,i,7)=f(j+1,i+1,7);
                    f_new(j,i,8)=f(j-1,i+1,8);

                    Rho_in = (f_new(j,i,1)+f_new(j,i,3)+f_new(j,i,5)+2*(f_new(j,i,4)+f_new(j,i,7)+f_new(j,i,8)))/(1-1.667*U_in);
                    f_new(j,i,2)=f_new(j,i,4)+U_in*Rho_in*2/3;
                    f_new(j,i,6)=0.5*(Rho_in*U_in+2*f_new(j,i,8)+f_new(j,i,5)-f_new(j,i,3));
                    f_new(j,i,9)=0.5*(Rho_in*U_in+2*f_new(j,i,7)+f_new(j,i,3)-f_new(j,i,5));
                elseif i==N_x % Right surface
                    f_new(j,i,1)=f(j,i,1);
                    f_new(j,i,2)=f(j,i-1,2);
                    f_new(j,i,3)=f(j+1,i,3);
                    f_new(j,i,5)=f(j-1,i,5);
                    f_new(j,i,6)=f(j+1,i-1,6);
                    f_new(j,i,9)=f(j-1,i-1,9);

                    f_new(j,i,4)=f_new(j,i-1,4);
                    f_new(j,i,7)=f_new(j,i-1,7);
                    f_new(j,i,8)=f_new(j,i-1,8);
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
    end
    
    %% Collision
    % Computing moments
    for j=1:N_y
        for i=1:N_x
            if Domain_ID(j,i)==1
                Rho(j,i)=sum(f_new(j,i,:));
                u(j,i)=(f_new(j,i,2)+f_new(j,i,6)+f_new(j,i,9)-f_new(j,i,4)-f_new(j,i,7)-f_new(j,i,8))/Rho(j,i);
                v(j,i)=(f_new(j,i,3)+f_new(j,i,6)+f_new(j,i,7)-f_new(j,i,5)-f_new(j,i,8)-f_new(j,i,9))/Rho(j,i);
            end
        end
    end
    % Computing f_eq
    for j=1:N_y
        for i=1:N_x
            if Domain_ID(j,i)==1
                f_eq(j,i,:)=eqm_d2q9(Rho(j,i),[u(j,i);v(j,i)]);
            end
        end
    end
    % Collision & Update
    f=f_new-1/Tau*(f_new-f_eq);

    u_t(:, :, t) = u;
    v_t(:, :, t) = v;
end

%% Visualization
figure;
contourf(flipud(Zone_ID==2),30)
title("Zone_ID==2")
axis equal tight

figure;
contourf(flipud(Zone_ID==1),30)
title("Zone_ID==1")
axis equal tight

figure;
quiver(flipud(u),flipud(v),10)
title("Velocity Plot")
axis equal tight

figure;
contourf(Rho,30)
title("Density Plot")
axis equal tight

dt = datetime('now', 'TimeZone', 'UTC');
customFormat = 'yyyy-MM-dd''T''HH-mm-ss''';
isoDashTimeString = string(dt, customFormat);

disp(['Custom ISO 8601 with dashes: ', isoDashTimeString]);
results_filename = 'results-MLS-' + isoDashTimeString + '.mat'
save(results_filename, "u", "v", "Rho")


%% Velocity Animation
u_anim = VideoWriter('u_animation.avi');
open(u_anim)

figure
ax = gca();
axis([0 N_x 0 N_x])
axis equal

[X,Y] = meshgrid(x, y); % Grid for quiver
disp("Size of coordinate grid, ")
disp(size(X))
disp("Size of velocity grid, ")
disp(size(u_t(:,:,t)))
for t = 1:T
    if t == 1
        h = quiver(ax, X, Y, u_t(:, :, t), v_t(:, :, t), 'LineWidth', 2);
        set(ax, 'XLimMode', 'manual', 'YLimMode', 'manual');
    else
        set(h, 'UData', u_t(t), 'VData', v_t(t));
    end
    
    drawnow();
    writeVideo(u_anim, getframe(ax));
end

close(u_anim);