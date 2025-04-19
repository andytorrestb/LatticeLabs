clc;
clear all;

%% Define physical and dimensional information
Rho_ref=2; % Reference density
U_lid=0.1; % This is the x-velocity on the moving lid

L=100; % The length of the entire computational domain
H=L; % The height of the entire computational domain
N_nd_x=50; % The total number of nodes in x-direction
N_nd_y=N_nd_x; % The total number of nodes in y-direction
dx=L/(N_nd_x-1); % The distance between the centroids of two horizontally neighboring CVs
dy=H/(N_nd_y-1); % The distance between the centroids of two vertically neighboring CVs

N_cv_x=N_nd_x-1; % The total number of CVs in x-direction
N_cv_y=N_nd_y-1; % The total number of CVs in y-direction

A=dx*dy; % The area (which is the volume in 2D) of each CV
n1=[1;0]; % n1 is the unit outward normal vector on face 1 (facing east)
n2=[0;-1]; % n2 is the unit outward normal vector on face 2 (facing south)
n3=[-1;0]; % n3 is the unit outward normal vector on face 3 (facing west)
n4=[0;1]; % n4 is the unit outward normal vector on face 4 (facing north)

%% Define Boltzmann related information
c_s=1/sqrt(3);
dt=0.1;
Tau=0.3;
Ksi=[0 1 0 -1  0 1 -1 -1  1;...
     0 0 1  0 -1 1  1 -1 -1 ];
w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];

%% Initialization
F=zeros(N_cv_y,N_cv_x,9); % The flux term in all direction at all CV centroids
f_eq=zeros(N_cv_y,N_cv_x,9); % Equilibrium PDF in all direction at all CV centroids
for i=1:N_cv_x
    for j=1:N_cv_y
        for k=1:9
            f_eq(j,i,k)=w(k)*Rho_ref*(1+Ksi(:,k)'*[0;0]/c_s^2+(Ksi(:,k)'*[0;0])^2/c_s^4/2-[0,0]*[0;0]/c_s^2/2);
        end
    end
end
f=f_eq; % PDF in all directions at all CV centroids
f_new=f; % New PDF in all directions at all CV centroids
f_neq=f-f_eq; % Non-equilibrium PDF in all direction at all CV centroids
u=zeros(N_cv_y,N_cv_x); % x-component of physical velocity at all CV centroids
v=zeros(N_cv_y,N_cv_x); % y-component of physical velocity at all CV centroids
Rho=ones(N_cv_y,N_cv_x)*Rho_ref; % The densities at all CV centroids

f_nd_eq=zeros(N_nd_y,N_nd_x,9); % Equilibrium PDF in all direction at all nodes
f_nd=f_nd_eq; % PDF in all direction at all nodes
f_nd_neq=f_nd-f_nd_eq; % Non-equilibrium PDF in all direction at all nodes
u_nd=zeros(N_nd_y,N_nd_x); % x-component of physical velocity at all nodes
v_nd=zeros(N_nd_y,N_nd_x); % y-component of physical velocity at all nodes
Rho_nd=ones(N_nd_y,N_nd_x)*Rho_ref; % The densities at all nodes

T=10000;
%% Solver
for t=1:T
    %% Boundary conditions-compute f_nd for all nodes on the boundaries
    for j=1:N_nd_y
        for i=1:N_nd_x
            if j==1  % This are the nodes on the top boundary
                if i==1 % This is the top-left corner node
                    %% f_neq
                    f_cv=f(j,i,:);
                    f_eq_cv=f_eq(j,i,:);
                    f_neq_cv=f_cv-f_eq_cv;
                    f_nd_neq(j,i,:)=f_neq_cv;
                    %% f_eq
                    Rho_nd(j,i)=Rho(j,i);
                    for k=1:9
                        f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i)*(1+Ksi(:,k)'*[U_lid;0]/c_s^2+...
                            (Ksi(:,k)'*[U_lid;0])^2/2/c_s^4-...
                            (U_lid^2)/2/c_s^2);
                    end
                    %% Total f
                    f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
                elseif i==N_nd_x % This is the top-right corner node
                    %% f_neq
                    f_cv=f(j,i-1,:);
                    f_eq_cv=f_eq(j,i-1,:);
                    f_neq_cv=f_cv-f_eq_cv;
                    f_nd_neq(j,i,:)=f_neq_cv;
                    %% f_eq
                    Rho_nd(j,i)=Rho(j,i-1);
                    for k=1:9
                        f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i)*(1+Ksi(:,k)'*[U_lid;0]/c_s^2+...
                            (Ksi(:,k)'*[U_lid;0])^2/2/c_s^4-...
                            (U_lid^2)/2/c_s^2);
                    end
                    %% Total f
                    f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
                else % This are all other nodes the top boundary
                    %% f_neq
                    f_cv_1=f(j,i-1,:);
                    f_cv_2=f(j,i,:);

                    f_eq_cv_1=f_eq(j,i-1,:);
                    f_eq_cv_2=f_eq(j,i,:);

                    f_neq_cv_1=f_cv_1-f_eq_cv_1;
                    f_neq_cv_2=f_cv_2-f_eq_cv_2;
                    f_nd_neq(j,i,:)=(f_neq_cv_1+f_neq_cv_2)/2;
                    %% f_eq
                    Rho_nd(j,i)=(Rho(j,i-1)+Rho(j,i))/2;
                    for k=1:9
                        f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i)*(1+Ksi(:,k)'*[U_lid;0]/c_s^2+...
                            (Ksi(:,k)'*[U_lid;0])^2/2/c_s^4-...
                            (U_lid^2)/2/c_s^2);
                    end
                    %% Total f
                    f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
                end
            elseif j==N_nd_y % This are the nodes on the bottom boundary
                if i==1 % This is the bottom-left corner node
                    %% f_neq
                    f_cv=f(end,i,:);
                    f_eq_cv=f_eq(end,i,:);
                    f_neq_cv=f_cv-f_eq_cv;
                    f_nd_neq(j,i,:)=f_neq_cv;
                    %% f_eq
                    Rho_nd(j,i)=Rho(end,i);
                    for k=1:9
                        f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i);
                    end
                    %% Total f
                    f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
                elseif i==N_nd_x % This is the bottom-right corner node
                    f_cv=f(end,i-1,:);
                    f_eq_cv=f_eq(end,i-1,:);
                    f_neq_cv=f_cv-f_eq_cv;
                    f_nd_neq(j,i,:)=f_neq_cv;
                    %% f_eq
                    Rho_nd(j,i)=Rho(end,i-1);
                    for k=1:9
                        f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i);
                    end
                    %% Total f
                    f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
                else % This are all other nodes on the bottom boundary
                    %% f_neq
                    f_cv_1=f(N_cv_y,i-1,:);
                    f_cv_2=f(end,i,:);

                    f_eq_cv_1=f_eq(end,i-1,:);
                    f_eq_cv_2=f_eq(end,i,:);

                    f_neq_cv_1=f_cv_1-f_eq_cv_1;
                    f_neq_cv_2=f_cv_2-f_eq_cv_2;
                    f_nd_neq(j,i,:)=(f_neq_cv_1+f_neq_cv_2)/2;
                    %% f_eq
                    Rho_nd(j,i)=(Rho(end,i-1)+Rho(end,i))/2;
                    for k=1:9
                        f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i);
                    end
                    %% Total f
                    f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
                end
            elseif i==1 % This are the nodes on the left boundary
                %% f_neq
                f_cv_1=f(j-1,i,:);
                f_cv_2=f(j,i,:);

                f_eq_cv_1=f_eq(j-1,i,:);
                f_eq_cv_2=f_eq(j,i,:);

                f_neq_cv_1=f_cv_1-f_eq_cv_1;
                f_neq_cv_2=f_cv_2-f_eq_cv_2;
                f_nd_neq(j,i,:)=(f_neq_cv_1+f_neq_cv_2)/2;
                %% f_eq
                Rho_nd(j,i)=(Rho(j-1,i)+Rho(j,i))/2;
                for k=1:9
                    f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i);
                end
                %% Total f
                f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
            elseif i==N_nd_x % This are the nodes on the right boundary
                %% f_neq
                f_cv_1=f(j-1,end,:);
                f_cv_2=f(j,end,:);

                f_eq_cv_1=f_eq(j-1,end,:);
                f_eq_cv_2=f_eq(j,end,:);

                f_neq_cv_1=f_cv_1-f_eq_cv_1;
                f_neq_cv_2=f_cv_2-f_eq_cv_2;
                f_nd_neq(j,i,:)=(f_neq_cv_1+f_neq_cv_2)/2;
                %% f_eq
                Rho_nd(j,i)=(Rho(j-1,end)+Rho(j,end))/2;
                for k=1:9
                    f_nd_eq(j,i,k)=w(k)*Rho_nd(j,i);
                end
                %% Total f
                f_nd(j,i,:)=f_nd_eq(j,i,:)+f_nd_neq(j,i,:);
            else     % All interior nodes -- DO NOTHING
                ;
            end
        end
    end
    %% F - the flux term in the governing equation
    for j=1:N_cv_y
        for i=1:N_cv_x
            if j==1  % This are the CVs touching the top boundary
                if i==1 % This is the CV touching the top-left corner node
                    %% Flux on face 1
                    F1=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i+1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n1<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc1=fu;
                    % Compute flux on the face
                    for k=1:9
                        F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                    end
                    %% Flux on face 2
                    F2=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j+1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n2<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc2=fu;
                    % Compute flux on the face
                    for k=1:9
                        F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                    end
                    %% Flux on face 3
                    F3=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(j,1,:))+squeeze(f_nd(j+1,1,:)))/2;% This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n3<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc3=fu;
                    % Compute flux on the face
                    for k=1:9
                        F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                    end
                    %% Flux on face 4 - This is the face on the top boundary
                    F4=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(1,i,:))+squeeze(f_nd(1,i+1,:)))/2;% This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n4<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc4=fu;
                    % Compute flux on the face
                    for k=1:9
                        F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                    end
                    %% Final flux term in the governing equation
                    F(j,i,:)=(F1+F2+F3+F4)/A;
                elseif i==N_cv_x % This is the CV touching the top-right corner node
                    %% Flux on face 1
                    F1=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(j,end,:))+squeeze(f_nd(j+1,end,:)))/2;% This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n1<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc1=fu;
                    % Compute flux on the face
                    for k=1:9
                        F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                    end
                    %% Flux on face 2
                    F2=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j+1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n2<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc2=fu;
                    % Compute flux on the face
                    for k=1:9
                        F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                    end
                    %% Flux on face 3
                    F3=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i-1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n3<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc3=fu;
                    % Compute flux on the face
                    for k=1:9
                        F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                    end
                    %% Flux on face 4 - This is the face on the top boundary
                    F4=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(1,i,:))+squeeze(f_nd(1,i+1,:)))/2;% This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n4<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc4=fu;
                    % Compute flux on the face
                    for k=1:9
                        F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                    end
                    %% Final flux term in the governing equation
                    F(j,i,:)=(F1+F2+F3+F4)/A;
                else % This are all other CVs touching the top boundary
                    %% Flux on face 1
                    F1=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i+1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n1<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc1=fu;
                    % Compute flux on the face
                    for k=1:9
                        F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                    end
                    %% Flux on face 2
                    F2=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j+1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n2<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc2=fu;
                    % Compute flux on the face
                    for k=1:9
                        F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                    end
                    %% Flux on face 3
                    F3=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i-1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n3<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc3=fu;
                    % Compute flux on the face
                    for k=1:9
                        F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                    end
                    %% Flux on face 4 - This is the face on the top boundary
                    F4=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(1,i,:))+squeeze(f_nd(1,i+1,:)))/2;% This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n4<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc4=fu;
                    % Compute flux on the face
                    for k=1:9
                        F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                    end
                    %% Final flux term in the governing equation
                    F(j,i,:)=(F1+F2+F3+F4)/A;
                end
            elseif j==N_cv_y % This are CVs touching the bottom boundary
                if i==1 % This is CV touching the bottom-left corner node
                    %% Flux on face 1
                    F1=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i+1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n1<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc1=fu;
                    % Compute flux on the face
                    for k=1:9
                        F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                    end
                    %% Flux on face 2 - this is the face on the bottom boundary
                    F2=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(end,i,:))+squeeze(f_nd(end,i+1,:)))/2;% This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n2<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc2=fu;
                    % Compute flux on the face
                    for k=1:9
                        F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                    end
                    %% Flux on face 3 - this is the face on the left boundary
                    F3=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(j,1,:))+squeeze(f_nd(j+1,1,:)))/2; % This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n3<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc3=fu;
                    % Compute flux on the face
                    for k=1:9
                        F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                    end
                    %% Flux on face 4
                    F4=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j-1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n4<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc4=fu;
                    % Compute flux on the face
                    for k=1:9
                        F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                    end
                    %% Final flux term in the governing equation
                    F(j,i,:)=(F1+F2+F3+F4)/A;
                elseif i==N_cv_x % This is CV touching the bottom-right corner node
                    %% Flux on face 1 - This is the face touching the right boundary
                    F1=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(j,end,:))+squeeze(f_nd(j+1,end,:)))/2; % This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n1<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc1=fu;
                    % Compute flux on the face
                    for k=1:9
                        F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                    end
                    %% Flux on face 2 - this is the face on the bottom boundary
                    F2=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(end,i,:))+squeeze(f_nd(end,i+1,:)))/2; % This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n2<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc2=fu;
                    % Compute flux on the face
                    for k=1:9
                        F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                    end
                    %% Flux on face 3
                    F3=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i-1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n3<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc3=fu;
                    % Compute flux on the face
                    for k=1:9
                        F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                    end
                    %% Flux on face 4
                    F4=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j-1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n4<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc4=fu;
                    % Compute flux on the face
                    for k=1:9
                        F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                    end
                    %% Final flux term in the governing equation
                    F(j,i,:)=(F1+F2+F3+F4)/A;
                else % This are all other CVs touching the bottom boundary
                    %% Flux on face 1
                    F1=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i+1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n1<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc1=fu;
                    % Compute flux on the face
                    for k=1:9
                        F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                    end
                    %% Flux on face 2 - this is the face on the bottom boundary
                    F2=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    f_nd_center=(squeeze(f_nd(end,i,:))+squeeze(f_nd(end,i+1,:)))/2; % This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                    fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n2<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc2=fu;
                    % Compute flux on the face
                    for k=1:9
                        F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                    end
                    %% Flux on face 3
                    F3=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j,i-1,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n3<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc3=fu;
                    % Compute flux on the face
                    for k=1:9
                        F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                    end
                    %% Flux on face 4
                    F4=zeros(9,1);
                    fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                    fd=squeeze(f(j-1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                    % Determine which is upwind or downwind from the
                    % perspective of face lattice velocity
                    for k=1:9
                        if Ksi(:,k)'*n4<=0
                            temp=fu(k,1);
                            fu(k,1)=fd(k,1);
                            fd(k,1)=temp;
                        end
                    end
                    % Flux scheme, this is where you change it to 2nd-order
                    % Lax-Wendroff scheme
                    fc4=fu;
                    % Compute flux on the face
                    for k=1:9
                        F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                    end
                    %% Final flux term in the governing equation
                    F(j,i,:)=(F1+F2+F3+F4)/A;
                end
            elseif i==1 % This are the CVs touching the left boundary
                %% Flux on face 1
                F1=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j,i+1,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n1<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc1=fu;
                % Compute flux on the face
                for k=1:9
                    F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                end
                %% Flux on face 2
                F2=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j+1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n2<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc2=fu;
                % Compute flux on the face
                for k=1:9
                    F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                end
                %% Flux on face 3 - this is the face on the left boundary
                F3=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                f_nd_center=(squeeze(f_nd(j,1,:))+squeeze(f_nd(j+1,1,:)))/2; % This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n3<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc3=fu;
                % Compute flux on the face
                for k=1:9
                    F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                end
                %% Flux on face 4
                F4=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j-1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n4<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc4=fu;
                % Compute flux on the face
                for k=1:9
                    F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                end
                %% Final flux term in the governing equation
                F(j,i,:)=(F1+F2+F3+F4)/A;
            elseif i==N_cv_x % This are the CVs touching the right boundary
                %% Flux on face 1 - this is the face on the right boundary
                F1=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                f_nd_center=(squeeze(f_nd(j,end,:))+squeeze(f_nd(j+1,end,:)))/2; % This is PDF calculated by averaging the PDFs at the two nodes bounding the current face
                fd=2*f_nd_center-fu; % This is PDF at the ghost cell center, which is also the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n1<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc1=fu;
                % Compute flux on the face
                for k=1:9
                    F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                end
                %% Flux on face 2
                F2=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j+1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n2<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc2=fu;
                % Compute flux on the face
                for k=1:9
                    F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                end
                %% Flux on face 3
                F3=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j,i-1,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n3<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc3=fu;
                % Compute flux on the face
                for k=1:9
                    F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                end
                %% Flux on face 4
                F4=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j-1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n4<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc4=fu;
                % Compute flux on the face
                for k=1:9
                    F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                end
                %% Final flux term in the governing equation
                F(j,i,:)=(F1+F2+F3+F4)/A;
            else     % All interior CVs
                %% Flux on face 1
                F1=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j,i+1,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n1<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc1=fu;
                % Compute flux on the face
                for k=1:9
                    F1(k,1)=fc1(k,1)*(Ksi(:,k)'*n1)*dy;
                end
                %% Flux on face 2
                F2=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j+1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n2<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc2=fu;
                % Compute flux on the face
                for k=1:9
                    F2(k,1)=fc2(k,1)*(Ksi(:,k)'*n2)*dx;
                end
                %% Flux on face 3
                F3=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j,i-1,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n3<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc3=fu;
                % Compute flux on the face
                for k=1:9
                    F3(k,1)=fc3(k,1)*(Ksi(:,k)'*n3)*dy;
                end
                %% Flux on face 4
                F4=zeros(9,1);
                fu=squeeze(f(j,i,:)); % This is the upwind CV from the perspective of face outward unit vector
                fd=squeeze(f(j-1,i,:)); % This is the downwind CV from the perspective of face outward unit vector
                % Determine which is upwind or downwind from the
                % perspective of face lattice velocity
                for k=1:9
                    if Ksi(:,k)'*n4<=0
                        temp=fu(k,1);
                        fu(k,1)=fd(k,1);
                        fd(k,1)=temp;
                    end
                end
                % Flux scheme, this is where you change it to 2nd-order
                % Lax-Wendroff scheme
                fc4=fu;
                % Compute flux on the face
                for k=1:9
                    F4(k,1)=fc4(k,1)*(Ksi(:,k)'*n4)*dx;
                end
                %% Final flux term in the governing equation
                F(j,i,:)=(F1+F2+F3+F4)/A;
            end
        end
    end
    %% Moments
    for j=1:N_cv_y
        for i=1:N_cv_x
            Rho(j,i)=sum(f(j,i,1:9));
            u(j,i)=(f(j,i,2)+f(j,i,6)+f(j,i,9)-...
                    f(j,i,4)-f(j,i,7)-f(j,i,8))/Rho(j,i);
            v(j,i)=(f(j,i,3)+f(j,i,6)+f(j,i,7)-...
                    f(j,i,5)-f(j,i,8)-f(j,i,9))/Rho(j,i);
        end
    end
    %% f_eq
    for j=1:N_cv_y
        for i=1:N_cv_x
            for k=1:9
                f_eq(j,i,k)=w(k)*Rho(j,i)*(1+Ksi(:,k)'*[u(j,i);v(j,i)]/c_s^2+...
                                             (Ksi(:,k)'*[u(j,i);v(j,i)])^2/2/c_s^4-...
                                             (u(j,i)^2+v(j,i)^2)/2/c_s^2);
            end
        end
    end
    %% Time marching
    f_new=((Tau-dt)/Tau)*f+(dt/Tau)*f_eq-dt*F;

    %% Update solution
    f=f_new;
end

%% Solution Visual Check
% y vs. u
u_sim=zeros(N_cv_y+2,1);
u_sim(1,1)=U_lid;
u_sim(end,1)=0;
for j=1:N_cv_y
    u_sim(1+j,1)=u(j,(N_cv_x-1)/2+1);
end
figure
plot(flipud(u_sim)/U_lid,[0,(dy/2:dy:H-dy/2),H]/H);

% v vs. x
v_sim=zeros(N_cv_x+2,1);
v_sim(1,1)=0;
v_sim(end,1)=0;
for i=1:N_cv_x
    v_sim(1+i,1)=v((N_cv_y-1)/2+1,i);
end
figure
plot([0,(dx/2:dx:L-dx/2),H]/H,v_sim/U_lid);

% Velocity quiver
figure
quiver(flipud(u),flipud(v),10)
axis equal tight

% Density contour
figure
contourf(flipud(Rho),30)
axis equal tight
