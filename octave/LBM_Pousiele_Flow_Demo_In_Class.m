%% Define parameters
% Domain
N_x=100;
N_y=50;
dx=1;
dy=1;
% LBM related
Ksi=[0 1 0 -1  0 1 -1  -1  1;...
     0 0 1  0 -1 1  1  -1 -1];
w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
c_s=1/sqrt(3);
Tau=1;
%% Initialization
Rho_in=1;
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
disp(f_new)


T=10000;
%% Solving
for t=1:T
    % Streaming
    for j=1:N_y
        for i=1:N_x
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

                    f_new(j,i,5)=f_new(j,i,3);
                    f_new(j,i,8)=(f_new(j,i,2)-f_new(j,i,4))/2+f_new(j,i,6);
                    f_new(j,i,9)=(f_new(j,i,4)-f_new(j,i,2))/2+f_new(j,i,7);
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

                U_in=1-(f_new(j,i,1)+f_new(j,i,3)+f_new(j,i,5)+2*(f_new(j,i,4)+f_new(j,i,7)+f_new(j,i,8)))/Rho_in;
                f_new(j,i,2)=f_new(j,i,4)+U_in*Rho_in*2/3;
                f_new(j,i,6)=f_new(j,i,8)+(f_new(j,i,5)-f_new(j,i,3))/2+U_in*Rho_in/6;
                f_new(j,i,9)=f_new(j,i,7)-(f_new(j,i,5)-f_new(j,i,3))/2+U_in*Rho_in/6;
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
    % Collsion
    % Computing moments
    for j=1:N_y
        for i=1:N_x
            Rho(j,i)=sum(f_new(j,i,:));
            u(j,i)=(f_new(j,i,2)+f_new(j,i,6)+f_new(j,i,9)-f_new(j,i,4)-f_new(j,i,7)-f_new(j,i,8))/Rho(j,i)
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

%% Visualization
% Analytical Solution
y_benchmark=0:0.01:1;
for i=1:length(y_benchmark)
    u_benchmark(i)=-4*(y_benchmark(i)^2-y_benchmark(i));
end
figure;
plot(y_benchmark,u_benchmark,"red")
% Simultion result
for j=1:N_y
    u_sim(j)=u(j,N_x-1);
end
hold on
plot((0:1:N_y-1)/(N_y-1),u_sim/max(u_sim),'blue');

figure
quiver(flipud(u),flipud(v),10)
axis equal tight

figure
contourf(Rho,30)
axis equal tight
