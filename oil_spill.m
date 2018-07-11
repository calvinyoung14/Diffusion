% advection-diffusion equation with no-flux Robin boundary conditions

clear;

% domain
xL = 0;
xR = 12;
yB = 0;
yT = 3;

% diffusion coefficient
D = 0.2;

% velocity field
vx = -0.8;
vy = -0.4;

% exact solution
c_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% boundary conditions determined by c_exact
c_bc = @(t,x,y) 0;

% source term determined by c_exact
eps = 0.1;
r_s = 0.1;
x_s = 10;
f1 = @(t,x,y) 0.5*(1-tanh((sqrt((x-x_s)^2 + y^2) - r_s)/eps)); %t<0.5
f2 = @(t,x,y) 0; %t>0.5

% initial conditions determined by c_exact
c_start = @(t,x,y) 0;

%no-flux conditions
g = @(t,x,y) 0;

% time interval
t_start = 0;
t_final = 10.0;
%t_final = 1;
% t_final = 4;
% t_final = 7;

Nx = 160;
Ny = 40;
% space discretization

x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);


y = linspace(yB, yT, Ny);
dy = (yT-yB)/(Ny-1);

% time-step 
dt = 0.1;

c_old = zeros(Nx*Ny,1);
c_new = zeros(Nx*Ny,1);
c_exa = zeros(Nx*Ny,1);

t = t_start;

for i = 1:Nx
    for j = 1:Ny
        p = (j-1)*Nx + i;  
        c_old(p) = c_start(t,x(i),y(j));
    end
end



% Create sparse matrix and allocate memory for right-hand side
A = sparse(Nx*Ny,Nx*Ny);
RHS = zeros(Nx*Ny,1);

% Calculate the matrix before the while-loop to save time
% internal points

aL = -D*dt/dx/dx;
aR = aL;
aC = 1 + 2*D*dt/dx/dx + 2*D*dt/dy/dy;
aT = -D*dt/dy/dy;
aB = aT;

for i = 2:Nx-1
    for j = 2:Ny-1
        p = (j-1)*Nx + i;      
        A(p,p) = aC;                                        %center
        A(p,p-1) = aL;                                      %left
        A(p,p+1) = aR;                                      %right
        A(p,p-Nx) = aB;                                     %bottom
        A(p,p+Nx) = aT;                                     %top
    end
end

% boundary points



% Dirichlet at the left end
for j = 1:Ny-1
    p = (j-1)*Nx + 1;
    A(p,p) = 1;
end

%Dirichlet at the right end
for j = 1:Ny-1
    p = (j-1)*Nx + Nx;
    A(p,p) = 1;
end

%Dirichlet at the top
for i = 1:Nx
    p = (Ny-1)*Nx + i;
    A(p,p) = 1;
end

% Robin at the bottom
for i = 2:Nx-1
    p = (1-1)*Nx + i;
    A(p,p) = aC + 2*(vy*dt/dy);
    A(p,p+1) = aR;
    A(p,p-1) = aL;
    A(p,p+Nx) = aT + aB;
end

index = 1;
while t < t_final
    time(index) = t;
    beach4(index) = c_old((1 - 1)*Nx + ceil(Nx/3));
     beach6(index) = c_old((1 - 1)*Nx + Nx/2);
     beach8(index) = c_old((1 - 1)*Nx + ceil(Nx*(2/3)));
    
    
    if t + dt > t_final
        dt = t_final-t;
            
        
        aL = -D*dt/dx/dx;
        aR = aL;
        aC = 1 + 2*D*dt/dx/dx + 2*D*dt/dy/dy;
        aT = -D*dt/dy/dy;
        aB = aT;
        % need to recalculate the matrix since dt has changed
        % internal points
        for i = 2:Nx-1
            for j = 2:Ny-1
                p = (j-1)*Nx + i;      
                A(p,p) = aC;                                        %center
                A(p,p-1) = aL;                                      %left
                A(p,p+1) = aR;                                      %right
                A(p,p-Nx) = aB;                                     %bottom
                A(p,p+Nx) = aT;                                     %top
            end
        end
        
        % Robin at the bottom
        for i = 2:Nx-1
            p = (1-1)*Nx + i;
            A(p,p) = aC + 2*(vy*dt/dy);
            A(p,p+1) = aR;
            A(p,p-1) = aL;
            A(p,p+Nx) = aT + aB;
        end
        
    end
    
  %second time
    % internal points
    for i = 2:Nx-1
        for j = 2:Ny-1
            p = (j-1)*Nx + i;
            if(t<0.5)
                RHS(p) = c_old(p) + dt*f1(t+dt,x(i),y(j)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old(p+Nx)-c_old(p))/dy;
            else
                RHS(p) = c_old(p) + dt*f2(t+dt,x(i),y(j)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old(p+Nx)-c_old(p))/dy;
            end
        end
    end
    
    % boundary points
    
    % Robin at bottom
    for i = 2:Nx-1
        p = (1-1)*Nx + i;
        if(t<0.5)
            RHS(p) = c_old(p) + dt*f1(t+dt,x(i),y(1)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old(p+Nx)-c_old(p))/dy - 2*dt*g(t+dt,x(i),y(1))/dy;
        else
            RHS(p) = c_old(p) + dt*f2(t+dt,x(i),y(1)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old(p+Nx)-c_old(p))/dy - 2*dt*g(t+dt,x(i),y(1))/dy;
        end
    end
    
    % Dirichlet at the left end
    for j = 1:Ny
        p = (j-1)*Nx + 1;
        RHS(p) = c_bc(t+dt,x(1),y(j));
    end

%Dirichlet at the right end
    for j = 1:Ny
        p = (j-1)*Nx + Nx;
        RHS(p) = c_bc(t+dt,x(Nx),y(j));
    end

%Dirichlet at the top
    for i = 1:Nx
        p = (Ny-1)*Nx + i;
        RHS(p) = c_bc(t+dt,x(i),y(Ny));
    end

    % solve system of equations
    c_new = A\RHS;
    
    c_old = c_new;
    
    contourPlot = zeros(Nx,Ny);
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx + i;  
            contourPlot(i,j) = c_new(p);
        end
    end
    
    contourf(x,y,contourPlot','LevelList',linspace(0,.05,100),'LineColor', 'none');
    colorbar;
    hold on;
    hold off;
    axis([xL, xR, yB, yT]);
    axis equal;
    pause(0.01);
    
    if(t>.99 && t<1.01)
        print('oil1','-dpng');
    end
    if(t>3.99 && t<4.01)
        print('oil4','-dpng');
    end
    if(t>6.99 && t<7.01)
        print('oil7','-dpng');
    end
    
    
    t = t+dt;
    index= index + 1;
end
figure;
graph = plot(time, beach4, time, beach6, time, beach8);
set(graph,'lineWidth', 2.0);
legend(graph,'beach4','beach6','beach8');