clear;

% domain
xL = -2;
xR = 2;
yB = -2.5;
yT = 2.5;

% thermal diffusivity
lambda = 1.5*(10^-3);

% time interval
t_start = 0;
t_final = 1500;

% boundary conditions ____calculated from t_exact
t_bc = @(t,x,y) min(20+80*(t/60),100);

% source term
f = @(t,x,y) 0;

% initial conditions  in degrees
c_start = @(x,y) 20;


% exact solution
t_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% space discretization
Nx = 80;
x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);

Ny = 100;
y = linspace(yB, yT, Ny);
dy = (yT-yB)/(Ny-1);

% time-step 
dt = 5;


c_old = zeros(1,Nx*Ny);
c_new = zeros(1,Nx*Ny);
c_exa = zeros(1,Nx*Ny);

for i = 1:Nx
    for j = 1:Ny
        p = (j-1)*Nx + i;  
        c_old(p) = c_start(i,j);
    end
    
end

t = t_start;

% Create sparse matrix and allocate memory for right-hand side
A = sparse(Nx*Ny,Nx*Ny);
RHS = zeros(Nx*Ny,1);

% Calculate the matrix before the while-loop to save time
% internal points
for i = 2:Nx-1
    for j = 2:Ny-1
        p = (j-1)*Nx + i;      
        A(p,p) = (1 + 2*(lambda*dt/(dx)^2) + 2*(lambda*dt/(dy)^2));          %center
        A(p,p-1) = -(lambda*dt/(dx)^2);                                      %left
        A(p,p+1) = -(lambda*dt/(dx)^2);                                      %right
        A(p,p-Nx) = -(lambda*dt/(dy)^2);                                     %bottom
        A(p,p+Nx) = -(lambda*dt/(dy)^2);                                     %top
    end
end

% boundary points
for i = 1:Nx
    if i == 1
        for j = 1:Ny
            p = (j-1)*Nx + i;
            A(p,p) = 1;
        end
    elseif i == Nx
        for j = 1:Ny
            p = (j-1)*Nx + i;
            A(p,p) = 1;
        end
    else
        A(i,i) = 1;
        p1 = (Ny-1)*Nx + i;
        A(p1,p1) = 1;
    end
end

    
while t < t_final
    
    if t + dt > t_final
        dt = t_final-t;
            
        % need to recalculate the matrix since dt has changed
        % internal points
        for i = 2:Nx-1
            for j = 2:Ny-1
                p = (j-1)*Nx + i;      
                A(p,p) = (1 + 2*(lambda*dt/(dx)^2) + 2*(lambda*dt/(dy)^2));          %center
                A(p,p-1) = -(lambda*dt/(dx)^2);                                      %left
                A(p,p+1) = -(lambda*dt/(dx)^2);                                      %right
                A(p,p-Nx) = -(lambda*dt/(dy)^2);                                     %bottom
                A(p,p+Nx) = -(lambda*dt/(dy)^2);                                     %top
            end
        end
    end
    
    % internal points
    for i = 2:Nx-1
        for j = 2:Ny-1
            p = (j-1)*Nx + i;  
            RHS(p) = c_old(p) + dt*f(t+dt,x(i),y(j));
        end
    end
    
    % boundary points
    for i = 1:Nx
        if i == 1
            for j = 1:Ny
                p = (j-1)*Nx + i;
                RHS(p) = t_bc(t+dt,x(i),y(j));
            end
        elseif i == Nx
            for j = 1:Ny
                p = (j-1)*Nx + i;
                RHS(p) = t_bc(t+dt,x(i),y(j));
            end
        else
            RHS(i) = t_bc(t+dt,x(i),y(1));
            p1 = (Ny-1)*Nx + i;
            RHS(p1) = t_bc(t+dt,x(i),y(Ny));
        end
    end
    
    % solve system of equations
    c_new = A\RHS;
    
    c_old = c_new;
    t = t+dt;
    contourPlot = zeros(Nx,Ny);
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx + i;  
            c_exa(p) = t_exact(t,x(i),y(j));
            contourPlot(i,j) = c_new(p);
        end
    end
    
    
    contourf(x,y,contourPlot','LevelList',linspace(-75.0,75.0,200),'LineColor', 'none');
    colorbar;
    hold on;
    hold off;
    axis([xL, xR, yB, yT]);
    axis equal;
    pause(0.01);
end