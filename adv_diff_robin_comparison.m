% Testing code for solving advection-diffusion equation with Robin boundary 
% conditions using a synthetic test constructed from exact solution 
% c = sin(x)*exp(-t)

clear;

% domain
xL = -0.5;
xR = 1.5;

% diffusion coefficient
D = 0.6;

% velocity field
vx = @(t,x) -0.8;

% exact solution
c_exact = @(t,x) sin(x)*exp(-t);

% time interval
t_start = 0;
t_final = 1;

% boundary conditions
gL = @(t) D*cos(xL)*exp(-t) - vx(t,xL) * c_exact(t,xL);
cR = @(t) c_exact(t,xR);

% source term
f = @(t,x) -sin(x)*exp(-t) + D*sin(x)*exp(-t) + vx(t,x)*cos(x)*exp(-t);

% initial conditions
c_start = @(x) c_exact(t_start, x);

% space discretization
Nx = 20;
x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);

% time-step 
% dt = 0.9*dx*dx/2/D;
dt = 0.1*dx;

c_old = zeros(Nx,1);
c_new = zeros(Nx,1);
c_exa = zeros(Nx,1);

for i = 1:Nx
    c_old(i) = c_start(x(i));
end

t = t_start;

% Create sparse matrix and allocate memory for right-hand side
A = sparse(Nx,Nx);
RHS = zeros(Nx,1);

% Calculate the matrix before the while-loop to save time
% internal points
for i = 2:Nx-1
    A(i,i) = 1+2*dt*D/dx/dx;
    A(i,i-1) = -dt*D/dx/dx;
    A(i,i+1) = -dt*D/dx/dx;
end

% boundary points

% Robin at the left end
A(1,1) = 1 + 2*dt*D/dx/dx + 2*dt*vx(t,x(1))/dx;
A(1,2) = -2*dt*D/dx/dx;

% Dirichlet at the right end
A(Nx,Nx) = 1;
    
while t < t_final
    
    if t + dt > t_final
        dt = t_final-t;
            
        % need to recalculate the matrix since dt has changed
        % internal points
        for i = 2:Nx-1
            A(i,i) = 1+2*dt*D/dx/dx;
            A(i,i-1) = -dt*D/dx/dx;
            A(i,i+1) = -dt*D/dx/dx;
        end
        
        A(1,1) = 1 + 2*dt*D/dx/dx + 2*dt*vx(t,x(1))/dx;
        A(1,2) = -2*dt*D/dx/dx;
    end
    
    % internal points
    for i = 2:Nx-1
        velo = vx(t, x(i));
        RHS(i) = c_old(i) + dt*f(t+dt,x(i)) -dt*velo*(c_old(i+1)-c_old(i))/dx;
    end
    
    % boundary points
    
    % Robin at xL
    velo = vx(t, x(1));
    RHS(1) = c_old(1) + dt*f(t+dt,x(1)) -dt*velo*(c_old(2)-c_old(1))/dx - 2*dt*gL(t+dt)/dx;
    
    % Dirichlet at xR
    RHS(Nx) = cR(t+dt);
    
    % solve system of equations
    c_new = A\RHS;
    
    c_old = c_new;
    t = t+dt;
    
    for i = 1:Nx
        c_exa(i) = c_exact(t,x(i));
    end
    
    plot(x,c_exa,'LineWidth',2);
    hold on
    plot(x,c_new,'o','LineWidth',1);
    hold off
    xlabel('x');
    ylabel('c');
    axis([xL xR -1 1]);
    legend('Exact', 'Numerical');
    pause(dt);
    
end

error = max(abs(c_new-c_exa))