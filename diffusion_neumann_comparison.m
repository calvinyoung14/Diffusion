% Testing code for solving diffusion equation with Neumann boundary 
% conditions using a synthetic test constructed from exact solution 
% c = sin(x)*exp(-t)

clear;

% domain
xL = -0.5;
xR = 1.5;

% diffusion coefficient
D = 0.6;

% exact solution
c_exact = @(t,x) sin(x)*exp(-t);

% time interval
t_start = 0;
t_final = 1;

% boundary conditions
gL = @(t)  D*cos(xL)*exp(-t);
gR = @(t) -D*cos(xR)*exp(-t);

% source term
f = @(t,x) -sin(x)*exp(-t) + D*sin(x)*exp(-t);

% initial conditions
c_start = @(x) c_exact(t_start, x);

% space discretization
Nx = 20;
x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);

% time-step 
% dt = 5*dx*dx/2/D;
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
A(1,1) = 1 + 2*dt*D/dx/dx;
A(1,2) = -2*dt*D/dx/dx;

A(Nx,Nx) = 1 + 2*dt*D/dx/dx;
A(Nx,Nx-1) = -2*dt*D/dx/dx;
    
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
        
        % boundary points
        A(1,1) = 1 + 2*dt*D/dx/dx;
        A(1,2) = -2*dt*D/dx/dx;
        
        A(Nx,Nx) = 1 + 2*dt*D/dx/dx;
        A(Nx,Nx-1) = -2*dt*D/dx/dx;
    end
    
    % internal points
    for i = 2:Nx-1
        RHS(i) = c_old(i) + dt*f(t+dt,x(i));
    end
    
    % boundary points
    RHS(1) = c_old(1) + dt*f(t+dt,x(i)) - 2*dt*gL(t+dt)/dx;
    RHS(Nx) = c_old(Nx) + dt*f(t+dt,x(i)) - 2*dt*gR(t+dt)/dx;
    
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