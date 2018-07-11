% diffusion equation with no-flux Neumann boundary conditions

clear;

% domain
xL = 0;
xR = 1;

% diffusion coefficient
D = 0.6;

% boundary conditions
gL = @(t) 0.;
gR = @(t) 0.;

% source term
f = @(t,x) 0;

% initial conditions
x0 = 0.5;
tau = 0.001;
c_start = @(x) 3*sqrt(4*pi*D*tau)/sqrt(4*pi*D*(tau)) * exp(-(x-x0)^2/(4*D*(tau)));

% time interval
t_start = 0;
t_final = 1.005;

% space discretization
Nx = 100;
x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);

% time-step 
dt = 5*dx*dx/2/D;
%dt = 0.01*dx;

c_old = zeros(1,Nx);
c_new = zeros(1,Nx);
c_exa = zeros(1,Nx);

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
        
    plot(x,c_new,'-','LineWidth',1);
    xlabel('x');
    ylabel('c');
    axis([xL xR 0 1]);
    pause(dt);
    
end