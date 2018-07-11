% advection-diffusion equation with no-flux Robin boundary conditions

clear;

% domain
xL = -1;
xR = 3;
yB = -1.5;
yT = 1.5;

% diffusion coefficient
D = 0.7;

% velocity field
vx = -0.8;
vy = -0.4;

% exact solution
c_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% boundary conditions determined by c_exact
c_bc = @(t,x,y) sin(x)*cos(y)*exp(-t);

% source term determined by c_exact
f = @(t,x,y) -sin(x)*cos(y)*exp(-t) + vx*cos(x)*cos(y)*exp(-t) - vy*sin(x)*sin(y)*exp(-t) + 2*D*sin(x)*cos(y)*exp(-t);

% initial conditions determined by c_exact
c_start = @(t,x,y) sin(x)*cos(y)*exp(-t);

%no-flux conditions
g = @(t,x,y) -D*sin(x)*sin(y)*exp(-t) - vy*c_exact(t,x,y);

% time interval
t_start = 0;
t_final = 1.0;

numSplits = 4;
Nx = 20;
Ny = 15;
for k = 1:numSplits
% space discretization

x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);


y = linspace(yB, yT, Ny);
dy = (yT-yB)/(Ny-1);

% time-step 
dt = 0.5*dx;

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

    
while t < t_final
    
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
        
%         %first time
%             % internal points
%         for i = 2:Nx-1
%             for j = 2:Ny-1
%                 p = (j-1)*Nx + i; 
%                 RHS(p) = c_old(p) + dt*f(t+dt,x(i),y(j)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old((j)*Nx + i)-c_old(p))/dy;
%             end
%         end
% 
%         % boundary points
% 
%         % Robin at bottom
%         for i = 2:Nx-1
%             p = (1-1)*Nx + i;
%             RHS(p) = c_old(p) + dt*f(t+dt,x(i),y(j)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old((1)*Nx + i)-c_old(p))/dy - 2*dt*g(t+dt,x(i),y(1))/dy;
%         end
% 
%         % Dirichlet at the left end
%         for j = 1:Ny-1
%             p = (j-1)*Nx + 1;
%             RHS(p) = c_bc(t+dt,x(1),y(j));
%         end
% 
%     %Dirichlet at the right end
%         for j = 1:Ny-1
%             p = (j-1)*Nx + Nx;
%             RHS(p) = c_bc(t+dt,x(Nx),y(j));
%         end
% 
%     %Dirichlet at the top
%         for i = 1:Nx
%             p = (Ny-1)*Nx + i;
%             RHS(p) = c_bc(t+dt,x(i),y(Ny));
%         end

        
    end
    
  %second time
    % internal points
    for i = 2:Nx-1
        for j = 2:Ny-1
            p = (j-1)*Nx + i; 
            RHS(p) = c_old(p) + dt*f(t+dt,x(i),y(j)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old(p+Nx)-c_old(p))/dy;
        end
    end
    
    % boundary points
    
    % Robin at bottom
    for i = 2:Nx-1
        p = (1-1)*Nx + i;
        RHS(p) = c_old(p) + dt*f(t+dt,x(i),y(1)) -dt*vx*(c_old(p+1)-c_old(p))/dx -dt*vy*(c_old(p+Nx)-c_old(p))/dy - 2*dt*g(t+dt,x(i),y(1))/dy;
        if i == 4 && t == 0
%             RHS(p)
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
    t = t+dt;
    
end

    for i = 1:Nx
            for j = 1:Ny
                p = (j-1)*Nx + i;  
                c_exa(p) = c_exact(t,x(i),y(j));
            end
    end
    
    ErrorList = zeros(Nx,Ny);
    tempMax = 0;
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx + i;  
            ErrorList(i,j) = abs(c_exa(p)-c_new(p));
            if(ErrorList(i,j) > tempMax)
                tempMax = ErrorList(i,j);
            end
        end
    end
    
    error(k) = tempMax;


Nx = Nx * 2;
Ny = Ny * 2;
end

for i = 2:numSplits
    order(i) = log(error(i-1)/error(i))/log(2);
end
for i = 1:numSplits
    fprintf('%g \t %g\n', error(i), order(i));
end