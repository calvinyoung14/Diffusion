clear;

% domain
xL = -2;
xR = 2;
yB = -2.5;
yT = 2.5;

% thermal diffusivity
lambda = .0015;

% time interval
t_start = 0;
t_final = 1500;

% boundary conditions ____calculated from t_exact
t_bc = @(t) min(20+80*(t/60),100);

% source term _____ calculated from t_exact
f = @(t,x,y) 0;

% initial conditions  _____calculated from t_exact
c_start = @(t,x,y) 20;


% exact solution
t_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

Nx = 80;
Ny = 100;
% space discretization

%Nx = 50;
%Nx = 100;
x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);


%Ny = 60;
%Ny = 120;
y = linspace(yB, yT, Ny);
dy = (yT-yB)/(Ny-1);


% time-step 
dt = 5;


c_old = zeros(1,Nx*Ny);
c_new = zeros(1,Nx*Ny);
c_exa = zeros(1,Nx*Ny);

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
for i = 2:Nx-1
    for j = 2:Ny-1
        p = (j-1)*Nx + i;
        A(p,p) = 1 + 2*lambda*dt/dx/dx + 2*lambda*dt/dy/dy;          %center
        A(p,p-1) = -lambda*dt/dx/dx;                                      %left
        A(p,p+1) = -lambda*dt/dx/dx;                                      %right
        A(p,p-Nx) = -lambda*dt/dy/dy;                                     %bottom
        A(p,p+Nx) = -lambda*dt/dy/dy;                                     %top
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
        p = (1-1)*Nx + i;
        A(p,p) = 1;
        p1 = (Ny-1)*Nx + i;
        A(p1,p1) = 1;
    end
end
  
index = 1;

while t < t_final
    time(index) = t;
    tempCenter(index) = c_old((Ny/2-1)*Nx + Nx/2);
    if t + dt > t_final
        dt = t_final-t;
            
        % need to recalculate the matrix since dt has changed
        % internal points
        for i = 2:Nx-1
            for j = 2:Ny-1
                p = (j-1)*Nx + i;      
                A(p,p) = 1 + 2*lambda*dt/dx/dx + 2*lambda*dt/dy/dy;          %center
                A(p,p-1) = -lambda*dt/dx/dx;                                      %left
                A(p,p+1) = -lambda*dt/dx/dx;                                      %right
                A(p,p-Nx) = -lambda*dt/dy/dy;                                     %bottom
                A(p,p+Nx) = -lambda*dt/dy/dy;                                     %top
            end
        end
        
        % internal points
    for i = 2:Nx-1
        for j = 2:Ny-1
            p = (j-1)*Nx + i;  
            RHS(p) = c_old(p) + dt*f(t+dt,x(i),y(j));
        end
    end
    
    %boundary points
    for i = 1:Nx
        if i == 1
            for j = 1:Ny
                p = (j-1)*Nx + i;
                RHS(p) = t_bc(t+dt);
            end
        elseif i == Nx
            for j = 1:Ny
                p = (j-1)*Nx + i;
                RHS(p) = t_bc(t+dt);
            end
        else
            p = (1-1)*Nx + i;
            RHS(p) = t_bc(t+dt);
            p1 = (Ny-1)*Nx + i;
            RHS(p1) = t_bc(t+dt);
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
                RHS(p) = t_bc(t+dt);
            end
        elseif i == Nx
            for j = 1:Ny
                p = (j-1)*Nx + i;
                RHS(p) = t_bc(t+dt);
            end
        else
            p = (1-1)*Nx + i;
            RHS(p) = t_bc(t+dt);
            p1 = (Ny-1)*Nx + i;
            RHS(p1) = t_bc(t+dt);
        end
    end
    
    % solve system of equations
    c_new = A\RHS;
    
    c_old = c_new;
    
    index = index + 1;
    
    contourPlot = zeros(Nx,Ny);
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx + i;  
            contourPlot(i,j) = c_new(p);
        end
    end
    
    contourf(x,y,contourPlot','LevelList',linspace(20.0,100.0,100),'LineColor', 'none');
    colorbar;
    hold on;
    hold off;
    axis([xL, xR, yB, yT]);
    axis equal;
    pause(0.01);
    
    if(t==0)
        print('snapshot0','-dpng');
    end
    if(t==200)
        print('snapshot200','-dpng');
    end
    if(t==400)
        print('snapshot400','-dpng');
    end
    if(t==600)
        print('snapshot600','-dpng');
    end
    
    t = t+dt;
end

 plot(time, tempCenter);

%     for i = 1:Nx
%         for j = 1:Ny
%             p = (j-1)*Nx + i;  
%             c_exa(p) = t_exact(t,x(i),y(j));
%         end
%     end

%     ErrorList = zeros(Nx,Ny);
%     tempMax = 0;
%     for i = 1:Nx
%         for j = 1:Ny
%             p = (j-1)*Nx + i;  
%             ErrorList(i,j) = abs(c_exa(p)-c_new(p));
%             if(ErrorList(i,j) > tempMax)
%                 tempMax = ErrorList(i,j);
%             end
%         end
%     end
%     
%     error(k) = tempMax;
% 
% Nx = Nx *2;
% Ny = Ny *2;

% for i = 2:numSplits
%     order(i) = log(error(i-1)/error(i))/log(2);
% end
% for i = 1:numSplits
%     fprintf('%g \t %g\n', error(i), order(i));
% end
