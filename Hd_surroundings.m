%{
    This code solves the poisson equation for the demag field
    using Dirichlet and Neuman border conditions, meshing both 
    the refrigerant and the air around it.

    Dif. eq.:
    Nabla \dot Nabla V = Nabla dot M 

    Border conditions for the refrigerant (Neumann)
    nabla V dot n = M dot n

    Border conditions for the container (Dirichlet)
    V = 0
%}

clear
close all
clc

% Dimensions and spacing
Lx = 1;
Ly = 1;
dx = 0.01;

padding = 100;

Lx_container = Lx + (2*padding)*dx; % Add a padding
Ly_container = Ly + (2*padding)*dx; % Add a padding

% Calculate number of grid points
nx = Lx_container / dx + 1;
ny = Ly_container / dx + 1;

% Initialize u
u = zeros(ny, nx);

% Define the Magnetization
% Mx = zeros(ny-2*padding,nx-2*padding); % Only of the magnetic material
Mx = zeros(ny,nx);
My = Mx;

% Mx(:) = 1; % Set all vectors point in XX
for i = 1:nx
    for j = 1:ny
        in_x = (i>padding) && (i<=nx-padding);
        in_y = (j>padding) && (j<=ny-padding);
        if in_x && in_y
            % ip = i-padding;
            % jp = j-padding;
            % Mx(jp,ip) = 1;
            Mx(j,i) = 1;
        end
    end
end

% Iterative solver parameters
max_iter = 1e5;
tol = 1e-6;

% Perform iterations to solve the interior points
for iter = 1:max_iter
    u_old = u;

    % Update interior points
    for i = 2:nx-1
        for j = 2:ny-1
            in_x = (i>padding) && (i<=nx-padding);
            in_y = (j>padding) && (j<=ny-padding);

            if in_x && in_y
                F = (Mx(j,i+1) - Mx(j,i) + ...
                     My(j+1,i) - My(j,i))/dx;
                % ip = i - padding; jp = j - padding;
                % F = (Mx(jp,ip+1) - Mx(jp,ip) + ...
                     % My(jp+1,ip) - My(jp,ip))/dx;
            else
                F = 0;
            end

            u(j, i) = 0.25 * (u_old(j+1, i) + u_old(j-1, i) + ...
                              u_old(j, i+1) + u_old(j, i-1) - ...
                              dx^2 * F);
        end
    end

    % Apply Dirichlet boundary conditions
    u(:,1)   = 0;
    u(:,end) = 0;
    u(1,:)   = 0;
    u(end,:) = 0;

    % Apply Neumann boundary conditions
    for j = padding+1:ny-padding
        u(j, padding+1) = u(j, padding+2) - dx * Mx(j,padding+1);          % Left boundary
        u(j, nx-padding) = u(j, nx-padding-1) + dx * Mx(j,nx-padding);     % Right boundary
    end
    
    for i = padding+1:nx-padding
        u(padding+1, i) = u(padding+2, i) - dx * My(padding+1,i);     % Bottom boundary
        u(ny-padding, i) = u(ny-padding-1, i) + dx * My(ny-padding,i); % Top boundary
    end

    % Check for convergence
    if max(max(abs(u - u_old))) < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
end

% Plot the solution
% x = 0:dx:Lx; y = 0:dx:Ly; [X,Y] = meshgrid(x,y);
x = 0:dx:Lx_container; y = 0:dx:Ly_container; [X,Y] = meshgrid(x,y);

figure;
surf(X,Y,u,'EdgeColor','interp')
% plot3(X,Y,u,'b.-')
title('Potential of demag field');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');

% Calculation of magnetic field Hd: Hd = -grad u
% dV/dx: u(j,i+1) - u(j,i)
% dV/dy: u(j+1,i) - u(j,i)

Hd = zeros(nx*ny,2);
% center nodes
for i = 2:nx-1
    for j = 2:ny-1
        Hd((i-1)*ny + j,1) = -(u(j,i+1)-u(j,i))/dx; % -dV/dx
        Hd((i-1)*ny + j,2) = -(u(j+1,i)-u(j,i))/dx; % -dV/dy
    end
end

% border nodes
for i = 2:nx-1
    
    Hd((i-1)*ny + 1,1) = -(u(1,i+1) - u(1,i))/dx;    % -dV/dx ; at y = 0
    Hd((i-1)*ny + 1,2) = -(u(2,i) - u(1,i))/dx;      % -dV/dy ; at y = 0 

    Hd((i-1)*ny + ny,1) = -(u(ny,i+1) - u(ny,i))/dx;  % -dV/dx ; at y = Ly
    Hd((i-1)*ny + ny,2) = -(u(ny,i) - u(ny-1,i))/dx;  % -dV/dy ; at y = Ly

end

for j = 1:ny-1

    Hd((1-1)*ny + j,1) = -(u(j,2) - u(j,1))/dx;      % -dV/dx ; at x = 0
    Hd((1-1)*ny + j,2) = -(u(j+1,1) - u(j,1))/dx;    % -dV/dy ; at x = 0 

    Hd((nx-1)*ny + j,1) = -(u(j,nx) - u(j,nx-1))/dx;  % -dV/dx ; at x = Lx
    Hd((nx-1)*ny + j,2) = -(u(j+1,nx) - u(j,nx))/dx;  % -dV/dy ; at x = Lx

end

% Corner nodes
Hd(1,1) = -(u(1,2)-u(1,1))/dx; % -dV/dx at (0,0)
Hd(1,2) = -(u(2,1)-u(1,1))/dx; % -dV/dy at (0,0)

Hd((nx-1)*ny + 1,1) = -(u(1,nx)-u(1,nx-1))/dx; % -dV/dx at (Lx,0)
Hd((nx-1)*ny + 1,2) = -(u(2,nx)-u(1,nx))/dx;   % -dV/dy at (Lx,0)

Hd(ny,1) = -(u(ny,2)-u(ny,1))/dx; % -dV/dx at (0,Ly)
Hd(ny,2) = -(u(ny,1)-u(ny-1,1))/dx; % -dV/dy at (0,Ly)

Hd((nx-1)*ny + ny,1) = -(u(ny,nx)-u(ny,nx-1))/dx; % -dV/dx at (Lx,Ly)
Hd((nx-1)*ny + ny,2) = -(u(ny,nx)-u(ny-1,nx))/dx; % -dV/dy at (Lx,Ly)

% Plot vector field
figure
quiver(X(:),Y(:),Hd(:,1),Hd(:,2))

% Plot |Hd|
Hd_norm = sqrt(sum(Hd.^2,2));
figure
scatter(X(:),Y(:),10,Hd_norm)
