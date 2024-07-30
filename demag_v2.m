%{
    This code solves the poisson equation for the demag field
    using Dirichlet and Neuman border conditions, meshing both 
    the refrigerant and the air around it.

    Dif. eq.:
    Nabla \dot Nabla V = Nabla dot M 

    Border conditions for the refrigerant (Neumann)
    (nabla V_int - nabla V_ext) dot n = M dot n

    Border conditions for the container (Dirichlet)
    V = 0
%}

clear
close all
clc

Hext = [1,0];

% Plot results
plotting = true;

% Dimensions and spacing
Lx = 1;
Ly = 0.5;
dx = 0.01;
dx_out = 0.05;

one_over_step = 1/(1/dx + 1/dx_out);

padding = 10;

Lx_container = Lx + (2*padding)*dx_out; % Add a padding
Ly_container = Ly + (2*padding)*dx_out; % Add a padding

% For the plot:
pgon = polyshape([padding*dx_out , padding*dx_out+Lx, padding*dx_out+Lx, padding*dx_out], [padding*dx_out , padding*dx_out, padding*dx_out+Ly , padding*dx_out+Ly]);


% Calculate number of grid points
nx = Lx/dx + 2*padding + 1; disp("nx: "+nx)
ny = Ly/dx + 2*padding + 1; disp("ny: "+ny)

x = [...
        0:dx_out:padding*dx_out, ...
        padding*dx_out + dx:dx:(Lx_container-padding*dx_out-dx), ...
        Lx_container-padding*dx_out:dx_out:Lx_container
    ];


y = [...
        0:dx_out:padding*dx_out, ...
        padding*dx_out + dx:dx:(Ly_container-padding*dx_out-dx), ...
        Ly_container-padding*dx_out:dx_out:Ly_container
    ];


% Initialize u
u = zeros(ny, nx);

% Define the Magnetization
% Mx = zeros(ny-2*padding,nx-2*padding); % Only of the magnetic material
Mx = zeros(ny,nx);
My = zeros(ny,nx);

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

tic
% Iterative solver parameters
max_iter = 1e5;
tol = 1e-6;

% Perform iterations to solve the interior points
for iter = 1:max_iter
    u_old = u;

    % Magnet interior nodes
    for i = padding+2:nx-padding-1        
        for j = padding+2:ny-padding-1
            
            F = (Mx(j,i) - Mx(j,i) + ...
                 My(j,i) - My(j,i))/dx;

            u(j,i) = 0.25 * (u_old(j+1,i) + u_old(j-1,i) + ...
                             u_old(j,i+1) + u_old(j,i-1) - ...
                             dx^2 * F);
        end % end of y loop
    end % end of x loop

    % Exterior nodes
    for j = 2:ny-1
        for i = 2:padding
            u(j,i) = 0.25 * (u_old(j+1,i) + u_old(j-1,i) + ...
                             u_old(j,i+1) + u_old(j,i-1));
        end

        for i = nx-padding+1:nx-1
            u(j,i) = 0.25 * (u_old(j+1,i) + u_old(j-1,i) + ...
                             u_old(j,i+1) + u_old(j,i-1));
        end
    end
    for i = 2:nx-1
        for j = 2:padding
            u(j,i) = 0.25 * (u_old(j+1,i) + u_old(j-1,i) + ...
                             u_old(j,i+1) + u_old(j,i-1));
        end

        for j = ny-padding+1:ny-1
            u(j,i) = 0.25 * (u_old(j+1,i) + u_old(j-1,i) + ...
                             u_old(j,i+1) + u_old(j,i-1));
        end
    end

    % Neumann border conditions
    % Border y = 0 of the magnet
    j = padding+1;
    for i = padding+1:nx-padding
        u(j,i) = (-My(j,i) + u(j+1,i)/dx + u(j-1,i)/dx_out )*one_over_step;
    end

    % Border y = Ly of the magnet
    j = ny-padding;
    for i = padding+1:nx-padding
        u(j,i) = (My(j,i) + u(j-1,i)/dx + u(j+1,i)/dx_out )*one_over_step;
    end

    % Border x = 0 of the magnet
    i = padding+1;
    for j = padding+1:ny-padding
        u(j,i) = (-Mx(j,i) + u(j,i+1)/dx + u(j,i-1)/dx_out )*one_over_step;
    end

    % Border x = Lx of the magnet
    i = nx-padding;
    for j = padding+1:ny-padding
        u(j,i) = (Mx(j,i) + u(j,i-1)/dx + u(j,i+1)/dx_out )*one_over_step;
    end

    % Dirichlet border conditions
    u(:,1) = 0;
    u(:,end) = 0;
    u(1,:) = 0;
    u(end,:) = 0;

    
    % Check for convergence
    if max(max(abs(u - u_old))) < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
end
clear u_old
toc

if plotting % Plot the solution
    [X,Y] = meshgrid(x,y);

    figure;
    surf(X,Y,u,'EdgeColor','interp')
    
    title('Potential of demag field');
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
end

% Calculation of magnetic field Hd: Hd = -grad u, only for the interior nodes
Hd_x = zeros(ny,nx);
Hd_y = zeros(ny,nx);

for i = padding+1:nx-padding
    ix2 = i+1;
    ix1 = i;

    if i == nx-padding
        ix2 = i;
        ix1 = i-1;
    end

    for j = padding+1:ny-padding
        iy2 = j+1;
        iy1 = j;

        if j == ny-padding
            iy2 = j;
            iy1 = j-1;
        end

        Hd_x(j,i) = -(u(j,ix2)-u(j,ix1))/dx; % -dV/dx
        Hd_y(j,i) = -(u(iy2,i)-u(iy1,i))/dx; % -dV/dy
    
    end

end

Hd = [Hd_x(:),Hd_y(:)];

if plotting % Plot vector field
    figure
    quiver(X(:),Y(:),Hd_x(:),Hd_y(:))
    hold on

    pgon.plot("FaceAlpha",0,"EdgeColor","r")

    % Plot |Hd|
    % figure
    % scatter(X(:),Y(:),10,Hd_norm); clear X Y
    % surf(X,Y,sqrt(Hd_x.^2+Hd_y.^2),'EdgeColor','interp'); % clear X Y
    % hold on
    % pgon.plot("FaceAlpha",0,"EdgeColor","r")
    % title("|Hd|")
    % colorbar

end
