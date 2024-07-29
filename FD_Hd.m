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

Hext = [1,0];

% Plot results
plotting = true;

% Dimensions and spacing
Lx = 1;
Ly = 0.5;
dx = 0.01;

padding = 100;

Lx_container = Lx + (2*padding)*dx; % Add a padding
Ly_container = Ly + (2*padding)*dx; % Add a padding

% For the plot:
pgon = polyshape([padding*dx , padding*dx+Lx, padding*dx+Lx, padding*dx], [padding*dx , padding*dx, padding*dx+Ly , padding*dx+Ly]);


% Calculate number of grid points
nx = Lx_container / dx + 1;
ny = Ly_container / dx + 1;

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
            My(j,i) = 1;
        end
    end
end

% Iterative solver parameters
max_iter = 1e5;
tol = 1e-6;

% Perform iterations to solve the interior points
for iter = 1:max_iter
    u_old = u;

    for i = 1:nx
        
        i1 = i;
        i2 = i+1;

        for j = 1:ny

            j1 = j;
            j2 = j+1;

            % Div M | magnet nodes
            if i > padding+1 && i < nx-padding && j > padding+1 && j < ny-padding
                
                if i == nx-padding
                    i2 = i;
                    i1 = i-1;
                end

                if j == ny-padding
                    j2 = j;
                    j1 = j-1;
                end

                F = (Mx(j,i2) - Mx(j,i1) + ...
                     My(j2,i) - My(j1,i))/dx;

                u(j, i) = 0.25 * (u_old(j+1, i) + u_old(j-1, i) + ...
                                  u_old(j, i+1) + u_old(j, i-1) - ...
                                  dx^2 * F);
                
                continue

            else % Div M | Air
                F = 0;
            end

            % Dirichlet border conditions
            if i==1 || i==nx || j == 1 || j == ny
                u(j,i) = 0;

                continue
            end

            % Neumann border conditions
            % (i==nx-padding)

            if (i==padding+1) && ( (j>padding) &&  (j<ny-padding) )
                u(j,i) = 0.5*(-Mx(j,i)*dx + u(j,i-1) + u(j,i+1) );
                
                continue
            end

            if (i==nx-padding) && ( (j>padding) &&  (j<ny-padding) )
                u(j,i) = 0.5*(Mx(j,i)*dx + u(j,i-1) + u(j,i+1) );

                continue
            end

            if (j==padding+1) && ( (i>padding) &&  (i<nx-padding) )
                u(j,i) = 0.5*(-My(j,i)*dx + u(j-1,i) + u(j+1,i) );
                
                continue
            end

            if (j==ny-padding) && ( (i>padding) &&  (i<nx-padding) )
                u(j,i) = 0.5*(My(j,i)*dx + u(j-1,i) + u(j+1,i) );
                
                continue
            end

            % Exterior nodes
            u(j, i) = 0.25 * (u_old(j+1, i) + u_old(j-1, i) + ...
                                  u_old(j, i+1) + u_old(j, i-1));
        end
    end
    
    % Check for convergence
    if max(max(abs(u - u_old))) < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
end
clear u_old

if plotting % Plot the solution
    % x = 0:dx:Lx; y = 0:dx:Ly; [X,Y] = meshgrid(x,y);
    x = 0:dx:Lx_container; y = 0:dx:Ly_container; [X,Y] = meshgrid(x,y); clear x y

    figure;
    surf(X,Y,u,'EdgeColor','interp')
    % plot3(X,Y,u,'b.-')
    title('Potential of demag field');
    xlabel('x');
    ylabel('y');
    zlabel('u(x,y)');
end

% Calculation of magnetic field Hd: Hd = -grad u
Hd_x = zeros(ny,nx);
Hd_y = zeros(ny,nx);

for i = 1:nx

    % Check if it is in a border
    ix1 = i;
    ix2 = i+1;
    if ix2 > nx
        ix2 = nx;
        ix1 = nx-1;
    end

    for j = 1:ny

        % Check if it is in a border
        iy1 = j;
        iy2 = j+1;
        if iy2 > ny
            iy2 = ny;
            iy1 = ny-1;
        end

        Hd_x(j,i) = -(u(j,ix2)-u(j,ix1))/dx; % -dV/dx
        Hd_y(j,i) = -(u(iy2,i)-u(iy1,i))/dx; % -dV/dy
    end
end

Hd = [Hd_x(:),Hd_y(:)];

% |Hd|
Hd_norm = sqrt(sum(Hd.^2,2));

if plotting % Plot vector field
    figure
    quiver(X(:),Y(:),Hd(:,1),Hd(:,2))
    hold on

    pgon.plot("FaceAlpha",0,"EdgeColor","r")

    % Plot |Hd|
    figure
    % scatter(X(:),Y(:),10,Hd_norm); clear X Y
    surf(X,Y,sqrt(Hd_x.^2+Hd_y.^2),'EdgeColor','interp'); clear X Y
    hold on
    pgon.plot("FaceAlpha",0,"EdgeColor","r")

    title("|Hd|")
    colorbar
end

% Save potential for comparing:
% aux = zeros(ny-2*padding,nx-2*padding);
% for i = 1:nx-2*padding
%     for j = 1:ny-2*padding
%         aux(j,i) = u(padding+j,padding+i);
%     end
% end
% save("potential.mat","aux")
