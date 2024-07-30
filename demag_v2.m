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
% close all
clc

% Simulation settings
max_iter = 1e5;
tol = 1e-6;
Hext = [1,0];

% Dimensions and spacing
Lx = 1;
Ly = 0.5;
dx = 0.01;
dx_out = 0.05;

padding = 10;

Lx_container = Lx + (2*padding)*dx_out; % Add a padding
Ly_container = Ly + (2*padding)*dx_out; % Add a padding

% Calculate number of grid points
nx = Lx/dx + 2*padding + 1; disp("nx: "+nx)
ny = Ly/dx + 2*padding + 1; disp("ny: "+ny)

% For the plot:
pgon = polyshape([padding*dx_out , padding*dx_out+Lx, padding*dx_out+Lx, padding*dx_out], [padding*dx_out , padding*dx_out, padding*dx_out+Ly , padding*dx_out+Ly]);

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

% Define the Magnetization
Mx = zeros(Ly/dx + 1,Lx/dx + 1);
My = zeros(Ly/dx + 1,Lx/dx + 1);;

Mx(:) = Hext(1); 
My(:) = Hext(2);

% for i = padding+1:nx-padding
%     ip = i-padding;
    
%     for j = padding+1:ny-padding
%         jp = j - padding;    

%         Mx(jp,ip) = dx*(ip-1); % Hext(1)
%         My(jp,ip) = Hext(2);
%     end
% end

tic
[Hd_x,Hd_y] =  Hdemag(Mx,My, ...
                      Lx,Ly, ...
                      nx,ny, ...
                      dx,dx_out, ...
                      padding,...
                      max_iter,...
                      tol);
toc

% Plot
[X,Y] = meshgrid(x(padding+1:nx-padding),y(padding+1:ny-padding));
figure
quiver(X(:),Y(:),Hd_x(:),Hd_y(:))
hold on

pgon.plot("FaceAlpha",0,"EdgeColor","r")


function [Hd_x,Hd_y] = Hdemag(Mx,My, ...
                              Lx,Ly, ...
                              nx,ny, ...
                              dx,dx_out, ...
                              padding,...
                              max_iter,...
                              tol)

    one_over_step = 1/(1/dx + 1/dx_out);

    % Initialize u
    u = zeros(ny, nx);

    % Perform iterations to solve the interior points
    for iter = 1:max_iter
        u_old = u;

        % Magnet interior nodes
        for i = padding+2:nx-padding-1        
            for j = padding+2:ny-padding-1

                % Div M
                F = (Mx(j-padding,i-padding) - Mx(j-padding,i-padding-1) + ...
                     My(j-padding,i-padding) - My(j-padding-1,i-padding))/dx;
                
                u(j,i) = 0.25 * (u_old(j+1,i) + u_old(j-1,i) + ...
                                 u_old(j,i+1) + u_old(j,i-1) - ...
                                 dx^2 * F);
            end
        end

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
            u(j,i) = (-My(j-padding,i-padding) + u(j+1,i)/dx + u(j-1,i)/dx_out )*one_over_step;
        end

        % Border y = Ly of the magnet
        j = ny-padding;
        for i = padding+1:nx-padding
            u(j,i) = (My(j-padding,i-padding) + u(j-1,i)/dx + u(j+1,i)/dx_out )*one_over_step;
        end

        % Border x = 0 of the magnet
        i = padding+1;
        for j = padding+1:ny-padding
            u(j,i) = (-Mx(j-padding,i-padding) + u(j,i+1)/dx + u(j,i-1)/dx_out )*one_over_step;
        end

        % Border x = Lx of the magnet
        i = nx-padding;
        for j = padding+1:ny-padding
            u(j,i) = (Mx(j-padding,i-padding) + u(j,i-1)/dx + u(j,i+1)/dx_out )*one_over_step;
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
    
    end % end of potential calculation
    clear u_old
    
    % Calculation of magnetic field Hd: Hd = -grad u, only for the interior nodes
    Hd_x = zeros(Ly/dx + 1,Lx/dx + 1);
    Hd_y = zeros(Ly/dx + 1,Lx/dx + 1);;
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

            Hd_x(j-padding,i-padding) = -(u(j,ix2)-u(j,ix1))/dx; % -dV/dx
            Hd_y(j-padding,i-padding) = -(u(iy2,i)-u(iy1,i))/dx; % -dV/dy
        end
    end % End of Hd calculation

end % End of Hdemag function
