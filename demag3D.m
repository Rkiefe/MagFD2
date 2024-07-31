%{
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

% H external
Hext = [1,0,0];

% Iterative solver parameters
max_iter = 1e5;
tol = 1e-6;

% Dimensions and spacing
Lx = 1;
Ly = 1;
Lz = 1;

padding = 15;

dx = 0.1;

Lx_container = Lx + (2*padding)*dx; % Add a padding
Ly_container = Ly + (2*padding)*dx; % Add a padding
Lz_container = Lz + (2*padding)*dx; % Add a padding

% Calculate number of grid points
nx = Lx/dx + 2*padding + 1; disp("nx: "+nx)
ny = Ly/dx + 2*padding + 1; disp("ny: "+ny)
nz = Lz/dx + 2*padding + 1; disp("nz: "+nz)

% Define the Magnetization
% Mx = zeros(ny-2*padding,nx-2*padding); % Only of the magnetic material
Mx = zeros(nz,ny,nx);
My = zeros(nz,ny,nx);
Mz = zeros(nz,ny,nx);

Mx(padding+1:nz-padding,padding+1:ny-padding,padding+1:nx-padding) = Hext(1);

x = 0:dx:Lx_container;
y = 0:dx:Ly_container;
z = 0:dx:Lz_container;

tic
u = zeros(nz,ny,nx);
for iter = 1:max_iter
    u_old = u;

    % Laplace u = div M
    for i = 2:nx-1
        for j = 2:ny-1
            for l = 2:nz-1
                % Div M, centered diff
                F = (Mx(l,j,i+1)-Mx(l,j,i-1) + ...
                     My(l,j+1,i)-My(l,j-1,i) + ...
                     Mz(l+1,j,i)-Mz(l-1,j,i))/(2*dx);

                u(l,j,i) = (u(l,j,i+1) + u(l,j,i-1) + ...
                            u(l,j+1,i) + u(l,j-1,i) + ...
                            u(l+1,j,i) + u(l-1,j,i) - F*dx^2)/6;
            end
        end
    end

    % Neuman border conditions
    % z-faces
    u(padding+1,:,:) = u(padding+2,:,:) -dx*Mz(padding+1,:,:);
    u(nz-padding,:,:) = u(nz-padding-1,:,:) +dx*Mz(nz-padding,:,:);

    % y-faces
    u(:,padding+1,:) = u(:,padding+2,:) -dx*My(:,padding+1,:);
    u(:,ny-padding,:) = u(:,ny-padding-1,:) +dx*My(:,ny-padding,:);

    % x-faces
    u(:,:,padding+1) = u(:,:,padding+2) -dx*Mx(:,:,padding+1);
    u(:,:,nx-padding) = u(:,:,nx-padding-1) +dx*Mx(:,:,nx-padding);
    
    % Dirichlet border conditions
    u(:,:,1) = 0;
    u(:,:,end) = 0;
    u(:,1,:) = 0;
    u(:,end,:) = 0;
    u(1,:,:) = 0;
    u(end,:,:) = 0;


    if max(abs(u(:)-u_old(:))) < tol
        fprintf('Converged after %d iterations.\n', iter);
        break
    end
end
clear u_old
toc

[Hz,Hy,Hx] = gradient(-u,dx);

[X,Y,Z] = meshgrid(x(padding+1:nx-padding),y(padding+1:ny-padding),z(padding+1:nz-padding));
px = Hx(padding+1:nz-padding,padding+1:ny-padding,padding+1:nx-padding);
py = Hy(padding+1:nz-padding,padding+1:ny-padding,padding+1:nx-padding);
pz = Hz(padding+1:nz-padding,padding+1:ny-padding,padding+1:nx-padding);

quiver3(X(:),Y(:),Z(:),px(:),py(:),pz(:))
