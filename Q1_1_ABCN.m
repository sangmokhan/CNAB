clc;
clear;
clear pressure_Poisson ABCN;

% Parameters
Nx = 32;                   % Number of grid points in x
Ny = 32;                   % Number of grid points in y
H = 1;                      % Cavity width/height
Lx = H;                     % Domain width
Ly = H;                     % Domain height
U0 = 1;                     % Top lid velocity
Re = 100;                   % Reynolds number
dt = 1 * 10^(-3);           % Time step
max_num_iter = 6.25/dt;       % Max time steps
eps = 1 * 10^(-8);          % Convergence criterion

% Grid initialization (staggered grid)
dx = Lx/Nx;
dy = Ly/Ny;

% Variable preallocation
u = zeros(Nx+1,Ny+2);
v = zeros(Nx+2,Ny+1);
p = zeros(Nx,Ny);
u_prev = zeros(Nx+1,Ny+2);
v_prev = zeros(Nx+2,Ny+1);
u_star = zeros(Nx+1,Ny+2);
v_star = zeros(Nx+2,Ny+1);

% Projection method
t = 0;

for i = 1:max_num_iter
    t = t + dt; u_prev = u; v_prev = v;
    % 1. Intermediate velocity calculation

    % Boundary conditions
        % Top
    u(:,end) = 2*U0 - u(:,end-1);
    v(:,end) = 0;
        % Bottom
    u(:,1) = -u(:,2); 
    v(:,1) = 0;
        % Left
    u(1,:) = 0;
    v(1,:) = -v(2,:);
        % Right
    u(end,:) = 0;
    v(end,:) = -v(end-1,:);

    % Viscous term calculation (central difference scheme)
    viscous_u_x = (u(1:end-2,2:end-1) - 2*u(2:end-1,2:end-1) + u(3:end,2:end-1))/dx^2; % (Nx-1) * Ny
    viscous_u_y = (u(2:end-1,1:end-2) - 2*u(2:end-1,2:end-1) + u(2:end-1,3:end))/dy^2; % (Nx-1) * Ny
    viscous_v_x = (v(1:end-2,2:end-1) - 2*v(2:end-1,2:end-1) + v(3:end,2:end-1))/dx^2; % Nx * (Ny-1)
    viscous_v_y = (v(2:end-1,1:end-2) - 2*v(2:end-1,2:end-1) + v(2:end-1,3:end))/dy^2; % Nx * (Ny-1)
    
    % Boundary condition on viscous term
    viscous_u_BC = zeros(size(viscous_u_y));
    viscous_u_BC(:,end) = 2*U0/dy^2;
%         viscous_v_BC = zeros(size(viscous_v_y));

    % Convective term calculation
        % Velocity interpolation for cell center & corner (linear interpolation)
    u_center = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
    u_corner = (u(:,1:end-1)+u(:,2:end))/2;
    v_center = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
    v_corner = (v(1:end-1,:)+v(2:end,:))/2;
    
    u_center_sq = u_center.^2;
    uv_corners = u_corner.*v_corner;
    v_center_sq = v_center.^2;
    
    convective_u = (u_center_sq(2:end,:) - u_center_sq(1:end-1,:))/dx + (uv_corners(2:end-1,2:end) - uv_corners(2:end-1,1:end-1))/dy;
    convective_v = (v_center_sq(:,2:end) - v_center_sq(:,1:end-1))/dy + (uv_corners(2:end,2:end-1) - uv_corners(1:end-1,2:end-1))/dx;
        
    if i == 1
        % Forward Euler method in time
        u_star(2:end-1,2:end-1) = u(2:end-1,2:end-1) + dt * (-convective_u + (viscous_u_x+viscous_u_y)/Re);
        v_star(2:end-1,2:end-1) = v(2:end-1,2:end-1) + dt * (-convective_v + (viscous_v_x+viscous_v_y)/Re);
    else
        % Adams-Bashforth method for convection, Crank-Nicolson for viscous
        b = u(2:end-1,2:end-1) + dt * (1/(2*Re)*(viscous_u_x+viscous_u_y+viscous_u_BC)-(3/2)*convective_u+(1/2)*convective_u_prev);
        u_star(2:end-1,2:end-1) = ABCN(u,b,dx,dy,dt,Re);

        b = v(2:end-1,2:end-1) + dt * (1/(2*Re)*(viscous_v_x+viscous_v_y)-(3/2)*convective_v+(1/2)*convective_v_prev);
        v_star(2:end-1,2:end-1) = ABCN(v,b,dx,dy,dt,Re);
    end
    
    % Storing previous convective terms
    convective_u_prev = convective_u; convective_v_prev = convective_v;

    % 2. Projection

    % RHS of pressure Poisson equation
    b = (1/dt) * ((u_star(2:end,2:end-1)-u_star(1:end-1,2:end-1))/dx + (v_star(2:end-1,2:end)-v_star(2:end-1,1:end-1))/dy);
    
    % Solve pressure Poisson equation
    p = pressure_Poisson(b,Nx,Ny,dx,dy);
    
    % Correction of intermediate velocity
    u(2:end-1,2:end-1) = u_star(2:end-1,2:end-1) - dt * (p(2:end,:)-p(1:end-1,:))/dx;
    v(2:end-1,2:end-1) = v_star(2:end-1,2:end-1) - dt * (p(:,2:end)-p(:,1:end-1))/dy;

    % Boundary conditions
        % Top
    u(:,end) = 2*U0 - u(:,end-1);
    v(:,end) = 0;
        % Bottom
    u(:,1) = -u(:,2); 
    v(:,1) = 0;
        % Left
    u(1,:) = 0;
    v(1,:) = -v(2,:);
        % Right
    u(end,:) = 0;
    v(end,:) = -v(end-1,:);

    % Convergence criterion using Courant number
    [div,courant] = divCourant(Nx,Ny,dt,dx,dy,u,v);

    [cmax] = converge(u,u_prev,v,v_prev);
    if mod(i,1000) == 0
        disp(['Iter = ',int2str(i),'; Courant = ',num2str(courant),'; max(change) = ',num2str(cmax)]);
    end
    if cmax < eps
        break;
    end
end

disp(' ');
disp(['Total iterations = ',int2str(i),'; Total time = ',num2str(t),'; final(Courant) = ',num2str(courant)]);

% Velocity at cell center
u_center = (u(1:end-1,2:end-1) + u(2:end,2:end-1))/2;
v_center = (v(2:end-1,1:end-1) + v(2:end-1,2:end))/2;

% Computing stream function and vorticity
[psi] = stream_function(Nx,Ny,dx,v);
[omega] = vorticity(Nx,Ny,dx,dy,u,v);

% Plot results
x_center = ((1:Nx)-0.5) * Lx/Nx;
y_center = ((1:Ny)-0.5) * Ly/Ny;
plot_results(x_center,y_center,u_center,v_center,psi,omega,p);

%% Model validation

x = (0:Nx)*Lx/Nx;
y = (0:Ny)*Ly/Ny;
x_validation = x([1 9 10 11 13 21 30 31 65 104 111 117 122 123 124 125 129]);
y_validation = y([1 8 9 10 14 23 37 59 65 80 95 110 123 124 125 126 129]);
u_vel_validation = (u(ceil(end/2),[1 8 9 10 14 23 37 59 65 80 95 110 123 124 125 126 129]) + u(ceil(end/2),[1 8 9 10 14 23 37 59 65 80 95 110 123 124 125 126 129]+1))/2;
v_vel_validation = (v([1 9 10 11 13 21 30 31 65 104 111 117 122 123 124 125 129],ceil(end/2)) + v([1 9 10 11 13 21 30 31 65 104 111 117 122 123 124 125 129]+1,ceil(end/2)))/2;
omega_validation = abs(omega([8 16 24 32 40 48 56 64 72 80 88 96 104 112 120],end));
