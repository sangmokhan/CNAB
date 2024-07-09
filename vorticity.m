function [omega] = vorticity(Nx,Ny,dx,dy,u,v)
    omega = zeros(Nx,Ny);

    for i = 1:Nx
        for j = 1:Ny
            omega(i,j) = (v(i+1,j)-v(i,j))/dx - (u(i,j+1)-u(i,j))/dy;
        end
    end
    
    % Boundaries
    omega(1,1) = 1/2 * (omega(1,2) + omega(2,1));
    omega(1,Ny) = 1/2 * (omega(2,Ny) + omega(1,Ny-1));
    omega(Nx,1) = 1/2 * (omega(Nx-1,1) + omega(Nx,2));
    omega(Nx,Ny) = 1/2 * (omega(Nx-1,Ny) + omega(Nx,Ny-1));
end