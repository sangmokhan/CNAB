function [div,courant] = divCourant(Nx,Ny,dt,dx,dy,u,v)
    div = zeros(Nx+1,Ny+1);
    
    for i = 2:Nx
        for j = 2:Ny
            % continuity
            div(i,j) = (u(i,j)-u(i-1,j))/dx + (v(i,j)-v(i,j-1))/dy;
        end
    end

    % Courant number calculation
    courant = 0;
    for i = 2:Nx
        for j = 2:Ny
            us = 1/2 * (u(i,j)-u(i-1,j));
            vs = 1/2 * (v(i,j)-v(i,j-1));
            courant_new = abs(us/dx) + abs(vs/dy);
            if courant_new > courant
                courant = courant_new;
            end
        end
    end

    courant = courant * dt;
end