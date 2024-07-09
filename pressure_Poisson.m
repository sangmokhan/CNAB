function p = pressure_Poisson(b,Nx,Ny,dx,dy)
    persistent Laplacian
    
    if isempty(Laplacian)
        Delta = (dx+dy)/2;
        
        % Second derivative matrix for x
        Dxx = diag(-4*ones(Nx,1))+diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1);
        Dxx = Dxx/Delta^2;
        
        % Second derivative matrix for y
        Dyy = diag(ones(Ny-1,1),1)+diag(ones(Ny-1,1),-1);
        Dyy = Dyy/Delta^2;
        
        % Row-wise concatenation
        Dxx_kronecker = kron(eye(Ny,Ny),Dxx);
        Dyy_kronecker = kron(Dyy,eye(Nx,Nx));
        Laplacian = Dxx_kronecker+Dyy_kronecker;
        
        % Pinning the pressure at a point.
        Laplacian(1,:) = 0;
        Laplacian(1,1) = 1;

        % Matrix decomposition for efficient solving of system of linear equations.
        Laplacian = decomposition(Laplacian);
    end
    
    b = b(:);
    % Pinning the pressure at a point.
    b(1) = 0;

    p_vec = Laplacian\b;
    p = reshape(p_vec,Nx,Ny);
end