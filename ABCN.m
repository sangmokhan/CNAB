function u_star = ABCN(u,b,dx,dy,dt,Re)
    persistent A B
    [Nx, Ny] = size(u(2:end-1,2:end-1));

    if isempty(A) && Nx < Ny
        % Second derivative matrix for x
        Dxx = diag(-2*ones(Nx,1))+diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1);
        Dxx = Dxx/dx^2;
        
        % Second derivative matrix for y
        Dyy = diag(-2*ones(Ny,1))+diag(ones(Ny-1,1),1)+diag(ones(Ny-1,1),-1);
        Dyy(1,1) = Dyy(1,1)-1; Dyy(end,end) = Dyy(end,end)-1;   % Boundary condition
        Dyy = Dyy/dy^2;
        
        % Row-wise concatenation
        Dxx_kronecker = kron(eye(Ny,Ny),Dxx);
        Dyy_kronecker = kron(Dyy,eye(Nx,Nx));
        Laplacian = Dxx_kronecker+Dyy_kronecker;

        A = eye(size(Laplacian)) - dt/(2*Re)*Laplacian;

        % Matrix decomposition for efficient solving of system of linear equations.
        A = decomposition(A);
    end

    if isempty(B) && Nx > Ny
        % Second derivative matrix for x
        Dxx = diag(-2*ones(Nx,1))+diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1);
        Dxx(1,1) = Dxx(1,1)-1; Dxx(end,end) = Dxx(end,end)-1;   % Boundary condition
        Dxx = Dxx/dx^2;
        
        % Second derivative matrix for y
        Dyy = diag(-2*ones(Ny,1))+diag(ones(Ny-1,1),1)+diag(ones(Ny-1,1),-1);
        Dyy = Dyy/dy^2;
        
        % Row-wise concatenation
        Dxx_kronecker = kron(eye(Ny,Ny),Dxx);
        Dyy_kronecker = kron(Dyy,eye(Nx,Nx));
        Laplacian = Dxx_kronecker+Dyy_kronecker;

        B = eye(size(Laplacian)) - dt/(2*Re)*Laplacian;

        % Matrix decomposition for efficient solving of system of linear equations.
        B = decomposition(B);
    end
    
    if Nx < Ny
        u_star = A\b(:);
    elseif Nx > Ny
        u_star = B\b(:);
    end
    u_star = reshape(u_star,Nx,Ny);
end