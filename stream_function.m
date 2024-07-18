% SangMok Han, Turbulence Lab, Yonsei University, August 2023
function [psi] = stream_function(Nx,Ny,dx,v)
    psi = zeros(Nx,Ny);

    for i = 2:Nx-1
        psi(i,2:Ny-1) = psi(i-1,2:Ny-1) - dx*v(i,2:Ny-1);
    end
end