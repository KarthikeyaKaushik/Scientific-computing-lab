function A = A_gen(Nx,Ny)
%A_GEN Generating the A matrix using Nx and Ny as paramters
hx = 1/(Nx+1);
hy = 1/(Ny+1);
hx2 = hx^2;
hy2 = hy^2;
%A = toeplitz([-4 1 zeros(1,Ny-2) 1 zeros(1,Nx*Ny-Ny-1)]);
A = toeplitz([-2*(1/hx2+1/hy2) 1/hy2 zeros(1,Ny-2) 1/hx2 zeros(1,Nx*Ny-Ny-1)]);
end

