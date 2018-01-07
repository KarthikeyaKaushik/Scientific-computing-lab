function [Asparse,b,soln_sparse,t] = sparse_solver(Nx,Ny)
%SPARSE_SOLVER Summary of this function goes here
%   Detailed explanation goes here
tic;
hx = 1/(Nx+1);
hy = 1/(Ny+1);
hx2 = hx^2;
hy2 = hy^2;

b = zeros(Ny,Nx);
for i = 1:Ny
    for j = 1:Nx
        b(i,j) = -2*pi^2*sin(pi*i/(1+Ny))*sin(pi*j/(1+Nx));
    end
end
b = reshape(b,Nx*Ny,1);

A = A_gen(Nx,Ny);
Asparse = sparse(A);
soln_sparse = Asparse\b;
soln_sparse = reshape(soln_sparse,Ny,Nx);
soln_sparse = padarray(soln_sparse,[1 1],0);
    
t = toc;

end

