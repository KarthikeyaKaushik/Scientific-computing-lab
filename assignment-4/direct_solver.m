function [A,b,soln,t] = direct_solver(Nx,Ny)
% Using the direct method

tstart = tic;
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

soln = A\b;
soln = reshape(soln,Ny,Nx);
soln = padarray(soln,[1 1],0);
    
t = toc(tstart);


end

