clc;
clear all;
close all;
Nx = 7;
Ny = 7;
hx = 1/(Nx+1);
hy = 1/(Ny+1);
hx2 = hx^2;
hy2 = hy^2;


A = A_gen(Nx,Ny);
b = zeros(Ny,Nx);

for i = 1:Ny
    for j = 1:Nx
        b(i,j) = -2*pi^2*sin(pi*i/(1+Ny))*sin(pi*j/(1+Nx));
    end
end

b = reshape(b,Nx*Ny,1);
soln = A\b;
soln = reshape(soln,Ny,Nx);
soln = padarray(soln,[1 1],0);

%% For Gauss Seidel solver :

b = reshape(b,Ny,Nx);
bnew = padarray(b,[1 1],0);
b = reshape(b,Nx*Ny,1);
T = zeros(Ny+2,Ny+2);
residual = 10;
while(residual > 10^-4)
    for i = 2:Ny+1
        for j = 2:Nx+1
            T(i,j) = ((T(i-1,j) + T(i+1,j))/hx2 + (T(i,j+1) + T(i,j-1))/hy2 - bnew(i,j))/(2*(1/hx2+1/hy2));
        end
    end
    error = zeros(Ny,Nx);
    for i = 2:Ny+1
        for j = 2:Nx+1
            error(i-1,j-1) = T(i,j)*(-2*(1/hx2+1/hy2)) + ((T(i-1,j) + T(i+1,j))/hx2 + (T(i,j+1) + T(i,j-1))/hy2);
        end
    end
    residual = sqrt((1/Nx*Ny)*sum((b - reshape(error,Nx*Ny,1)).^2));
end




