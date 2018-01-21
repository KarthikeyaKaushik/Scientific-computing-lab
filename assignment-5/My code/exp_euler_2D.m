function Xout = exp_euler_2D(X, dt, nx,ny)
hx = 1/(nx+1);
hy = 1/(ny+1);
% The equation implemented for Explicit Euler looks like:
% Xout(i,j) = X(i,j) + dt*Tt .... "Tt" is our heat equation
Xout = zeros(nx+2, ny+2);
a = dt/hx^2;
b = dt/hy^2;
c = 1 - 2*(dt/hx^2) - 2*(dt/hy^2);
% as this is just matrix-matrix addition/subtraction so it can be done in
% one line
Xout(2:end-1, 2:end-1) = a*(X(1:end-2, 2:end-1) + X(3:end, 2:end-1))+...
                         b*(X(2:end-1, 3:end) + X(2:end-1, 1:end-2))+...
                         c*(X(2:end-1, 2:end-1));
end