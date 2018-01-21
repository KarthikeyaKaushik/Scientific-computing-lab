function Xout = implicit_Euler_2D(X, dt, nx, ny)
hx = 1/(nx+1);
hy = 1/(ny+1);
Xout = zeros(nx+2, ny+2);
a = dt/hx^2;
b = dt/hy^2;
c = 1 + 2*(dt/hx^2) + 2*(dt/hy^2);
residual = 5;
N = nx*ny;
n=0;
error = zeros(nx, ny);
while residual>10^-6 && n < 10000  %arbitrary number of iteration to get out
                                   %of the loop in unstable case
    %Iteration for Gauss Seidel 
    % Equation for Implicit Euler:
    % Xout = X + dt*Tt    Tt is out heat equation with time n+1, X is at n
    for i = 2:nx+1
        for j = 2:ny+1
        % by rearranging the implicit euler equation we get this
            Xout(i, j) = (X(i,j))/c + (a/c)*(Xout(i+1, j) + Xout(i-1, j))+...
                    (b/c)*(Xout(i, j+1) + Xout(i, j-1));
        end
    end
    % Calculation residual at every point. It Should as close to zero as
    % possible
    for i = 2:nx+1
        for j=2:ny+1
            % Moving every term in the heat equation by making it implicit
            % on the R.H.S we get the following and it should be close to
            % zeros
            error(i-1, j-1) = X(i, j) + (a)*(Xout(i+1, j) + Xout(i-1, j))+...
                (b)*(Xout(i, j+1) + Xout(i, j-1)) - (c)*(Xout(i, j)); 
        end
    end
    residual = sqrt((1/N)*(sum(sum(error.^2))));
    n=n+1;
end
end
