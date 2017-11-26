function y = ex2(f,dt,tend,yi) 

%y is the output 
%f is the function, depending on which method we are using 
%dt is the timestep size, for example 1 second, 
%tend is the end time 
%yi is y0, or y initial

t = 0 ; %current time "t", each loop t = t+dt 
y = zeros(tend/dt, 1); %y is a vector (1 column with zeros)
counter = 1; %the entry on the vector and the times loops are being run

while t <= tend
    if t == 0
        y(counter)=yi;  
        counter = counter +1; 
        t = t+dt; 
    else
    y(counter) = y(counter-1) + f(y(counter-1))*dt; 
    counter = counter +1; 
    t = t + dt; 
    end
end
end


%un finished product, only tested against Euler's method, 
%but if this function is method in-dependent, then any method should work. 