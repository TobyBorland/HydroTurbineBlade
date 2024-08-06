function [F] = Euler(X,DFDX,Fo)
// X and DFDX are equally spaced vectors
// Fo is a scalar, the initial value of the unknown function F(X)

// h is the stepsize
h = X(2) - X(1);

f = zeros(length(X),1);  // Preallocating for speed
f(1) = Fo; // Initial value of the unknown function f(X)

for n = 1:(length(X)-1)
    f(n+1) = f(n) + h*DFDX(n);
end

F = f;
endfunction
