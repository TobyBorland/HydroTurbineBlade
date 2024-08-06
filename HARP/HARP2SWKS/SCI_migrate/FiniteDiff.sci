function [DFDX] = FiniteDiff(X,F)

// h is the stepsize
// h = X(2) - X(1);

dfdx = zeros(length(X),1); // Pre-allocate for speed

//  // Forward Difference
//  for n = 1
//      dfdx(n,1) = (F(n+1) - F(n)) / h;
//  end
//  
//  // Central Difference
//  for n = 2
//      dfdx(n,1) = (F(n+1) - F(n-1)) / (2*h);
//  end
//  
//  // Five-Stencil Finite Difference
//  for n = 3:(length(X)-2)
//      dfdx(n,1) = (F(n-2) - 8*F(n-1) + 8*F(n+1) -F(n+2))/(12*h);
//  end
//  
//  // Central Difference
//  for n = length(X) - 1
//      dfdx(n,1) = (F(n+1) - F(n-1)) / (2*h);
//  end
//  
//  // Backward Difference
//  for n = length(X)
//      dfdx(n,1) = (F(n) - F(n-1)) / h;
//  end

// Forward Difference
for n = 1:length(X)-1
    h = X(n+1) - X(n);
    dfdx(n,1) = (F(n+1) - F(n)) / h;
end

// Backward Difference
for n = length(X)
    h = X(n) - X(n-1);
    dfdx(n,1) = (F(n) - F(n-1)) / h;
end

DFDX = dfdx;
endfunction
