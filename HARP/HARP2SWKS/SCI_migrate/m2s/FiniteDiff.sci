function [DFDX] = FiniteDiff(X,F)

// Output variables initialisation (not found in input variables)
DFDX=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


//h is the stepsize
//h = X(2) - X(1);

dfdx = zeros(max(size(X)),1);//Pre-allocate for speed

// %Forward Difference
// for n = 1
//     dfdx(n,1) = (F(n+1) - F(n)) / h;
// end
// 
// %Central Difference
// for n = 2
//     dfdx(n,1) = (F(n+1) - F(n-1)) / (2*h);
// end
// 
// %Five-Stencil Finite Difference
// for n = 3:(length(X)-2)
//     dfdx(n,1) = (F(n-2) - 8*F(n-1) + 8*F(n+1) -F(n+2))/(12*h);
// end
// 
// %Central Difference
// for n = length(X) - 1
//     dfdx(n,1) = (F(n+1) - F(n-1)) / (2*h);
// end
// 
// %Backward Difference
// for n = length(X)
//     dfdx(n,1) = (F(n) - F(n-1)) / h;
// end

//Forward Difference
for n = 1:max(size(X))-1
  h = mtlb_s(mtlb_e(X,n+1),mtlb_e(X,n));
  dfdx(n,1) = mtlb_s(mtlb_e(F,n+1),mtlb_e(F,n))/h;
end;

//Backward Difference
for n = max(size(X))
  h = mtlb_s(mtlb_e(X,n),mtlb_e(X,n-1));
  dfdx(n,1) = mtlb_s(mtlb_e(F,n),mtlb_e(F,n-1))/h;
end;

DFDX = dfdx;
endfunction
