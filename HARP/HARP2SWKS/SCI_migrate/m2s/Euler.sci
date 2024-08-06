function [F] = Euler(X,DFDX,Fo)

// Output variables initialisation (not found in input variables)
F=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//X and DFDX are equally spaced vectors
//Fo is a scalar, the initial value of the unknown function F(X)

//h is the stepsize
h = mtlb_s(mtlb_e(X,2),mtlb_e(X,1));

f = zeros(max(size(X)),1);//Preallocating for speed
f = mtlb_i(f,1,Fo);//Initial value of the unknown function f(X)

for n = 1:max(size(X))-1
  f = mtlb_i(f,n+1,mtlb_a(f(n),h*mtlb_e(DFDX,n)));
end;

F = f;
endfunction
