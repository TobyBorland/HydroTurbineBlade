function [y,%poly] = CubeFit2(x1,x2,fx1,fx2,dfx1,dfx2,X)

// Output variables initialisation (not found in input variables)
y=[];
%poly=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);



A = [3*(x1^2),2*x1,1,0;
     x1^3,x1^2,x1,1;
     3*(x2^2),2*x2,1,0;
     x2^3,x2^2,x2,1];
b = [dfx1;fx1;dfx2;fx2];
poly = A\b;

//x=x1:(x2-x1)/100:x2;
y = mtlb_a(mtlb_a(mtlb_a(%poly(1)*(X .^3),%poly(2)*(X .^2)),%poly(3)*X),%poly(4));
//plot(x,y);
endfunction
