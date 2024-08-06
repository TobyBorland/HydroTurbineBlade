function [X] = Bezier(C,N)

// Output variables initialisation (not found in input variables)
X=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


// Bezier.m plots N points on the Bezier curve with control points in C
// credit for this algorithm goes to Dr. Craig L. Zirbel at Bowling Green
// State University (zirbel@bgsu.edu)

//Inputs: C is expected to be a two column vector, N is a scalar

//Transpose C:
C = mtlb_t(C);

n = max(size(C(1,:)));// number of control points
a = mtlb_imp(0,1/mtlb_s(N,1),1);// parameter values along the Bezier curve
I = (0:n-1)'*ones(1,N);// matrix of index values
P = ones(n,1)*a;// matrix of probabilities
// !! L.16: Matlab toolbox(es) function binopdf not converted, original calling sequence used
A = binopdf(I,n-1,P);// weights for control points
X = (C*A)';// points on the Bezier curve



endfunction
