function [X] = Bezier(C,N)

// Bezier.m plots N points on the Bezier curve with control points in C
// credit for this algorithm goes to Dr. Craig L. Zirbel at Bowling Green
// State University (zirbel@bgsu.edu)

// Inputs: C is expected to be a two column vector, N is a scalar

// Transpose C:
C = C';

n = length(C(1,:));            % number of control points
a = 0:(1/(N-1)):1;             % parameter values along the Bezier curve
I = (0:(n-1))'*ones(1,N);      % matrix of index values
P = ones(n,1)*a;               % matrix of probabilities
A = binopdf(I,n-1,P);          % weights for control points
X = [C*A]';                    % points on the Bezier curve



