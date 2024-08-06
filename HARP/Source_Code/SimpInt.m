function [Output] = SimpInt(Vector1,Vector2,a,b)
%This function performs numerical integration via the Simpsons's method.
%Input:     Vector1: vector of the "x" values to be integrated over
%           Vector2: vector of the "y(x)" values 
%           a: index value of Vector1 to begin the integration
%           b: index value of Vector1 to end the integration
%
%Output:    Output: scalar, the integral of Vector2 over the interval
%                   Vector1 from a to b

X = Vector1(a:b);
Y = Vector2(a:b);

%Integration via Simpsons Rule
    sum1 = 0;
    sum2 = 0;
    for n = 2:2:(length(X) - 1)
        sum1 = sum1 + Y(n);
    end
    for n = 3:2:(length(X) - 2)
        sum2 = sum2 + Y(n);
    end
    h = (X(end) - X(1))/(length(X)-1);
    Output = (h/3)*(Y(1) + 4*sum1 + 2*sum2 + Y(end));
    
%End of Function