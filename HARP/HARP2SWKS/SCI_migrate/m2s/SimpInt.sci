function [Output] = SimpInt(Vector1,Vector2,a,b)

// Output variables initialisation (not found in input variables)
Output=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This function performs numerical integration via the Simpsons''s method.
//Input:     Vector1: vector of the """"x"""" values to be integrated over
//           Vector2: vector of the """"y(x)"""" values 
//           a: index value of Vector1 to begin the integration
//           b: index value of Vector1 to end the integration
// 
//Output:    Output: scalar, the integral of Vector2 over the interval
//                   Vector1 from a to b

X = mtlb_e(Vector1,mtlb_imp(a,b));
Y = mtlb_e(Vector2,mtlb_imp(a,b));

//Integration via Simpsons Rule
sum1 = 0;
sum2 = 0;
for n = 2:2:max(size(X))-1
  sum1 = mtlb_a(sum1,mtlb_e(Y,n));
end;
for n = 3:2:max(size(X))-2
  sum2 = mtlb_a(sum2,mtlb_e(Y,n));
end;
h = mtlb_s(mtlb_e(X,$),mtlb_e(X,1))/(max(size(X))-1);
Output = (h/3)*mtlb_a(mtlb_a(mtlb_a(mtlb_e(Y,1),4*sum1),2*sum2),mtlb_e(Y,$));

//End of Function
endfunction
