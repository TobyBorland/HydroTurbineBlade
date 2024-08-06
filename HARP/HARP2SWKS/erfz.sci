
function [e] = erfz(z)

// Output variables initialisation (not found in input variables)
e=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

if isempty(z) then
  e = z;
elseif isreal(z,0) then
  e = erf(z);
else
  R = real(z);
  I = imag(z);
  %v0_3 = size(z);  e = ones(%v0_3(1),%v0_3(2)) .*. %nan;
  n = (abs(R)<%inf) + (2*(abs(I)<%inf));
  // 
  // n  ->>  e
  // 
  // 0      nan
  // 1      nan
  // 2      erf(R)
  // 3      erf(R) + parts(R,I)
  // 
  k = n >= uint8(2);
  e(1,k) = matrix(erf(R(k)),1,-1);
  n = n == uint8(3) & I ~= uint8(0);

  if mtlb_any(n(:)) then
    I = mtlb_e(I,n);
    R = mtlb_e(R,n);
    R = parts(R,abs(I));
    k = I < uint8(0);
    R = mtlb_i(R,k,conj(mtlb_e(R,k)));
    e = mtlb_i(e,n,mtlb_a(e(n),R));
  end;
end;


//Evaluate all partial functions E(z) to G(z).
// 
//Because of the comparably large overhead for comparison
//and sub-indexing, this MATLAB implementation calculates
//simply all sub-functions regardless of the significance
//for the result.
// 
// R      Real part
// I      Imaginary part
// e      Total result
// 
endfunction

function [e] = parts(R,I)

// Output variables initialisation (not found in input variables)
e=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

R2 = R .*R;
e2iRI = exp(0 - %i*2*R.*I);
E = mtlb_s(1,e2iRI) ./((2*%pi)*R);
E = mtlb_i(E,~R,0);
F = 0;
Hr = 0;
Hi = 0;
N = sqrt(1-4*log(%eps/2));
for n = 1:ceil(N)
  H = (n*n)/4;
  H = exp(-H) ./mtlb_a(H,R2);
  F = mtlb_a(F,H);
  H = H .*exp(-n*I);
  Hi = mtlb_a(Hi,(n/2)*H);
  Hr = mtlb_a(Hr,H);
end;
e = exp(-R2) .*mtlb_s(mtlb_a(E,(R .*F)/%pi),(e2iRI .*(R.*Hr + %i*Hi))/(2*%pi));
clear("E","F","H*");
R3 = mtlb_a(R2,log(2*%pi));
Gr = 0;
Gi = 0;
M = 2*I;
n = mtlb_max(1,floor(mtlb_s(M,N)));
%v0 = mtlb_s(mtlb_a(M,N),n);M = ceil(mtlb_max(%v0,firstnonsingleton(%v0)));
for m = mtlb_imp(0,M)
  n1 = n/2;
  n2 = n1 .*n1;
  G = exp(mtlb_s(mtlb_s(mtlb_s(n.*I,n2),R3),log(mtlb_a(n2,R2))));
  Gi = mtlb_s(Gi, n1.*G);
  Gr = mtlb_a(Gr, G);
  n = n+1;
end;
e = mtlb_s(e,e2iRI .*(R.*Gr + %i*Gi));
n = mtlb_logic(R,"==",mtlb_uint8(0));
e = mtlb_i(e,n,mtlb_a(mtlb_e(e, n),(0 + %i*mtlb_e(I, n)./%pi)));
endfunction
