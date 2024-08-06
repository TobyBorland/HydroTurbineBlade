function [y,N,algo] = erf_(x,N,algo)

// Output variables initialisation (not found in input variables)
y=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// error function for a real or complex array
// 
// Note: The Fresnel integrals C(x) and S(x) can be computed for real x as
//   C(x)+1i*S(x) == ((1+1i)/2)*erf_((sqrt(pi)*(1-1i)/2)*x)
// (See demo example at the end of this file.)
// 
// syntax:
//   y = erf_(x);
//   [y,N,algo] = erf_(x,N,algo);
// 
// inputs:
// 
//   x: numeric array, can be complex
// 
// optional inputs:
// 
//   N: positive integer - scalar or size-matched to x, or [] (optional; unspecified or
//   [] defaults to 1000)
//   maximum number of terms to use in series or continued-fraction expansion.  N can be
//   specified independently for each x element
// 
//   algo: integers 1, 2, or 3 - scalar or size-matched to x, or [] (optional; default
//   is [])
//   Set algo = 1 to force use of Taylor series expansion (good for small arguments).
//   Set algo = 2 to force use of erf continued-fraction expansion (good for arguments
//   with large real part).
//   Set algo = 3 to force use of Dawson''s-function continued-fraction expansion (good
//   for arguments with large imaginary part).
//   Leave algo unspecified or [] to make the choice automatically.
//   algo can be specified independently for each x element.
// 
// outputs:
// 
//   y: numeric array, size-matched to x
//   y = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt.
// 
// optional outputs:
// 
//   N: integer array, size-matched to x
//   number of terms used in series or continued-fraction expansion for each x element
// 
//   algo: integer array of values 1, 2, or 3; size-matched to x
//   This indicates which algorithm was used for each x. 1: Taylor series, 2: erf
//   continued fraction, 3: Dawson''s-function continued fraction
// 
// 
// Author: Kenneth C. Johnson, kjinnovation@earthlink.net, kjinnovation.com
// Version: November 1, 2011
// 
// 
// BSD Copyright notice:
// 
// Copyright 2011 by Kenneth C. Johnson (kjinnovation.com)
// 
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that this entire copyright notice is duplicated in all such copies.
// 
// This software is provided """"as is"""" and without any express or implied warranties,
// including, without limitation, the implied warranties of merchantibility and fitness
// for any particular purpose.
// 

if %nargin<3 then
  algo = [];
  if %nargin<2 then
    N = [];
  end;
end;

size_ = mtlb_size(x);

if isempty(N) then
  N = 1000;
end;
if sum(length(N))==1 then
  N = repmat(N,size_);
elseif ~isequal(mtlb_size(N),size_) then
  error("erf_:validation, Size mismatch between args N and x.")
end;

if isempty(algo) then
  // Algorithm selection for erf_(x):
  // Apply algo 1 (Taylor series) for |real(x)|<=1, |imag(x)|<=6.
  // Apply algo 3 (Dawson''s-function continued fraction) for |real(x)|<=1, |imag(x)|>6.
  // Apply algo 2 (erf continued fraction) for |real(x)|>1.
  // See test code at the bottom of this file for checking algorithm continuity across
  // algo boundaries.
  algo = ones(mtlb_size(x));
  algo(abs(imag(x))>6) = 3;
  algo(abs(real(x))>1) = 2;
else
  if sum(length(algo))==1 then
    algo = repmat(algo,size_);
  elseif ~isequal(mtlb_size(algo),size_) then
    error("erf_:validation, Size mismatch between args algo and x.")
  end;
  if ~and(algo==1 | algo==2 | algo==3) then
    error("erf_:validation, Invalid arg: algo values must be 1, 2, or 3 .")
  end;
end;

x = x(:);
N = N(:);
algo = algo(:);
y = zeros(mtlb_size(x));
sgn = ones(mtlb_size(x));
sgn(real(x)<0) = -1;
x = sgn .*x;// Only apply algorithm with abs(angle(x))<=pi/2; use erf(x)==-erf(-x).

warn = %f;// for possible non-convergence
i = find(algo==1);
if ~isempty(i) then
  // Apply series expansion to x(i) (see http://mathworld.wolfram.com/Erf.html, Eq 6).
  n = 0;  // summation index
  Ni = N(i);
  term = x(i);
  xx = term .*term;
  term = (2/sqrt(%pi))*term;  // n-th summation term
  psum = term;  // partial sum up to n-th term
  prec = abs(psum)*%eps;  // estimated precision of partial sum
  while %t
    n = n+1;
    term = ((-(2*n-1)/(n*(2*n+1)))*term) .*xx;
    psum = psum+term;
    test = (abs(term)<=prec);
    if and(test) then
      y(i(test)) = psum(test);
      N(i(test)) = n;
      i(test) = [];
      if isempty(i)
        break
      end
      Ni(test) = [];
      term(test) = [];
      psum(test) = [];
      prec(test) = [];
      xx(test) = [];
    end;
    test = (n==Ni);
    if and(test)
      warn = true;
      y(i(test)) = psum(test);
      N(i(test)) = n;
      i(test) = [];
      if isempty(i)
        break
      end
      Ni(test) = [];
      term(test) = [];
      psum(test) = [];
      prec(test) = [];
      xx(test) = [];
    end
    prec = max(prec,abs(psum)*%eps);
  end
end


i = find(algo==2);
if ~isempty(i) then
  // Apply erf continued fraction to x(i).  Use Eq 6.9.4 (second form) in Numerical
  // Recipes in C++, 2nd Ed. (2002).  Apply Lentz''s method (Sec. 5.2) to top-level
  // denominator in Eq 6.9.4.
    n = 0;
    Ni = N(i);
    b = x(i);
    b = 2*b.*b+1;
    f = b;
    C = f;
    D = zeros(mtlb_size(i));
    while %t
        n = n+1;
        a = -(2*n-1)*(2*n);
        b = b+4;
        D = b + a*D;
        D(D==0) = %eps^2;
        D = 1./D;
        C = b+a./C;
        C(C==0) = %eps^2;
        Delta = C.*D;
        f = f.*Delta;
        test = (abs(Delta-1)<=%eps);
        if and(test)
            y(i(test)) = f(test);
            N(i(test)) = n;
            i(test) = [];
            if isempty(i)
                break
            end
            Ni(test) = [];
            b(test) = [];
            f(test) = [];
            C(test) = [];
            D(test) = [];
        end
        test = (n==Ni);
        if and(test)
            warn = %t;
            y(i(test)) = f(test);
            N(i(test)) = n;
            i(test) = [];
            if isempty(i)
                break
            end
            Ni(test) = [];
            b(test) = [];
            f(test) = [];
            C(test) = [];
            D(test) = [];
        end
    end
    i = find(algo==2);
    y(i) = 1-exp(-x(i).^2).*((2/sqrt(%pi))*x(i)./y(i));
end

i = find(algo==3);
if ~isempty(i) then
  // Apply Dawson''s-function (F) continued fraction to x(i).  See
  // http://mathworld.wolfram.com/DawsonsIntegral.html, Eq 14.
  // erf(x) = i*erfi(-i*x) = (2i/sqrt(pi))*exp(-x^2)*F(-i*x)
  // Apply Lentz''s method to top-level denominator in Eq 14.
  n = 0;
  Ni = N(i);
  xx = x(i);
  xx = -xx.*xx; // (-1i*x)^2
  b = 1+2*xx;
  f = b;
  C = f;
  D = zeros(mtlb_size(i));
  while true
    n = n+1;
    a = -4*n*xx;
    b = 2+b;
    D = b+a.*D;
    D(D==0) = %eps^2;
    D = 1./D;
    C = b+a./C;
    C(C==0) = %eps^2;
    Delta = C.*D;
    f = f.*Delta;
    test = (abs(Delta-1)<=%eps);
    if any(test)
      y(i(test)) = f(test);
      N(i(test)) = n;
      i(test) = [];
      if isempty(i)
        break
      end
      Ni(test) = [];
      xx(test) = [];
      b(test) = [];
      f(test) = [];
      C(test) = [];
      D(test) = [];
    end
    test = (n==Ni);
    if and(test)
      warn = %t;
      y(i(test)) = f(test);
      N(i(test)) = n;
      i(test) = [];
      if isempty(i)
        break
      end
      Ni(test) = [];
      xx(test) = [];
      b(test) = [];
      f(test) = [];
      C(test) = [];
      D(test) = [];
    end
  end
  i = find(algo==3);
  // Current state: y(i)==(-1i*x(i))./F(-1i*x(i)).  Convert F(-1i*x) to erfi(x).
  y(i) = (2/sqrt(pi))*exp(-x(i).^2).*x(i)./y(i);
end

  

y = sgn .*y;
y = matrix(y,size_);
N = matrix(mtlb_double(N),size_);
algo = matrix(mtlb_double(algo),size_);

if warn then
  warning("erf_:convergence Possible non-convergence in erf_");
end;

//disp('x');disp(x);
//disp('y');disp(y);
//pause
return;

//--------------------------------------------------------------------------------------
// Demo example: Cornu Spiral (http://en.wikipedia.org/wiki/Cornu_spiral)
//FresnelCS = @(x) ((1+%i)/2)*erf_((sqrt(pi)*(1-%i)/2)*x); % Fresnel integrals
//CS = FresnelCS(-10:.001:10);
//figure, plot(real(CS),imag(CS)), axis equal


//--------------------------------------------------------------------------------------
// Test code for checking algorithm discontinuity across algo boundaries:

// boundary between algo 1 and algo 2:
xdel(winsid())
x = (0:0.02:3).';y = [0,6];%v0 = mtlb_size(y);%v1 = mtlb_size(x);z = mtlb_a(ones(%v0(1),%v0(2)) .*. x,(1*%i)*ones(%v1(1),%v1(2)) .*. y);
// !! L.303: Unknown function erf_ not converted, original calling sequence used.
[e,N] = erf_(z,[],1);// !! L.303: Unknown function erf_ not converted, original calling sequence used.
[e_,N_] = erf_(z,[],2);d = log10(abs((2*mtlb_s(e,e_)) ./mtlb_a(e,e_)));
// !! L.304: Matlab function figure not yet converted.
// !! L.304: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(x,N),// !! L.304: Matlab function ylim not yet converted.
plot([1,1],mtlb(ylim),":k")
xlabel("x"),ylabel("N, algo 1"),// !! L.305: Matlab function legend not yet converted, original calling sequence used.
legend("y = 0","y = 6")
// !! L.306: Matlab function figure not yet converted.
// !! L.306: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(x,N_),// !! L.306: Matlab function ylim not yet converted.
plot([1,1],mtlb(ylim),":k")
xlabel("x"),ylabel("N, algo 2"),// !! L.307: Matlab function legend not yet converted, original calling sequence used.
legend("y = 0","y = 6")
// !! L.308: Matlab function figure not yet converted.
// !! L.308: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(x,d),// !! L.308: Matlab function ylim not yet converted.
plot([1,1],mtlb(ylim),":k")
xlabel("x"),ylabel("log10(algo 1-2 diff)"),// !! L.309: Matlab function legend not yet converted, original calling sequence used.
legend("y = 0","y = 6")
title("erf\_(x+1i*y), algo 1 vs 2 (boundary at |x|==1, |y|<=6)")

// boundary between algo 3 and algo 2:
xdel(winsid())
x = (0:0.02:3).';y = [6,10];%v0 = mtlb_size(y);%v1 = mtlb_size(x);z = mtlb_a(ones(%v0(1),%v0(2)) .*. x,(1*%i)*ones(%v1(1),%v1(2)) .*. y);
// !! L.315: Unknown function erf_ not converted, original calling sequence used.
[e,N] = erf_(z,[],3);// !! L.315: Unknown function erf_ not converted, original calling sequence used.
[e_,N_] = erf_(z,[],2);d = log10(abs((2*mtlb_s(e,e_)) ./mtlb_a(e,e_)));
// !! L.316: Matlab function figure not yet converted.
// !! L.316: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(x,N),// !! L.316: Matlab function ylim not yet converted.
plot([1,1],mtlb(ylim),":k")
xlabel("x"),ylabel("N, algo 3"),// !! L.317: Matlab function legend not yet converted, original calling sequence used.
legend("y = 6","y = 10")
// !! L.318: Matlab function figure not yet converted.
// !! L.318: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(x,N_),// !! L.318: Matlab function ylim not yet converted.
plot([1,1],mtlb(ylim),":k")
xlabel("x"),ylabel("N, algo 2"),// !! L.319: Matlab function legend not yet converted, original calling sequence used.
legend("y = 6","y = 10")
// !! L.320: Matlab function figure not yet converted.
// !! L.320: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(x,d),// !! L.320: Matlab function ylim not yet converted.
plot([1,1],mtlb(ylim),":k")
xlabel("x"),ylabel("log10(algo 3-2 diff)"),// !! L.321: Matlab function legend not yet converted, original calling sequence used.
legend("y = 6","y = 10")
title("erf\_(x+1i*y), algo 3 vs 2 (boundary at |x|==1, |y|>6)")

// boundary between algo 1 and algo 3:
xdel(winsid())
y = (0:0.02:10).';x = [0,1];%v0 = mtlb_size(y);%v1 = mtlb_size(x);z = mtlb_a(ones(%v0(1),%v0(2)) .*. x,(1*%i)*ones(%v1(1),%v1(2)) .*. y);
// !! L.327: Unknown function erf_ not converted, original calling sequence used.
[e,N] = erf_(z,[],1);// !! L.327: Unknown function erf_ not converted, original calling sequence used.
[e_,N_] = erf_(z,[],3);d = log10(abs((2*mtlb_s(e,e_)) ./mtlb_a(e,e_)));
// !! L.328: Matlab function figure not yet converted.
// !! L.328: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(y,N),// !! L.328: Matlab function ylim not yet converted.
plot([6,6],mtlb(ylim),":k")
xlabel("y"),ylabel("N, algo 1"),// !! L.329: Matlab function legend not yet converted, original calling sequence used.
legend("x = 0","x = 1")
// !! L.330: Matlab function figure not yet converted.
// !! L.330: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(y,N_),// !! L.330: Matlab function ylim not yet converted.
plot([6,6],mtlb(ylim),":k")
xlabel("y"),ylabel("N, algo 3"),// !! L.331: Matlab function legend not yet converted, original calling sequence used.
legend("x = 0","x = 1")
// !! L.332: Matlab function figure not yet converted.
// !! L.332: Matlab function figure not yet converted, original calling sequence used.
figure,set(gca(),"auto_clear","off"),plot(y,d),// !! L.332: Matlab function ylim not yet converted.
plot([6,6],mtlb(ylim),":k")
xlabel("y"),ylabel("log10(algo 1-3 diff)"),// !! L.333: Matlab function legend not yet converted, original calling sequence used.
legend("x = 0","x = 1")
title("erf\_(x+1i*y), algo 1 vs 3 (boundary at |x|<=1, |y|==6)")

endfunction
