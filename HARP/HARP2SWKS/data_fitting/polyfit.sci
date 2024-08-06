function [paj, ycalc, statpar] = polyfit(x, y, n, Mcomp, GraphChk, pChar)
//--------------------------------------------------------------------------//
// polyfit: polynomial regression                                           //
// Author: Javier I. Carrero, jicarrerom@unal.edu.co                        // 
// Last modification date: 2009-12-14                                       //
//--------------------------------------------------------------------------//
// In its default behavior obtains a scilab polynomial object that best fit //
// the y(x) data in the form                                                //
//                                                                          //
// y = a0 + (a1*x) + (a2*x^2) + ... + (an*x^n )                             //
//                                                                          //
// ai values are optimal in the sense of minimizing the sum of square       // 
// residuals. When invoked with the Mcomp option polyfit returns the same   // 
// ai coefficients in the Matlab's form, a vector [an .. a1 a0].            //
//--------------------------------------------------------------------------// 
// Input arguments                                                          //
// * x: independent variable values, a (m,1) vector                         //
// * y: dependent variable values, a (m,1) vector                           //
// * n: integer, order of the polynomial                                    //
// * Mcomp (optional), specifies if output paj takes the Matlab form, to    //
//   invoke this option simply write the argument as MComp = 'Y'            //
// * GraphChk (optional) if present produces a graphic comparing the        //
//   tabulated values of the input (as circles) against the calculated      //
//   values (as a line). To invoke this option simply write the argument as //
//   GraphChk = 'Y'                                                         //
// * pChar (optional) a character to be used in the output polynomial paj,  // 
//   by default PolLetter is 'x', but it can be changed, for example        //
//   invoking the function wiht pChar = 'z'                                 //
//--------------------------------------------------------------------------//
// Output argumens                                                          //  
// * paj in the default form is a scilab polynomial of n degree that        //
//   adjust the input data y(x). If the Mcomp option was used it becomes    //
//   a (n+1, 1) matrix containing the values of ai, paj = [an ... a2 a1 a0] //
// * ycalc is a column vector with m elements corresponding to the y values //
//   calculated with the m sets of x1, x2, ..., xn                          //
// * statpar is a vector with statistical parameters,                       //
//   statpar = [St Sr stdv r2] where                                        //
//   - St: sum of the square differences between yi and yavg, being yavg    //
//     the average of y                                                     // 
//   - Sr: sum of the m residuals, being each residual (yi-yi_calc)^2       //
//     being ycalc the value obtained from the set x1i, x2i, ..., xni       // 
//   - stdv: standard deviation, defined as ( St / m-1 ) ^ 0.5              //
//   - r2: correlation coefficient, defined as r2 = ( St - Sr ) / St        //
//   - Syx: standard error, defined as Syx = (Sr/(m-(n+1)))^1/2             // 
//--------------------------------------------------------------------------//
// Mathematical background:                                                 //
// Based on the minimization of residuals, see polyfit documentation (pdf)  //
// and Chapra and Canale's "Numerical Methods for Engineers,                //
// 5th ed., ch. 17 (McGraw-Hill, 2005)                                      //
//--------------------------------------------------------------------------//  

//--------------------------------------------------------------------------//  
// Argument checking: x, y must be column vectors with equal lengths        //
//--------------------------------------------------------------------------//
termino = %f

[numrows numcols] = size(x)

if ( ~( ( numrows == 1 ) | ( numcols == 1 ) ) ) then
  printf(' Inconsistent argument: x must be a vector. ')
  termino = %t
else
  if ( numrows == 1 ) then
    x = x'
  end,  
end

[numrows numcols] = size(y)

if ( ~( ( numrows == 1 ) | ( numcols == 1 ) ) ) then
  printf(' Inconsistent argument: x must be a vector. ')
  termino = %t
else
  if ( numrows == 1 ) then
    y = y'
  end,  
end

NumPoints = length(x)

if ( length(y) <> NumPoints ) then
  printf(' Inconsistent arguments: x, and y must have the same length. ')
  termino = %t  
end,

if ( ( n < 0 ) | ( ( NumPoints - 1 ) < n  ) ) then
  printf(' Inconsistent argument: n must be in the range (0, m], where m=length(x). ')
  termino = %t
end

for i=1:NumPoints-1
   for j = i+1:NumPoints
      if ( abs( x(i) - x(j) ) <= %eps ) then
        printf(' Fail: there are repeated x values. ')
        termino = %t
      end
   end
end

if ( termino ) then
  abort
end,

if exists('Mcomp') then
  M_output = %t
else
  M_output = %f
end,

if exists('GraphChk') then
  DoGraph = %t
else
  DoGraph = %f
end,

if ( ~exists('pChar') ) then
  pChar = 'x'
end,

//--------------------------------------------------------------------------//
// Calculation procedure: since use of powers of x up to n produces         //
// ill-conditioned matrices QR decomposition must be used to solve the      //
// linear system of equations. As the polynomial has powers up to n there   //
// are n+1 unknown variables: a0, a1, ..., an stored in lstcoef             //
// System is solved using the upper diagonal R starting from n+1 row down   //
// to row 1                                                                 //
//--------------------------------------------------------------------------//
Z = zeros(NumPoints, n+1);
Z(:,1) = 1

for i=1:n
  Z(:,i+1) = x .^ i
end,

A = Z' * Z
B = Z' * y

[Q, R] = qr(A)

Qtb = Q' * B

lstcoef = zeros(n+1,1)
lstcoef(n+1) = Qtb(n+1) / R(n+1,n+1)

for j=n:-1:1
  lstcoef(j) = ( Qtb(j) - ( R(j,j+1:n+1) * lstcoef(j+1:n+1) ) ) / R(j,j)
end,

//--------------------------------------------------------------------------//
// Main output                                                              //
//--------------------------------------------------------------------------//
poli  = poly(lstcoef, pChar, 'coeff')
ycalc = horner(poli, x)

if ( M_output ) then
  paj = flipdim(lstcoef',2)
else
  paj = poli
end,

//------------------------------------------------------------//  
// Statistics                                                 //
//------------------------------------------------------------// 
avg_y = sum( y ) / NumPoints
St    = sum( ( y - avg_y ) .^ 2 )
Sr    = sum( ( ycalc - y ) .^ 2 )
stdv  = sqrt( St / ( NumPoints - 1 ) )
r2    = ( St - Sr ) / St

if ( NumPoints > (n + 1) ) then
  Syx = sqrt( Sr / ( NumPoints - ( n + 1 ) ) )
else
  Syx = %inf
end,

statpar = [St Sr stdv r2 Syx]

//------------------------------------------------------------//  
// Graphic checking                                           //
//------------------------------------------------------------//  
if (DoGraph) then
  Ntest = 100
  xtest = linspace(min(x), max(x), Ntest)
  ytest = horner(poli, xtest)
  
  clf()
  plot(x, y, 'o', xtest, ytest, '-')
  ejes = get("current_axes")
  ejes.x_label.text = pChar
  ejes.y_label.text = 'y'
end,

endfunction


