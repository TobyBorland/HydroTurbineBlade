function [geom,iner,cpmo] = polygeom(x,y)

// Output variables initialisation (not found in input variables)
geom=[];
iner=[];
cpmo=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//POLYGEOM Geometry of a planar polygon
// 
//   POLYGEOM( X, Y ) returns area, X centroid,
//   Y centroid and perimeter for the planar polygon
//   specified by vertices in vectors X and Y.
// 
//   [ GEOM, INER, CPMO ] = POLYGEOM( X, Y ) returns
//   area, centroid, perimeter and area moments of 
//   inertia for the polygon.
//   GEOM = [ area   X_cen  Y_cen  perimeter ]
//   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
//     u,v are centroidal axes parallel to x,y axes.
//   CPMO = [ I1     ang1   I2     ang2   J ]
//     I1,I2 are centroidal principal moments about axes
//         at angles ang1,ang2.
//     ang1 and ang2 are in radians.
//     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv

// H.J. Sommer III - v02.05.14 - tested under MATLAB v5.2
// this function was posted on the MATLAB File Exchange

// sample data
// x = [ 2.000  0.500  4.830  6.330 ]'';
// y = [ 4.000  6.598  9.098  6.500 ]'';
// 3x5 test rectangle with long axis at 30 degrees
// area=15, x_cen=3.415, y_cen=6.549, perimeter=16
// Ixx=659.561, Iyy=201.173, Ixy=344.117
// Iuu=16.249, Ivv=26.247, Iuv=8.660
// I1=11.249, ang1=30deg, I2=31.247, ang2=120deg, J=42.496
// 
// credit for this function goes to:
// H.J. Sommer III, Ph.D., Professor of Mechanical Engineering, 337 Leonhard Bldg
// The Pennsylvania State University, University Park, PA  16802
// (814)863-8997  FAX (814)865-9693  hjs1@psu.edu  www.me.psu.edu/sommer/

// begin function POLYGEOM

// check if inputs are same size
if ~isequal(size(x),size(y)) then
  error("X and Y must be the same size");
end;

// number of vertices
// !! L.45: Matlab function shiftdim not yet converted, original calling sequence used.
[x,ns] = shiftdim(x);
// !! L.46: Matlab function shiftdim not yet converted, original calling sequence used.
[y,ns] = shiftdim(y);
[n,c] = size(x);

// temporarily shift data to mean of vertices for improved accuracy
xm = mean(x,"m");
ym = mean(y,"m");
x = mtlb_s(x,xm*ones(n,1));
y = mtlb_s(y,ym*ones(n,1));

// delta x and delta y
dx = mtlb_s(x([2:n,1]),x);
dy = mtlb_s(y([2:n,1]),y);

// summations for CW boundary integrals
A = mtlb_sum(mtlb_s(y .*dx,x .*dy))/2;
Axc = mtlb_sum(mtlb_a(mtlb_a(mtlb_s(((6*x) .*y) .*dx,((3*x) .*x) .*dy),((3*y) .*dx) .*dx),(dx .*dx) .*dy))/12;
Ayc = mtlb_sum(mtlb_s(mtlb_s(mtlb_s(((3*y) .*y) .*dx,((6*x) .*y) .*dy),((3*x) .*dy) .*dy),(dx .*dy) .*dy))/12;

Ixx = mtlb_sum(mtlb_s(mtlb_s(mtlb_s(mtlb_s(mtlb_s((((2*y) .*y) .*y) .*dx,(((6*x) .*y) .*y) .*dy),(((6*x) .*y) .*dy) .*dy),(((2*x) .*dy) .*dy) .*dy),(((2*y) .*dx) .*dy) .*dy),((dx .*dy) .*dy) .*dy))/12;

Iyy = mtlb_sum(mtlb_a(mtlb_a(mtlb_a(mtlb_a(mtlb_s((((6*x) .*x) .*y) .*dx,(((2*x) .*x) .*x) .*dy),(((6*x) .*y) .*dx) .*dx),(((2*y) .*dx) .*dx) .*dx),(((2*x) .*dx) .*dx) .*dy),((dx .*dx) .*dx) .*dy))/12;

Ixy = mtlb_sum(mtlb_s(mtlb_a(mtlb_s(mtlb_a(mtlb_s((((6*x) .*y) .*y) .*dx,(((6*x) .*x) .*y) .*dy),(((3*y) .*y) .*dx) .*dx),(((3*x) .*x) .*dy) .*dy),(((2*y) .*dx) .*dx) .*dy),(((2*x) .*dx) .*dy) .*dy))/24;
P = mtlb_sum(sqrt(mtlb_a(dx .*dx,dy .*dy)));

// check for CCW versus CW boundary
if mtlb_logic(A,"<",0) then
  A = -A;
  Axc = -Axc;
  Ayc = -Ayc;
  Ixx = -Ixx;
  Iyy = -Iyy;
  Ixy = -Ixy;
end;

// centroidal moments
xc = Axc/A;
yc = Ayc/A;
Iuu = mtlb_s(Ixx,(A*yc)*yc);
Ivv = mtlb_s(Iyy,(A*xc)*xc);
Iuv = mtlb_s(Ixy,(A*xc)*yc);
J = mtlb_a(Iuu,Ivv);

// replace mean of vertices
x_cen = mtlb_a(xc,xm);
y_cen = mtlb_a(yc,ym);
Ixx = mtlb_a(Iuu,(A*y_cen)*y_cen);
Iyy = mtlb_a(Ivv,(A*x_cen)*x_cen);
Ixy = mtlb_a(Iuv,(A*x_cen)*y_cen);

// principal moments and orientation
I = [Iuu,-Iuv;
     -Iuv,Ivv];
[eig_vec,eig_val] = spec(I);
I1 = eig_val(1,1);
I2 = eig_val(2,2);
// ! L.101: real(eig_vec(2,1)) may be replaced by:
// !    --> eig_vec(2,1) if eig_vec(2,1) is real
// ! L.101: real(eig_vec(1,1)) may be replaced by:
// !    --> eig_vec(1,1) if eig_vec(1,1) is real
ang1 = atan(real(eig_vec(2,1)),real(eig_vec(1,1)));
// ! L.102: real(eig_vec(2,2)) may be replaced by:
// !    --> eig_vec(2,2) if eig_vec(2,2) is real
// ! L.102: real(eig_vec(1,2)) may be replaced by:
// !    --> eig_vec(1,2) if eig_vec(1,2) is real
ang2 = atan(real(eig_vec(2,2)),real(eig_vec(1,2)));

// return values
geom = [A,x_cen,y_cen,P];
iner = [Ixx,Iyy,Ixy,Iuu,Ivv,Iuv];
cpmo = [I1,ang1,I2,ang2,J];

// end of function POLYGEOM

endfunction
