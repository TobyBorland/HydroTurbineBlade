
function [fs,fc,ff,fg] = fresnel(x,xc)

// Output variables initialisation (not found in input variables)
fs=[];
fc=[];
ff=[];
fg=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//FRESNEL(x) Calculate Fresnel sin and cos integrals
//%
//         /x
// fs(x) = |   sin(pi/2*x^2) dx 
//         /0
// 
// 
//         /x
// fc(x) = |   cos(pi/2*x^2) dx 
//         /0
// 
//  d fs
//  ---- =  sin(pi/2*x^2)
//  d x
// 
// 
//  d fc
//  ---- =  cos(pi/2*x^2)
//  d x
// 
//  for |x| < 1.8 uses a powers series about x=0         
//                accuracy is better than 5e-16
//  for |x| > 1.8 uses a minimax rational approximation
//                based on the auxilliary functions f and g
//                accuracy is better than 1e-9
// 
// The accuracy can be checked by numerically calculating
// dfs/dx and dfc/dx. For all x this is better than 2e-8 
// 
// To test the accuracy call with no arguments
// Approximations and 20 digit accurate results were generated with maple 6
// 
// Definitions for the Fresnel integrals and auxialliary functions
// f and g are taken from 
//      Handbook of Mathematical Functions
//      Abramowitz and Stegun (9th printing 1970)
// 
// uses power series for x<xc=1.8                   accurate to machine precision
// uses minimax rational approximation for x>xc     accurate to about 1e-9
// 
// 
// Released under the GPL
// jnm11@amtp.cam.ac.uk
// J. N. McElwaine
// Department of Applied Mathematcis and Theoretical Physics
// University of Cambridge 
// Silver Street
// UK
// Version 1.0 October 2001
// 
// The accuracy should be improved for larger values of x


if %nargin==0 then
  // !! L.56: Unknown function test_fresnel not converted, original calling sequence used.
  test_fresnel;
  return;
end;


xc = 1.8;// Cut point for switching between power series and rational approximation

sgnx = sign(x);
x = abs(x);

f1 = mtlb_find(mtlb_logic(x,"<=",xc));
f2 = mtlb_find(mtlb_logic(x,">",xc));

s = (%pi/2)*(mtlb_e(x,f1) .^2);
ss = s .^2;
is = sqrt(2/%pi) ./mtlb_e(x,f2);

%v0 = size(x);fs = ones(%v0(1),%v0(2)) .*. %nan;
fc = fs;
ff = fs;
fg = fs;

// Approximations for x < 1.6 accurate to eps

fs(1,f1) = matrix((mtlb_e(x,f1) .*s) .*mtlb_a(1/3,mtlb_a(-0.0238095238095,mtlb_a(0.0007575757576,mtlb_a(-0.0000132275132,mtlb_a(0.0000001450385,mtlb_a(-0.0000000010892,mtlb_a(0.0000000000059,mtlb_a(-2.466827010D-14,mtlb_a(8.032735012D-17,mtlb_a(-2.107855191D-19,mtlb_a(4.551846759D-22,mtlb_a(-8.230149299D-25,mtlb_a(1.264107899D-27,mtlb_a(-1.669761793D-30,1.916942862D-33 .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss),1,-1);

fc(1,f1) = matrix(mtlb_e(x,f1) .*mtlb_a(1,mtlb_a(-0.1,mtlb_a(0.0046296296296,mtlb_a(-0.0001068376068,mtlb_a(0.0000014589169,mtlb_a(-0.0000000131225,mtlb_a(0.0000000000835,mtlb_a(-0.0000000000004,mtlb_a(1.448326464D-15,mtlb_a(-4.221407289D-18,mtlb_a(1.002516493D-20,mtlb_a(-1.977064754D-23,mtlb_a(3.289260349D-26,mtlb_a(-4.678483516D-29,5.754191644D-32 .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss) .*ss),1,-1);

ff(1,f2) = matrix(((sqrt(2/%pi)*is) .*mtlb_a(0.3663490102576,mtlb_a(0.8469566622219,mtlb_a(-3.730127348735,mtlb_a(10.520237918558,mtlb_a(-10.472617402498,4.4169634834911 .*is) .*is) .*is) .*is) .*is)) ./mtlb_a(0.7326980266120,mtlb_a(1.6939102288853,mtlb_a(-7.4599994789665,mtlb_a(21.032436583849,mtlb_a(-20.269535575486,8.9995024877629 .*is) .*is) .*is) .*is) .*is),1,-1);

fg(1,f2) = matrix(((sqrt(2/%pi)*(is .^3)) .*mtlb_a(0.1093569232008,mtlb_a(0.7202553305554,mtlb_a(-2.5630322339412,mtlb_a(7.2404268469721,mtlb_a(-7.0473933910698,2.9315207450904 .*is) .*is) .*is) .*is) .*is)) ./mtlb_a(0.4374277218534,mtlb_a(2.8810049959848,mtlb_a(-10.250672312139,mtlb_a(28.912791022418,mtlb_a(-25.740131167525,15.14513436371 .*is) .*is) .*is) .*is) .*is),1,-1);

s = (%pi/2)*(mtlb_e(x,f2) .^2);
cx = cos(s);
sx = sin(s);
fs = mtlb_i(fs,f2,1/2-ff(f2) .*cx-fg(f2) .*sx);
fc = mtlb_i(fc,f2,1/2+ff(f2) .*sx-fg(f2) .*cx);

fs = fs .*sgnx;
fc = fc .*sgnx;


endfunction

function [] = test_fresnel()

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


x1 = 0:0.1:1.6;
s1 = [0;'0.0005235895476;'0.0041876091617;'0.0141169980066;'0.0333594326606;'0.0647324328600;'0.1105402073594;'0.1721364578635;'0.2493413930539;'0.3397763443931;'0.4382591473904;'0.5364979110968;'0.6234009185462;'0.6863332855347;'0.7135250773634;'0.6975049600821;'0.6388876835094];

c1 = [0;'0.0999975326271;'0.1999210575945;'0.2994009760520;'0.3974807591724;'0.4923442258714;'0.5810954469917;'0.6596523519045;'0.7228441718964;'0.7648230212733;'0.7798934003768;'0.7638066660620;'0.7154377229231;'0.6385504547270;'0.5430957835463;'0.4452611760398;'0.3654616834405];

x2 = [1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.5,4.5,5,6,7,8,9,10,15,20,50];

// map(Fresnelf[1.61.71.81.92.2.22.42.62.83.03.54.55.6.7.8.9.10.15.20.50.]);
f1 = [0.1921938959194;'0.1820080012233;'0.1727448359130;'0.1643008528219;'0.1565843216363;'0.1430201832028;'0.1315175579491;'0.1216656017292;'0.1131482205164;'0.1057207892977;'0.0907655583153;'0.0706835395885;'0.0636311887040;'0.0530392387631;'0.0454670925470;'0.0397857856070;'0.0353661274681;'0.0318300214151;'0.0212205316744;'0.0159154640740;'0.0063661974141];

// map(Fresnelg[1.61.71.81.92.2.22.42.62.83.03.54.55.6.7.8.9.10.15.20.50.]);

g1 = [0.0212568488870;'0.0181740929177;'0.0156278255443;'0.0135129446179;'0.0117465939247;'0.0090108056826;'0.0070410079902;'0.0055940626785;'0.0045112596247;'0.0036870010326;'0.0023401756317;'0.0011078332747;'0.0008086180829;'0.0004685321445;'0.0002952105466;'0.0001978196228;'0.0001389543703;'0.0001013057945;'0.0000300201903;'0.0000126650277;'0.0000008105693];

//map(FresnelS,[1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3.0,3.5,4.5,5,6,7,8,9,10,15,20,50]);

s2 = [0.6388876835094;'0.5491959403216;'0.4509387692676;'0.3733473178170;'0.3434156783637;'0.4557046121247;'0.6196899649457;'0.5499893231527;'0.3915284435432;'0.4963129989674;'0.4152480119724;'0.4342729750487;'0.4991913819171;'0.4469607612369;'0.4997047894534;'0.4602142143930;'0.4998610456297;'0.4681699785849;'0.4999699798097;'0.4840845359260;'0.4936338025859];

//map(FresnelC,[1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3.0,3.5,4.5,5,6,7,8,9,10,15,20,50]);

c2 = [0.3654616834405;'0.3238268760039;'0.3336329272216;'0.3944705348915;'0.4882534060753;'0.6362860449033;'0.5549614058564;'0.3889374961920;'0.4674916516989;'0.6057207892977;'0.532572435028;'0.5260259150535;'0.5636311887040;'0.4995314678555;'0.5454670925470;'0.4998021803772;'0.5353661274681;'0.4998986942055;'0.5212205316744;'0.4999873349723;'0.4999991894307];

[fs,fc,ff,fg] = fresnel(x1);

subplot(3,1,1)
plot(x1,fs-s1,x1,fc-c1);
title("Fresnel integrals errors for small x");
subplot(3,1,2)
[fs,fc,ff,fg] = fresnel(x2);
plot(x2,fs-s2,x2,fc-c2);
title("Fresnel integrals errors for large x");

subplot(3,1,3)
x = linspace(0,10,1000);
dx = %eps+x*sqrt(%eps);
[sp,cp] = fresnel(x+dx);
[sm,cm] = fresnel(x-dx);
plot(x,cos((%pi/2)*(x .^2))-(cp-cm) ./(2*dx));
plot(x,sin((%pi/2)*(x .^2))-(sp-sm) ./(2*dx));
title("Differential of Fresnel integrals errors for all x");


return;
endfunction
