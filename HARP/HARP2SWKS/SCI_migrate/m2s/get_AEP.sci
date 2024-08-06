function [AEP] = get_AEP(PvsV,pUvars)

// Output variables initialisation (not found in input variables)
AEP=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This funtion calculate the annual energy production
//Inputs:    PvsV: a two column vector, first column is the flow speeds, second column is the power values
//           pUvars: vector which contains the appropriate variables needed for the selected probability distribution, pUvars  is defined in Main.m
// 
//Outputs:   AEP: scalar, the AEP value

global("user_pU_interp","ProbDist")

U = PvsV(:,1);
P = PvsV(:,2);

if mtlb_logic(ProbDist,"==",1) then
  //Rayleigh Distribution
  U_mean = pUvars;
  pU = ((%pi ./2) .*(U ./(U_mean .^2))) .*exp(-(%pi ./4) .*((U ./U_mean) .^2));

elseif mtlb_logic(ProbDist,"==",2) then
  //Weibull Distribution
  k = mtlb_e(pUvars,1);
  c = mtlb_e(pUvars,2);
  pU = ((k ./c) .*((U ./c) .^mtlb_s(k,1))) .*exp(-(U ./c) .^k);

elseif mtlb_logic(ProbDist,"==",3) then
  //User Defined Distribution, this was already read in and interpolated in the HARP_Opt.m file
  //Make sure the Power curve and flow distribution vectors are the same length
  pU = mtlb_e(user_pU_interp,1:max(size(U)));
elseif mtlb_logic(ProbDist,"==",0) then
  AEP = 0;
  return;
end;

//Average Power
// !! L.34: Unknown function SimpInt not converted, original calling sequence used.
Pavg = SimpInt(U,P .*pU,1,max(size(U)));

//Annual Energy Production
Hours = 8760;
AEP = Pavg*Hours;

//End of Function
endfunction
