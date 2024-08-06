function [Cav_Error] = Cavitation_Check(RElm,LocVel,CpMin,NumCol)

// Output variables initialisation (not found in input variables)
Cav_Error=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This function checks to see if cavitation occurs at ANY element along the
//blade while the blade is rotated to the 12-o-clock position


global("NumSeg","ShaftTilt","PreCone","Patm","Pv","rho","HubHt","WatDepth","CavSF","SpdSt","SpdEnd","SpdDel")

g = 9.80665;//m/s2 acceleration due to gravity
NumCases = round(mtlb_a(mtlb_s(SpdEnd,SpdSt)/SpdDel,1));

Cav_Error = zeros(1,NumCol);
//absolute pressure calculated at the element''s highest point 
//of rotation (e.g. closest to the free surface).
d2r = %pi/180;
h = mtlb_s(mtlb_s(WatDepth,HubHt),(RElm*cos(ShaftTilt*d2r))*cos(PreCone*d2r));
Pabs = mtlb_a(Patm,(rho*g) .*h);

Pabs_rep = zeros((NumCases*NumSeg)*NumCol,1);
for n = mtlb_imp(1,NumCases*NumCol)
  Pabs_rep = mtlb_i(Pabs_rep,mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),Pabs);
end;

//Using the classical definition for now.  I was using the newly derived
//definition of the cavitation number, but it''s slightly faster to use the
//the classical definition, else Sigma must be calculated in a FOR loop
//which is slower.  WT_Perf will soon include the updated cavitation number
//equations.  The difference between the two equations is on the order of
//only a couple of percent.

Sigma = (2 .*mtlb_s(Pabs_rep,Pv)) ./(rho .*(LocVel .^2));//the classical definition
Cav = mtlb_logic(mtlb_a(Sigma ./CavSF,CpMin),"<",0);
for n = mtlb_imp(1,NumCol)
  // ! L.33: abs(Cav(mtlb_imp(mtlb_a(((n-1)*NumSeg)*NumCases,1),(NumSeg*NumCases)*n))) may be replaced by:
  // !    --> Cav(mtlb_imp(mtlb_a(((n-1)*NumSeg)*NumCases,1),(NumSeg*NumCases)*n)) if Cav(mtlb_imp(mtlb_a(((n-1)*NumSeg)*NumCases,1),(NumSeg*NumCases)*n)) is Real.

  if mtlb_any(abs(Cav(mtlb_imp(mtlb_a(((n-1)*NumSeg)*NumCases,1),(NumSeg*NumCases)*n)))) then
    Cav_Error(1,n) = 1;
  end;
end;


//End of Function
endfunction
