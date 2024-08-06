function [Strain,Iner,ShellArea] = getStrains(ST,STmax,CHORD,Mtan,Mnorm,Scaled_AF_Coordinates)

// Output variables initialisation (not found in input variables)
Strain=[];
Iner=[];
ShellArea=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


global("NumAFPoints","STdel","SFstruct","Einput")

//For the inside shell, need to make the coordinates in CCW direction for
//the polygeom.m function to work properly
%v0 = Scaled_AF_Coordinates(:,1);xx = %v0($:-1:1,:);
%v0 = Scaled_AF_Coordinates(:,2);yy = %v0($:-1:1,:);

//separate into upper and lower coordinates
Xlow = mtlb_e(xx,mtlb_imp(1,NumAFPoints/2));
Xup = mtlb_e(xx,mtlb_imp(mtlb_s(NumAFPoints,1),-1,NumAFPoints/2));
Ylow = mtlb_e(yy,mtlb_imp(1,NumAFPoints/2));
Yup = mtlb_e(yy,mtlb_imp(mtlb_s(NumAFPoints,1),-1,NumAFPoints/2));

// !! L.16: Matlab function atand not yet converted, original calling sequence used.
AngUp = atand(-1 ./FiniteDiff(Xup,Yup));
// !! L.17: Matlab function atand not yet converted, original calling sequence used.
AngLow = atand(-1 ./FiniteDiff(Xlow,Ylow));
//Make sure surface normal vectors are pointing towards inside of airfoil
for n = mtlb_imp(1,NumAFPoints/2)
  if mtlb_logic(mtlb_e(AngUp,n),">",0) then
    AngUp = mtlb_i(AngUp,n,mtlb_s(mtlb_e(AngUp,n),180));
  end;
  if mtlb_logic(mtlb_e(AngLow,n),"<",0) then
    AngLow = mtlb_i(AngLow,n,mtlb_a(mtlb_e(AngLow,n),180));
  end;
end;

%v02 = %f;if mtlb_logic(ST,"<=",STmax) then %v02 = mtlb_logic(ST,">",mtlb_s(STmax,STdel));end;
if %v02 then
  Xshell = Scaled_AF_Coordinates(:,1);
  Yshell = Scaled_AF_Coordinates(:,2);
  //fprintf(1,''reached STmax\n'');
else
  //Calculate the inner x,y coordinates based on a shell thickness measured normal to the surface    
  // !! L.34: Matlab function cosd not yet converted, original calling sequence used.
  InnerXup = mtlb_a(Xup,ST .*cosd(AngUp));
  // !! L.35: Matlab function sind not yet converted, original calling sequence used.
  InnerYup = mtlb_a(Yup,ST .*sind(AngUp));
  // !! L.36: Matlab function cosd not yet converted, original calling sequence used.
  InnerXlow = mtlb_a(Xlow,ST .*cosd(AngLow));
  // !! L.37: Matlab function sind not yet converted, original calling sequence used.
  InnerYlow = mtlb_a(Ylow,ST .*sind(AngLow));

  //Find any intersections of upper and lower inside surfaces
  // !! L.40: Unknown function intersections not converted, original calling sequence used.
  [Xo,Yo,I,J] = intersections(InnerXup,InnerYup,InnerXlow,InnerYlow,1);

  if isempty(Xo) then
    //Blade is completely solid material
    innerX = [];
    innerY = [];
  else
  
    innerX = [mtlb_e(Xo,1);InnerXlow(mtlb_imp(ceil(mtlb_e(J,1)),floor(mtlb_e(J,$))));mtlb_e(Xo,$);InnerXup(mtlb_imp(floor(mtlb_e(I,$)),-1,ceil(mtlb_e(I,1))));mtlb_e(Xo,1)];
  
    innerY = [mtlb_e(Yo,1);InnerYlow(mtlb_imp(ceil(mtlb_e(J,1)),floor(mtlb_e(J,$))));mtlb_e(Yo,$);InnerYup(mtlb_imp(floor(mtlb_e(I,$)),-1,ceil(mtlb_e(I,1))));mtlb_e(Yo,1)];
  end;

  Xshell = [Scaled_AF_Coordinates(:,1);innerX];
  Yshell = [Scaled_AF_Coordinates(:,2);innerY];
end;

//calculate the new values for Area, Inertia, and Centroids for the shell    
// !! L.58: Unknown function polygeom not converted, original calling sequence used.
[shellGEOM,shellINER] = polygeom(Xshell,Yshell);
Area = mtlb_e(shellGEOM,1);
Xcen_sh = mtlb_e(shellGEOM,2);
Ycen_sh = mtlb_e(shellGEOM,3);
Iuu_sh = mtlb_e(shellINER,4);
Ivv_sh = mtlb_e(shellINER,5);

//NOTE: for stress and strain calculations, the upper and lower points
//only consider the normal moment component and the LE and TE only
//consider the tangential moment component. In reality, both moment
//components (normal AND tangential) will contribute to the
//stress/strain at these points.

//Calculate strains  
cLE = Xcen_sh;
cTE = mtlb_s(CHORD,Xcen_sh);
%v0 = Scaled_AF_Coordinates(mtlb_imp(1,NumAFPoints/2),2);cUpper = abs(mtlb_s(mtlb_max(%v0,firstnonsingleton(%v0)),Ycen_sh));
%v0 = Scaled_AF_Coordinates(mtlb_imp(NumAFPoints/2,mtlb_s(NumAFPoints,1)),2);cLower = abs(mtlb_s(mtlb_min(%v0,firstnonsingleton(%v0)),Ycen_sh));
StrainLE = (((SFstruct .*Mtan) .*1000) .*cLE) ./((1000000000 .*Einput) .*Ivv_sh);//(m/m)
StrainTE = (-((SFstruct .*Mtan) .*1000) .*cTE) ./((1000000000 .*Einput) .*Ivv_sh);//(m/m)
StrainUpper = (((SFstruct .*Mnorm) .*1000) .*cUpper) ./((1000000000 .*Einput) .*Iuu_sh);//(m/m)
StrainLower = (-((SFstruct .*Mnorm) .*1000) .*cLower) ./((1000000000 .*Einput) .*Iuu_sh);//(m/m)

//function outputs
Strain.LE = StrainLE;
Strain.TE = StrainTE;
Strain.Upper = StrainUpper;
Strain.Lower = StrainLower;
Iner.uu = Iuu_sh;
Iner.vv = Ivv_sh;
ShellArea = Area;

endfunction
