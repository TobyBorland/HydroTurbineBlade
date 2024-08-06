function [ShapeError,RElm,TWIST,CHORD,PERCENT_THICKNESS,DIMENSIONAL_THICKNESS,R_CHORD_CP,CHORD_CP,R_TWIST_CP,TWIST_CP,THICK_CP] = Define_Blade_Shape(x)

// Output variables initialisation (not found in input variables)
ShapeError=[];
RElm=[];
TWIST=[];
CHORD=[];
PERCENT_THICKNESS=[];
DIMENSIONAL_THICKNESS=[];
R_CHORD_CP=[];
CHORD_CP=[];
R_TWIST_CP=[];
TWIST_CP=[];
THICK_CP=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);



global("NumSeg","RotorRad","HubRad","Thickness_values","ThickMethod","CircleRoot","minRootChord","maxRootChord","RootTranSt","RootTranEnd","RadialSpacing")

//=========================================================================%

//This section defines the blade geometery
ShapeError = 0;//initial value

//Define the Twist and Chord control points
Twist_CP(1:length(mtlb_t(mtlb_e(x,1:5))),1) = mtlb_t(mtlb_e(x,1:5));
Chord_CP(1:length(mtlb_t(mtlb_e(x,6:10))),1) = mtlb_t(mtlb_e(x,6:10));

//Define the radial locations of the control points: using equal or cosine spacing
usecosine = 1;
if usecosine==0 then //Using equal distance spacing for the design variables
 if mtlb_logic(CircleRoot,"==",0) then
   rad_CP = mtlb_t(mtlb_linspace(HubRad,RotorRad,5));
 elseif mtlb_logic(CircleRoot,"==",1) then
   rad_CP = mtlb_t(mtlb_linspace(RootTranEnd*RotorRad,RotorRad,5));
 end;
elseif usecosine==1 then //Using half cosine spacing for the design variables
 if mtlb_logic(CircleRoot,"==",0) then
   // !! L.25: Unknown function hcosspace not converted, original calling sequence used.
   rad_CP = hcosspace(HubRad,RotorRad,5,0);
 elseif mtlb_logic(CircleRoot,"==",1) then
   // !! L.27: Unknown function hcosspace not converted, original calling sequence used.
   rad_CP = hcosspace(RootTranEnd*RotorRad,RotorRad,5,0);
 end;
end;

//Define the radia locations of the WT_Perf blade elements: radial or cosine spacing
if mtlb_logic(RadialSpacing,"==",0) then //equal spacing of blade elements  
 //Define the radius of the center of the blade elements
 delElm = mtlb_s(RotorRad,HubRad)/NumSeg;
 RElm = mtlb_imp(mtlb_a(HubRad,delElm/2),delElm,mtlb_a(mtlb_a(HubRad,delElm/2),mtlb_s(NumSeg,1)*delElm))'; //Values at the center of each element (as defined in WT_Perf)

elseif mtlb_logic(RadialSpacing,"==",1) then //cosine spacing of blade elements
 // !! L.38: Unknown function hcosspace not converted, original calling sequence used.
 RElm_cs_EndPts = hcosspace(HubRad,RotorRad,mtlb_a(NumSeg,1),3); //cosine spacing values at the boundaries of each element
 //calculate the centers of the elements for cosine spacing
 RElm = zeros(NumSeg,1);
 for n = mtlb_imp(1,NumSeg)
   RElm(n,1) = mtlb_a(mtlb_e(RElm_cs_EndPts,n+1),mtlb_e(RElm_cs_EndPts,n))/2;
 end;
end;

if mtlb_logic(CircleRoot,"==",0) then //Non-circular root
 //From the control points, interpolate using a Bezier curve
 B_Chord = Bezier([rad_CP,Chord_CP],NumSeg);
 B_Twist = Bezier([rad_CP,Twist_CP],NumSeg);
 //Now interpolate to spacing used in RElm
 CHORD = interp1(B_Chord(:,1),B_Chord(:,2),RElm);
 TWIST = interp1(B_Twist(:,1),B_Twist(:,2),RElm);
 R_CHORD_CP = rad_CP;
 R_TWIST_CP = rad_CP;
 CHORD_CP = Chord_CP;
 TWIST_CP = Twist_CP;

 //Now define the percent thickness distribution
 %v1_2 = %f; if mtlb_logic(ThickMethod,"==",2) then %v1_2 = max(size(Thickness_values))>1;end;
 if %v03 then //Piecewise Constant
  THICK_CP(1:length(RotorRad*mtlb_e(x,11:10+max(size(Thickness_values))-1)),1) = RotorRad*mtlb_e(x,11:10+max(size(Thickness_values))-1);
  THICK_CP = THICK_CP;
  //create a step function
  Thickness(1:max(size(RElm)),1) = mtlb_e(Thickness_values,$);
  a = 1; //intial index value
  for n = 1:max(size(THICK_CP))
    for k = a:max(size(RElm))
      if mtlb_logic(RElm(k),"<",THICK_CP(n)) then
        Thickness(k,1) = mtlb_e(Thickness_values,n);
      elseif mtlb_logic(RElm(k),">",THICK_CP(n)) then
        a = k;
        break
      end;
    end;
  end;
  Thickness = mtlb_i(Thickness,$+1:max(size(RElm)),mtlb_e(Thickness_values,$));
  PERCENT_THICKNESS = Thickness;
  DIMENSIONAL_THICKNESS = (CHORD .*PERCENT_THICKNESS) ./100;
 elseif %v1_2 then //Piecewise Linear
  THICK_CP(1:length(RotorRad*mtlb_e(x,11:10+max(size(Thickness_values)))),1) = RotorRad*mtlb_e(x,11:10+max(size(Thickness_values)));
  THICK_CP = THICK_CP;
  //Make sure that none of the control points are equal, so the interpolation does not crash
  for n = 2:max(size(THICK_CP))
    if mtlb_logic(THICK_CP(n),"<=",THICK_CP(n-1)) then
      THICK_CP = mtlb_i(THICK_CP,n,THICK_CP(n-1)+0.0001);
    end;
  end;
  Thickness = interp1(THICK_CP,Thickness_values,RElm);
  Thickness = mtlb_i(Thickness,mtlb_logic(RElm,"<",THICK_CP(1)),mtlb_e(Thickness_values,1));
  Thickness = mtlb_i(Thickness,mtlb_logic(RElm,">",THICK_CP($)),mtlb_e(Thickness_values,$));
  PERCENT_THICKNESS = Thickness;
  DIMENSIONAL_THICKNESS = (CHORD .*PERCENT_THICKNESS) ./100;
 else
   THICK_CP = mtlb_e(Thickness_values,1);
   PERCENT_THICKNESS(mtlb_imp(1,NumSeg),1) = mtlb_e(Thickness_values,1);
   DIMENSIONAL_THICKNESS = (CHORD .*PERCENT_THICKNESS) ./100;
 end;

elseif mtlb_logic(CircleRoot,"==",1) then //Circular root
 ShapeError = 1; //set error = 1 initially

 //Define new control points for circular root chord
 CP_cr1 = RootTranSt*RotorRad;
 CP_cr2 = RootTranEnd*RotorRad;

 if isempty(minRootChord) then %v2_1 = mtlb_e(x,6:10); minChord = mtlb_min(%v2_1,firstnonsingleton(%v2_1));else minChord = minRootChord;end;
 if isempty(maxRootChord) then %v3_1 = mtlb_e(x,6:10); maxChord = mtlb_max(%v3_1,firstnonsingleton(%v3_1));else maxChord = maxRootChord;end;
 if mtlb_logic(minChord,">",maxChord) then minChord = maxChord;end;
 delChord = 0.01; //m

 RtChord = mtlb_imp(mtlb_e(minChord,1),delChord,mtlb_e(maxChord,1))';

 for c = 1:max(size(RtChord))
   RootChord = RtChord(c);
 
   BlendMethod = 2;  //I''ve experimented with different ways of blending the root region, currently this is hardcoded in and not an option for the user to change
   if BlendMethod==1 then
     rad_CP_cr = [CP_cr1;mtlb_a(CP_cr1,CP_cr2)/2;rad_CP];
     Chord_CP_cr = [RootChord;mtlb_a(RootChord,mtlb_e(Chord_CP,1))/2;Chord_CP];
     R_CHORD_CP = rad_CP_cr;
     CHORD_CP = Chord_CP_cr;
   elseif BlendMethod==2 then
     b = 0.7;
     d1 = mtlb_s(CP_cr2,b*mtlb_s(CP_cr2,CP_cr1));
     d2 = mtlb_a(CP_cr1,b*mtlb_s(CP_cr2,CP_cr1));
     Chord_rad_CP_cr = [CP_cr1;d1;d2;rad_CP];
     Chord_CP_cr = [RootChord;RootChord;mtlb_e(Chord_CP,1);Chord_CP];
     R_CHORD_CP = Chord_rad_CP_cr;
     CHORD_CP = Chord_CP_cr;
   end;
 
   B_Chord = Bezier([Chord_rad_CP_cr,Chord_CP_cr],NumSeg);  //returns equal spaced vector
   //Now interpolate to spacing used in RElm
   Chord = interp1(B_Chord(:,1),B_Chord(:,2),RElm);
   Chord = mtlb_i(Chord,isnan(Chord),RootChord);  //only changes values near the root
   CHORD = Chord;
 
   //Now find the location of max chord
   %v4_2 = Chord;  r_ChordMax = mtlb_e(RElm,mtlb_logic(Chord,"==",max(%v4_2,firstnonsingleton(%v4_2))));
   r_ChordMax = r_ChordMax($);  //incase multiple maximums exist, choose the last one
   //         hold on; subplot(4,1,1); plot(R_CHORD_CP,CHORD_CP,''sk'',RElm,CHORD,''x-r''); legend(''Design Variables'',''Chord'');
   %v53 = %f;  if mtlb_logic(r_ChordMax,">",mtlb_e(rad_CP,1)) then %v53 = c~=max(size(RtChord));end;
   if %v53 then
     continue;  //skip to next iteration with larger chord value to see if this makes a feasible blade shape
   end;
   Twist_rad_CP = [r_ChordMax;mtlb_e(rad_CP,2:$)];
   R_TWIST_CP = Twist_rad_CP;
   TWIST_CP = Twist_CP;
   B_Twist = Bezier([R_TWIST_CP,TWIST_CP],NumSeg);  //returns equal spaced vector
   //Now interpolate to spacing used in RElm
   Twist = interp1(B_Twist(:,1),B_Twist(:,2),RElm);
   %v6_2 = Twist;  Twist = mtlb_i(Twist,isnan(Twist),max(%v6_2,firstnonsingleton(%v6_2)));
   TWIST = Twist;
   //         hold on; subplot(4,1,2); plot(R_TWIST_CP,TWIST_CP,''sk'',RElm,TWIST,''x-b''); legend(''Design Variables'',''Pre-Twist'');
 
 
   //Now define the percent thickness distribution
   if max(size(Thickness_values))>1 then
     if mtlb_logic(ThickMethod,"==",1) then //Piecewise Constant
      THICK_CP(1:length(RotorRad*mtlb_e(x,11:10+max(size(Thickness_values))-1)),1) = RotorRad*mtlb_e(x,11:10+max(size(Thickness_values))-1);
      //create a step function
      Thickness(1:max(size(RElm)),1) = mtlb_e(Thickness_values,$);
      a = 1; //intial index value
      for n = 1:max(size(THICK_CP))
        for k = a:max(size(RElm))
          if mtlb_logic(RElm(k),"<",THICK_CP(n)) then
            Thickness(k,1) = mtlb_e(Thickness_values,n);
          elseif mtlb_logic(RElm(k),">",THICK_CP(n)) then
            a = k;
            break
          end;
        end;
      end;
      //Thickness(end+1:length(RElm)) = Thickness_values(end);
      Thickness = mtlb_i(Thickness,mtlb_logic(RElm,"<=",RotorRad*RootTranSt),100);
      dimT = (CHORD .*Thickness) ./100; //Dimensional Thickness
     
     elseif mtlb_logic(ThickMethod,"==",2) then //Piecewise Linear
      THICK_CP(1:length(RotorRad*mtlb_e(x,11:10+max(size(Thickness_values)))),1) = RotorRad*mtlb_e(x,11:10+max(size(Thickness_values)));
      //Make sure that none of the control points are equal, so the interpolation does not crash
      for n = 2:max(size(THICK_CP))
        if mtlb_logic(THICK_CP(n),"<=",THICK_CP(n-1)) then
          THICK_CP = mtlb_i(THICK_CP,n,THICK_CP(n-1)+0.0001);
        end;
      end;
      Thickness = interp1(THICK_CP,Thickness_values,RElm);
      Thickness = mtlb_i(Thickness,mtlb_logic(RElm,"<",THICK_CP(1)),mtlb_e(Thickness_values,1));
      Thickness = mtlb_i(Thickness,mtlb_logic(RElm,">",THICK_CP($)),mtlb_e(Thickness_values,$));
      Thickness = mtlb_i(Thickness,mtlb_logic(RElm,"<=",RotorRad*RootTranSt),100);
      dimT = (CHORD .*Thickness) ./100; //Dimensional Thickness
     end;
   else //only one airfoil was used
    THICK_CP = r_ChordMax;
    Thickness(mtlb_logic(RElm,"<=",RotorRad*RootTranSt),1) = 100;
    Thickness(mtlb_logic(RElm,">",RotorRad*RootTranSt),1) = mtlb_e(Thickness_values,1);
    dimT = (CHORD .*Thickness) ./100; //Dimensional Thickness
   end;
 
   //For the Piecewise Linear case, if the location of Max Chord is inboard of the 1st design
   //variable for thickness, blend to the location of the 1st design variable.  However, for 
   //the Piecewise Constant case, we will blend the location of max chord.
   if mtlb_logic(ThickMethod,"==",2) then F = %f;else F = %t;end;
   if mtlb_logic(r_ChordMax,">",THICK_CP(1))+F then
     t_i = mtlb_logic(RElm,"<=",r_ChordMax) & mtlb_logic(RElm,">=",RotorRad*RootTranSt);
     x2 = r_ChordMax;
     fx2 = mtlb_e(dimT,mtlb_logic(RElm,"==",r_ChordMax));
     // !! L.205: Unknown function FiniteDiff not converted, original calling sequence used.
     DFDX = FiniteDiff(RElm,dimT);
     dfx2 = mtlb_e(DFDX,mtlb_logic(RElm,"==",r_ChordMax));
   else
     t_i = mtlb_logic(RElm,"<=",THICK_CP(1)) & mtlb_logic(RElm,">=",RotorRad*RootTranSt);
     // !! L.209: Unknown function FindInd not converted, original calling sequence used.
     x2 = mtlb_e(RElm,FindInd(RElm,THICK_CP(1)));
     // !! L.210: Unknown function FindInd not converted, original calling sequence used.
     fx2 = mtlb_e(dimT,FindInd(RElm,THICK_CP(1)));
     // !! L.211: Unknown function FiniteDiff not converted, original calling sequence used.
     DFDX = FiniteDiff(RElm,dimT);
     // !! L.212: Unknown function FindInd not converted, original calling sequence used.
     dfx2 = mtlb_e(DFDX,FindInd(RElm,THICK_CP(1)));
   end;
 
   t_vals = mtlb_e(RElm,t_i);
   x1 = RotorRad .*RootTranSt;
   fx1 = RootChord;
   dfx1 = 0;
   dimT = mtlb_i(dimT,t_i,CubeFit2(x1,x2,fx1,fx2,dfx1,dfx2,t_vals));
   Thickness_Corrected = (100 .*dimT) ./Chord;
   PERCENT_THICKNESS = Thickness_Corrected;
   DIMENSIONAL_THICKNESS = (CHORD .*PERCENT_THICKNESS) ./100;
 
   //             hold on; subplot(4,1,3); plot(THICK_CP,Thickness_values,''sk'',RElm,PERCENT_THICKNESS,''x-g''); legend(''Design Variables'',''% Thickness (t/c)'');
   //             hold on; subplot(4,1,4); plot(RElm,DIMENSIONAL_THICKNESS,''x-k''); legend(''Dimensional Thickness (t)'');
 
 
   //Now calculate the slope of the non-dimensional AND dimensional thickness, and make sure
   //that BOTH are monotonically decreasing
   // !! L.230: Unknown function FiniteDiff not converted, original calling sequence used.
   ddimTdr = FiniteDiff(RElm,DIMENSIONAL_THICKNESS);
   // !! L.231: Unknown function FiniteDiff not converted, original calling sequence used.
   dTdr = FiniteDiff(RElm,PERCENT_THICKNESS);
   tol = 0.0001;
   // ! L.233: abs(mtlb_logic(ddimTdr,">",tol)) may be replaced by:
   // !    --> mtlb_logic(ddimTdr,">",tol) if mtlb_logic(ddimTdr,">",tol) is Real.
   // ! L.233: abs(mtlb_logic(dTdr,">",tol)) may be replaced by:
   // !    --> mtlb_logic(dTdr,">",tol) if mtlb_logic(dTdr,">",tol) is Real.
   %v73 = %f;  if mtlb_logic(mtlb_any(abs(mtlb_logic(ddimTdr,">",tol))),"==",%f) then %v73 = mtlb_logic(mtlb_any(abs(mtlb_logic(dTdr,">",tol))),"==",%f);end;
   if %v73 then
     ShapeError = 0;
     //break %end the loop, the dimensional & percent thickness are now both monotonically decreasing, hurrah
     return;
   end;
 
 end;

 //     if ShapeError == 1;     %c == length(RtChord)
 //         %We tried all the possible root cord values, but could not make a valid blade shape
 //             %ShapeError = 1;
 //             RElm = [];
 //             TWIST = [];
 //             CHORD = [];
 //             PERCENT_THICKNESS = [];
 //             DIMENSIONAL_THICKNESS = [];
 //             %Print a warning:
 //             %fprintf(''Warning: Could not satisfy constraints for monotonically decreasing\nthickness distributions. Please check for reasonable inputs under\nthe """"Circular Root"""" section, continuing optimization...\n'');
 //             return
 //     end

end;

// if plot == 1;
// =========================================================================%
// Plots
// figure(''name'',''Initial Population'');
// subplot(4,1,1); plot(R_CHORD_CP,CHORD_CP,''sk'',RElm,CHORD,''x-r''); legend(''Design Variables'',''Chord'');
// subplot(4,1,2); plot(R_TWIST_CP,TWIST_CP,''sk'',RElm,TWIST,''x-b''); legend(''Design Variables'',''Pre-Twist'');
// if ThickMethod == 2
// subplot(4,1,3); plot(THICK_CP,Thickness_values,''sk'',RElm,PERCENT_THICKNESS,''x-g''); legend(''Design Variables'',''% Thickness (t/c)'');
// else
//     for k = 1:length(Thickness_values)-1;
//         tcp(k,1) = 0.5*(Thickness_values(k)+Thickness_values(k+1));
//     end
// subplot(4,1,3); plot(THICK_CP,tcp,''sk'',RElm,PERCENT_THICKNESS,''x-g''); legend(''Design Variables'',''% Thickness (t/c)'');
// end
// subplot(4,1,4); plot(RElm,DIMENSIONAL_THICKNESS,''x-k''); legend(''Dimensional Thickness (t)'');
// =========================================================================%
// end
endfunction
