function [UserInputError] = Error_Check
// Checks for errors in the user input and displays an error message.

global NumBlade NumSeg RotorRad HubHt WatDepth Type HubRad Prated SpdSt SpdEnd SpdDel...
       OmgMin OmgMax rho KinVisc Pv Patm CavSF ProbDist U_mean Weib_k Weib_c NumGens PopSize EliteCount...
       ShaftTilt PreCone CrossFrc Thickness_values family TwistLB TwistUB ChordLB ChordUB...
       ThickMethod ThickLB ThickUB CircleRoot minRootChord maxRootChord RootTranSt RootTranEnd...
       StructuralOpt Einput allowableStrain matDensity SFstruct STmin STdel ParetoFraction OpenFile_Error;
     
UserInputError = 0; // Initially there are no errors
fprintf(1,'\n');
//  NumBlade
if NumBlade < 1 | abs(round(NumBlade) - NumBlade) > 1e-6; // Checks if integer 
   disp('NumBlade must be an integer > 0.');
   UserInputError = 1; fprintf(1,'\n');
end

// NumSeg
if NumSeg < 5 | abs(round(NumSeg) - NumSeg) > 1e-6; // Checks if integer 
   disp('NumSeg must be an integer >= 5.');
   UserInputError = 1; fprintf(1,'\n');
end

// RotorRad
if RotorRad < 0;
    disp('RotorRad must be > 0');
    UserInputError = 1; fprintf(1,'\n');
elseif RotorRad > HubHt & Type == 1;
    disp('The rotor radius is larger than the hub height. The blades will dig into the ground.');
    disp('Increase the hub height or decrease the rotor radius.');
    UserInputError = 1; fprintf(1,'\n');
elseif ((2*RotorRad > WatDepth) |(HubHt + RotorRad > WatDepth)) & Type == 2;
    disp('Either the rotor diameter is larger than the water depth, or the hub height is too tall.');
    disp('The blades will either extend beyond the free surface of the water, or dig into the ground.');
    UserInputError = 1; fprintf(1,'\n');
end

// HubRad
if HubRad < 0 | HubRad >= RotorRad;
    disp('The hub radius must be >= 0, and the hub radius must be less than the rotor radius.');
    UserInputError = 1; fprintf(1,'\n');
end

// Prated
if Prated <= 0 
    disp('The rated power must be greater than 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// ShaftTilt
if ShaftTilt < -180 | ShaftTilt > 180 
    disp('ShaftTilt must be -180 <= ShaftTilT <= 180 (deg).');
    UserInputError = 1; fprintf(1,'\n');
end

// PreCone
if PreCone < -90 | PreCone > 90 
    disp('PreCone must be -90 <= PreCone <= 90 (deg).');
    UserInputError = 1; fprintf(1,'\n');
end

// SpdSt
if SpdSt < 0
    disp('SpdSt must be >= 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// SpdEnd
if SpdEnd <= SpdSt
    disp('SpdEnd must be > SpdSt.');
    UserInputError = 1; fprintf(1,'\n');
end
    
// SpdDel
int = (SpdEnd - SpdSt)/SpdDel;
if SpdDel <= 0
    disp('SpdDel must be > 0.');
    UserInputError = 1; fprintf(1,'\n');
elseif SpdDel >= (SpdEnd - SpdSt);
    disp('SpdDel must be < (SpdEnd - SpdSt).');
    UserInputError = 1; fprintf(1,'\n');
elseif abs(round(int) - int) > 1e-6; // Checks if integer
    disp('The flow speeds are not evenly spaced.')  
    disp('Choose a new value for SpdSt or SpdEnd or SpdDel which results in an evenly spaced interval.');
    UserInputError = 1; fprintf(1,'\n');
end

// OmgMin
if OmgMin <= 0
    disp('OmgMin must be > 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// OmgMax
if OmgMax < OmgMin
    disp('OmgMax must be >= OmgMin.');
    UserInputError = 1; fprintf(1,'\n');
elseif OmgMax <= 0;
    disp('OmgMax must be > 0');
    UserInputError = 1; fprintf(1,'\n');
elseif OmgMin == 0 & OmgMax == 0;
    disp('Increase OmgMax to a non-zero value.');
    UserInputError = 1; fprintf(1,'\n');
end

// Rho
if rho <= 0
    disp('The fluid density Rho must be > 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// KinVisc
if KinVisc <= 0;
    disp('The kinematic fluid viscosity KinVisc must be > 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// Pv         
if Pv < 0;
    disp('Pv must be >= 0.');
    UserInputError = 1; fprintf(1,'\n');
elseif Pv >= Patm;
    disp('Pv must be < Patm.');
    UserInputError = 1; fprintf(1,'\n');
end

// Patm
if Patm <= 0
    disp('Patm must be > 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// CavSF
if CavSF < 0
    disp('CavSF must be >= 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// U_mean
if ProbDist == 1 & U_mean <= 0
    disp('The Rayleigh distribution mean flow speed must be > 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// k and c
if ProbDist == 2 & (Weib_k <= 0 | Weib_c <= 0)
    disp('The Weibull shape factor, k, and the scale factor, c, must both be > 0.');
    UserInputError = 1; fprintf(1,'\n');
end

// NumGens
if NumGens < 1;
    disp('NumGens must be >= 1.');
    UserInputError = 1; fprintf(1,'\n');
end

// PopSize and EliteCount
if PopSize <= 0 | PopSize <= EliteCount
    disp(['PopSize must be > 0, and also PopSize must be > ' num2str(EliteCount)'.']);
    UserInputError = 1; fprintf(1,'\n');
end

// CrossFrc
if CrossFrc < 0 | CrossFrc > 1;
    disp('CrossFrc must be between 0 and 1 (inclusive).');
    UserInputError = 1; fprintf(1,'\n');
end

//  Thickness_values
if any(Thickness_values<=0)
    disp('The available % Thick values must be > 0');
    UserInputError = 1; fprintf(1,'\n');
end

// need to make sure we can open all the desired airfoil files, and
// check if they even exist
if CircleRoot == 1
    ThickVals = [100 Thickness_values];
else
    ThickVals = Thickness_values;
end

for n = 1:length(ThickVals);
   if ThickVals(n) >= 99.95;
        Thick_Suffix = '_1000';
      elseif ThickVals(n) >= 9.95 & ThickVals(n) < 99.95;
        Thick_Suffix = ['_0' num2str(10*ThickVals(n),'%3.0f')];
      else
        Thick_Suffix = ['_00' num2str(10*ThickVals(n),'%2.0f')];
   end;
   
   filename1 =  ['Input_Files\Airfoil_Data\' family Thick_Suffix '.dat'];
   fid1 = fopen(filename1);
   if fid1 == -1
       disp(['Could not find the airfoil file ' filename1]);
       disp('Make sure you have named your airfoil files according to the proper convention,');
       disp('and that you have not used an illegal family name (FoilFam). Also, if you are');
       disp('using a circular root, make sure you have an airfoil file for a circle and it');
       disp('has the correct filename.');
       UserInputError = 1; fprintf(1,'\n');
       
   else // if the file was able to be opened, find out how many Reynolds tables and get the Stall Angle
   NumTables(n,1) = cell2mat(textscan(fid1,'%f','HeaderLines',3)); // number of tables for Reynolds number in the file
   frewind(fid1);
   StallAoA(n,1) = cell2mat(textscan(fid1,'%f','HeaderLines',5));
   end
   
   if StructuralOpt == 1; // check if the profile file exists
   filename2 =  ['Input_Files\Airfoil_Data\' family Thick_Suffix '.prof'];
   fid2 = fopen(filename2);
       if fid2 == -1
       disp(['Could not find the airfoil file ' filename2]);
       disp('For stuctural optimization problems, you need to have an airfoil file for each');
       disp('airfoil being used, where each file contains the x-y coordinates for the airfoil.');
       disp('Make sure you have named your airfoil files according to the proper convention,');
       disp('and that you have not used an illegal family name (FoilFam). Also, if you are');
       disp('using a circular root, make sure you have an airfoil file for a circle and it');
       disp('has the correct filename.');
       UserInputError = 1; fprintf(1,'\n');        
       end
   end
   
end

// make sure all airfoil files have the same number of Reynolds tables
if fid1 ~= -1;
    if length(NumTables)>1
        for n = 2:length(NumTables)
            if NumTables(n) ~= NumTables(n-1)
                disp('Airfoil files have mismatching number of Reynolds number tables.');
                disp('Make sure each airfoil file has the same number of Reynolds number tables.');
                UserInputError = 1; fprintf(1,'\n');
                break
            end
        end
    end

    // make sure the Stall AoA have been entered, display a warning if a value of
    // 0 is found, otherwise assume the user has input the correct stall AoA.
    // Do not check the circle airfoil file, is will always have StallAoA=0.
    if CircleRoot == 0; N = 1; elseif CircleRoot == 1; N = 2; end;    
    for n = N:length(StallAoA)
        if StallAoA(n) == 0
            disp('One or more airfoil files have a value of 0 entered for the');
            disp('stall angle of attack. If this is not the correct value (or you');
            disp('forgot to input a value for stall angle of attack), please stop');
            disp('HARP_Opt and enter the correct value in the airfoil file.');
            disp('HARP_Opt is continuing to run ...');
        end
    end
end

//  filename_main // no spaces or special charcters are allowed in this string

// Make sure all bounds for Twist and Chord have 5 values
if length(TwistLB)<5 | length(TwistUB)<5 | length(ChordLB)<5 | length(ChordUB)<5
    disp('Please enter 5 values in each field for TwistLB, TwistUB, ChordLB, and ChordUB.');
    UserInputError = 1; fprintf(1,'\n');
    return
end

//  TwistLB and TwistUB
if any(abs(TwistLB)>90) | any(abs(TwistUB)>90)
    disp('Twist values must be between -90 and 90 degrees. Change the upper & lower twist bounds to be within this range.');
    UserInputError = 1; fprintf(1,'\n');
end

//  ChordLB and ChordUB
if any(ChordLB<=0) | any(ChordUB<=0)
    disp('Cannot have a chord length = 0 m. Change the upper & lower chord bounds to non-zero numbers.');
    UserInputError = 1; fprintf(1,'\n');
end

// Make sure all bounds for Twist and Chord are monotonically decreasing
for n = 2:5
    if TwistLB(n) > TwistLB(n-1) | TwistUB(n) > TwistUB(n-1) | ChordLB(n) > ChordLB(n-1) | ChordUB(n) > ChordUB(n-1)
    disp('The values for TwistLB, TwistUB, ChordLB, and ChordUB must all be');
    disp('monotonically decreasing (or be equal) from the blade hub to tip.');
    UserInputError = 1; fprintf(1,'\n');
    end
end 
    
//  ThickLB and ThickUB
if any(ThickLB<0) | any(ThickLB>1) | any(ThickUB<0) | any(ThickUB>1)
    disp('Bounds for the airfoil/hydrofoil positions must be between 0 and 1.');
    disp('Values of NaN are also accepted, refer to the HARP_Opt user''s guide');
    disp('for instructions on how to correctly specify bounds for the airfoils.');
    UserInputError = 1; fprintf(1,'\n');
end
// make sure the correct number of bounds were entered for the //  thickness distributions
if length(Thickness_values) == 1
    if isempty(ThickLB)==0 | isempty(ThickUB)==0
    disp('You have only input a single airfoil/hydrofoil, bounds are not needed.');
    disp('Please clear the (r/R) LB and (r/R) UB fields.');
    UserInputError = 1; fprintf(1,'\n');
    end
elseif length(Thickness_values) > 1
  if ThickMethod == 1; // Piecewise Constant
    if (length(ThickLB)~=length(Thickness_values)-1 | length(ThickUB)~=length(Thickness_values)-1) & (isempty(ThickLB)==0 | isempty(ThickUB)==0)
    disp('For the Piecewise Constant percent thickness distribution, you need to input');
    disp('N-1 (r/R) upper bounds and N-1 (r/R) lower bounds for N airfoils (N>=1).');
    disp('Refer to the HARP_Opt user''s guide for instructions on how to correctly');
    disp('specify bounds for the aifoils.');
    UserInputError = 1; fprintf(1,'\n');
    end
  elseif ThickMethod == 2; // Piecewise Linear
    if (length(ThickLB)~=length(Thickness_values) | length(ThickUB)~=length(Thickness_values)) & (isempty(ThickLB)==0 | isempty(ThickUB)==0)
    disp('For the Piecewise Linear percent thickness distribution, you need to input');
    disp('N (r/R) upper bounds and N (r/R) lower bounds for N airfoils (N>=1).');
    disp('Refer to the HARP_Opt user''s guide for instructions on how to correctly');
    disp('specify bounds for the aifoils.');
    UserInputError = 1; fprintf(1,'\n');
    end  
  end   
end

// If thickness bounds are entered, make sure the bounds do not overlap
if all(ThickLB==0) & all(ThickUB==1)
    // do nothing, the user left the fields blank and then the bounds were
    // set automatically to 0 and 1
else // the user has entered bounds other than 0 and 1, make sure they don't overlap
    if ~isempty(ThickLB) & ~isempty(ThickUB) & length(Thickness_values) > 1
       for n = 1:length(ThickLB)-1
           if  any(ThickUB(n) > ThickLB(n+1:end)) | any(ThickLB(n) > ThickUB(n+1:end))           
           disp('The bounds for (r/R) LB and (r/R) UB are overlapping, make sure');
           disp('these bounds do not overlap. Refer to the HARP_Opt user''s guide for');
           disp('instructions on how to correctly specify bounds for the airfoils.');
           UserInputError = 1; fprintf(1,'\n');
           break
           end
       end
       
       for n = 1:length(ThickLB)
           if  ThickLB(n) > ThickUB(n)           
           disp('The bounds for (r/R) LB and (r/R) UB are overlapping, make sure');
           disp('these bounds do not overlap. Refer to the HARP_Opt user''s guide for');
           disp('instructions on how to correctly specify bounds for the airfoils.');
           UserInputError = 1; fprintf(1,'\n');
           break
           end
       end
       
      for n = 2:length(ThickLB)
           if  ThickLB(n) < ThickUB(n-1) | ThickUB(n) < ThickUB(n-1)           
           disp('The bounds for (r/R) LB and (r/R) UB are overlapping, make sure');
           disp('these bounds do not overlap. Refer to the HARP_Opt user''s guide for');
           disp('instructions on how to correctly specify bounds for the airfoils.');
           UserInputError = 1; fprintf(1,'\n');
           break
           end
       end
       
    end 
end

if CircleRoot == 1
   if minRootChord<0
       disp('RtChordMin must be > 0. Alternatively, you can leave this value blank');
       disp('and HARP_Opt will set RtChordMin to the smallest chord value of the');
       disp('current individual (meaning this values changes every iteration).');
       UserInputError = 1; fprintf(1,'\n');
   end
   
   if maxRootChord < minRootChord
       disp('RtChordMax must be > RtChordMin. Alternatively, you can leave this value');
       disp('blank and HARP_Opt will set RtChordMax to the largest chord value of the');
       disp('current individual (meaning this values changes every iteration).');
       UserInputError = 1; fprintf(1,'\n'); 
   end
   
   if RootTranSt<HubRad/RotorRad | RootTranSt>=RootTranEnd
       disp('RtTranSt must be >= HubDia/RotorDia and < RtTranEnd');
       UserInputError = 1; fprintf(1,'\n');
   end
   
   if RootTranEnd > 1 | RootTranEnd < RootTranSt
       disp('RtTranEnd must be <= 1 and > RtTranSt.');
       UserInputError = 1; fprintf(1,'\n');
   end    
end

if StructuralOpt == 1;
    if Einput <= 0
        disp('E must be > 0.');
        UserInputError = 1; fprintf(1,'\n');
    end
    
    if allowableStrain <= 0
        disp('MaxStrain must be > 0.');
        UserInputError = 1; fprintf(1,'\n');
    end
    
    if matDensity <= 0
        disp('MatDensity must be > 0.');
        UserInputError = 1; fprintf(1,'\n');
    end
    
    if SFstruct <= 0
        disp('SF must be > 0.');
        UserInputError = 1; fprintf(1,'\n');
    end
    
    if STmin <= 0
        disp('STmin must be > 0.');
        UserInputError = 1; fprintf(1,'\n');
    end
    
    if STdel <= 0
        disp('STdel must be > 0.');
        UserInputError = 1; fprintf(1,'\n');
    end
    
    if ParetoFraction <= 0 | ParetoFraction >= 1
       disp('ParetoFrac must be between 0 and 1')
       UserInputError = 1; fprintf(1,'\n');
    end        
end

// make sure the flow distribution file was opened correctly
if ProbDist == 3 & OpenFile_Error == 1
    disp('The file containing user-defined flow distribution data could not be opened.')
    UserInputError = 1; fprintf(1,'\n');
end
