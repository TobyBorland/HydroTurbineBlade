function [F] = Main(x)

// Output variables initialisation (not found in input variables)
F=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This is the main function of HARP_Opt, and is the function which is called
//upon by the genetic algorithm.  Its purpose is simple: accept an
//individual from the genetic algorithm and then return its fitness value.
//However this function must perform many different tasks in a very specific
//sequence in order to do this.  This function performs the following tasks:
//   -Builds the WT_Perf input files for either Parametric or Combined Case
//    analysis
//   -calls WT_Perf to evaluate the WT_Perf input files
//   -reads in output data from the WT_Perf output files
//   -determines the optimal rotor speeds for variable speed cases
//   -determines the optimal blade pitch for variable pitch cases
//   -writes data to the *Detailed_GA_Output.dat file


//=========================================================================%



global("SpdSt","SpdEnd","SpdDel","OmgMin","OmgMax","Prated","rho","Patm","Pv","HubHt","WatDepth","Type","NumSeg","HubRad","Thickness_values","RotorRad","RootDir","filename_main","NumCases","SpdCtrl","PitCtrl","ProbDist","U_mean","Weib_k","Weib_c","CavSF","NumVars","ThickMethod","StructuralOpt","RecordFailures");

// !! L.22: Matlab function persistent not yet converted, original calling sequence used.
persistent("FailCount")
//=========================================================================%

//=========================================================================%
// This section only used for debugging, comment out when using the GUI
// x=[25 10 4 3 2 0.7 0.4 0.3 0.2 0.1 0.4 0.8 32]''; %ThickMethod=1; SpdCtrl=0
// x=[25 10 4 3 2 0.7 0.4 0.3 0.2 0.1 0.4 0.8]''; %ThickMethod=1; SpdCtrl=1
// x=[25 10 4 3 2 0.7 0.4 0.3 0.2 0.1 0.4 0.6 0.8 32]''; %ThickMethod=2; SpdCtrl=0
// % x=[25 10 4 3 2 0.7 0.4 0.3 0.2 0.1 0.4 0.6 0.8]''; %ThickMethod=2; SpdCtrl=1
// filename_main = ''Debugging'';
// family = ''FFA'';
// Thickness_values = [30 24 21]'';
// NumSeg = 40;
// Prated = 72;    %kW
// RotorRad = 2.5; %m
// HubRad = 0.2;   %m
// HubHt = 5;      %m
// WatDepth = 10;  %m
// SpdSt = 0.2;    %m/s
// SpdEnd = 3.0;   %m/s
// SpdDel = 0.1;   %m/s
// OmgMin = 10;    %rpm
// OmgMax = 60;    %rpm
// rho = 1025;     %kg/m^3
// Patm = 101325;  %Pa
// Pv = 2500;      %Pa
// CavSF = 1.0;    %Safety factor for cavitation number    
// Type = 2;       %1 for Wind Turbine, 2 for Hydrokinetic Turbine
// Correct_3D = 0; %0 for no-3Dcorrections, 1 for Selig-Du Lift & Drag, 2 for Selig-Du Lift & Eggars Drag
// SpdCtrl = 0;    %0 for fixed rotor speed, 1 for variable rotor speed
// PitCtrl = 0;    %0 for fixed blade pitch, 1 for variable blade pitch
// ProbDist = 3;   %0 don''t calc. AEP, 1 Rayleigh, 2 Weibull, 3 Custom
// ThickMethod = 1; %1 for Piecewise Constant, 2 for Piecewise Linear
//=========================================================================%

//Set errors = 0 initially
ShapeError = 0;
WTP_Error = 0;
BED_Error = 0;
Cav_Error = 0;

NumCases = round(mtlb_a(mtlb_s(SpdEnd,SpdSt)/SpdDel,1));
if isnan(NumCases) then NumCases = 1;end;

V = mtlb_imp(SpdSt,SpdDel,SpdEnd)';

//Define the appropriate variables for the AEP calculation
%v02 = %t;if ~mtlb_logic(ProbDist,"==",0) then %v02 = mtlb_logic(ProbDist,"==",3);end;
if %v02 then //Either dont calculate AEP or use custom distribution
 pUvars = [];
elseif mtlb_logic(ProbDist,"==",1) then //Rayleigh Distribution
 pUvars = U_mean;
elseif mtlb_logic(ProbDist,"==",2) then //Weibull Distribution
 pUvars = [Weib_k,Weib_c];
end;

if mtlb_logic(SpdCtrl,"==",0) then //Fixed Speed: rotor speed becomes the last variable
 flag = 1; //WT_Perf will use Parametric Analysis
 CC = []; //The combined case variable is empty since this is a Parametric Analysis
 Omg = mtlb_e(x,NumVars); //Rotor speed is the last variable
else
  flag = 2;  //WT_Perf will use Combined Case Analysis
  Omg = 1;  //A dummy value is used for the Variable Speed Case
end;


//=========================================================================%
//This section defines the blade geometery
[ShapeError,RElm,Twist,Chord,Thickness] = Define_Blade_Shape(x);
if ShapeError==1 then
  F(1,1) = %inf;
  if mtlb_logic(StructuralOpt,"==",1) then F(1,2) = %inf;end;
  //fprintf(''Shape Error\n'');
  if mtlb_logic(RecordFailures,"==",0) then return;end;
end;
//=========================================================================%

//========================  Variable Speed cases ==========================%
// This section of the code is for the variable speed cases.  
// First the Cp vs TSR curve is calculated to determine the optimal TSR, and
// then the optimal rotor speeds are calculated as a function of flow speed.

if mtlb_logic(SpdCtrl,"==",1) then
  //Build the input file to calculate the Cp vs. TSR curve
  //Build_Input(flag,RElm,Twist,Chord,Thickness,SpdSt,SpdEnd,SpdDel,Omg,PitCtrl,CC)
  // !! L.106: Unknown function Build_Input not converted, original calling sequence used.
  Build_Input(0,RElm,Twist,Chord,Thickness,1,10,0.05,Omg,0,[])

  cd(RootDir+"\Output_Files\"+filename_main);  //Change the working directory to where the WT_Perf executable is located
  // !! L.109: Matlab function evalc not yet converted, original calling sequence used.
  evalc("!Wt_Perf_Jan1 "+filename_main+".wtp");
  cd(RootDir);  //Change the working directory back to the root directory

  // Reads in the Cp vs TSR curve from WT_Perf output file
  // Don''t check for WT_Perf errors here because it''s still possible that the
  // Power Curve may not produce errors, even if the Cp vs. TSR curve does
  fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".oup","rt");
  // !! L.116: Unknown function Read_WTP_Output not converted, original calling sequence used.
  CpvsTSR = Read_WTP_Output(fid,1,0);
  mclose(fid);


  // Finds the optimal Tip Speed Ratio (TSR at which the Power Coefficient is maximum)
  // Truncates the Cp values to only 3 decimal points and then choose the
  // TSR that results in the slowest rotor speed
  e = 10 .^3;
  CpvsTSR(:,2) = fix(CpvsTSR(:,2) .*e) ./e;  //truncates to only 3 decimal points
  %v0_1 = CpvsTSR(:,2);  [CpMax,Cp_ind] = mtlb_max(%v0_1,firstnonsingleton(%v0_1));
  TSR_Opt = CpvsTSR(Cp_ind,1);

  OmgCC = ((TSR_Opt .*V) .*30)/(%pi .*RotorRad);
  for n = mtlb_imp(1,NumCases)
    if mtlb_logic(OmgCC(n,1),"<",OmgMin) then
      OmgCC(n,1) = OmgMin;
    elseif mtlb_logic(OmgCC(n,1),">",OmgMax) then
      OmgCC(n,1) = OmgMax;
    end;
    //TSRCC(n,1) = OmgCC(n)*RotorRad*pi/(30*V(n));
    //ii = FindInd(CpvsTSR(:,1),TSRCC(n));
    //PwrCC(n,1) = 0.5 .*rho.*pi.*RotorRad^2 .*V(n).^3 .*CpvsTSR(ii,2)./1000;
  end;
  //OmgCC(OmgCC<OmgMin) = OmgMin; The FOR loop above actually performs faster than these two simple lines of code...interesting
  //OmgCC(OmgCC>OmgMax) = OmgMax;

  PitCC(mtlb_imp(1,NumCases),1) = 0;  //Set the Pitch = 0 for all flow speeds, the optimal blade pitches are determined after the optimization is complete
  CC = [V,OmgCC,PitCC];  //Define the Combined Case analysis for WT_Perf

  //================== Variable Speed Fixed Pitch case ==============%
  //Variable Speed Fixed Pitch case
  %v13 = %f;  if mtlb_logic(SpdCtrl,"==",1) then %v13 = mtlb_logic(PitCtrl,"==",0);end;
  if %v13 then
    //Builds the WT_Perf input file to calculate the Power Curve (P vs. V) 
    // !! L.149: Unknown function Build_Input not converted, original calling sequence used.
    Build_Input(flag,RElm,Twist,Chord,Thickness,SpdSt,SpdEnd,SpdDel,Omg,0,CC)
    cd(RootDir+"\Output_Files\"+filename_main);  //Change the working directory to where the WT_Perf executable is located
    // !! L.151: Matlab function evalc not yet converted, original calling sequence used.
    evalc("!Wt_Perf_Jan1 "+filename_main+".wtp");
    cd(RootDir);  //Change the working directory back to the root directory
  
    //Reads in the Power vs. Velocity curve from the WT_Perf output file
    fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".oup","rt");
    // !! L.156: Unknown function Read_WTP_Output not converted, original calling sequence used.
    [PvsV,WTP_Error] = Read_WTP_Output(fid,1,1);
    mclose(fid);
    if mtlb_logic(WTP_Error,"==",1) then
      F = mtlb_i(F,1,%inf);
      if mtlb_logic(StructuralOpt,"==",1) then F = mtlb_i(F,2,%inf);end;
      //             fprintf(''WT_Perf Error line159\n'');
      if mtlb_logic(RecordFailures,"==",0) then return;end;
    end;
  
    //Find Vrated,  and then change the rotor speed to remain constant above Vrated
    //ii = PvsV(:,2)>=Prated; %The FOR loop is faster than this line of code      
    ii(mtlb_imp(1,NumCases),1) = %f;
    for n = mtlb_imp(1,NumCases)
      if mtlb_logic(PvsV(n,2),">=",Prated) then
        ii = mtlb_i(ii,n,%t);
        Vrated_ind = n;
        Vrated = V(n);
        V_reg3 = V(n:$);
        break
      end;
    end;
  
    //If the power exceeded Prated, regulate power by reducing the rotor speed
    // ! L.179: abs(ii) may be replaced by:
    // !    --> ii if ii is Real.
  
    if mtlb_any(abs(ii)) then
      //Find the required power coefficients and corresponding rotor speeds to maintain constant rated power above Vrated
      req_Cp = (1000*Prated) ./((((0.5*rho)*%pi)*(RotorRad^2)) .*(V_reg3 .^3));
    
      for n = 1:max(size(req_Cp))
        Cp_Diff = abs(CpvsTSR(mtlb_imp(1,Cp_ind),2)-req_Cp(n));
        %v2_1 = Cp_Diff;  [Diff_min,jj] = min(%v2_1,firstnonsingleton(%v2_1));
        req_TSR = CpvsTSR(jj,1);
        Omg_reg3 = ((req_TSR*V_reg3(n))*30)/(%pi*RotorRad);
        if mtlb_logic(Omg_reg3,"<",OmgMin) then
          Omg_reg3 = OmgMin;
        elseif mtlb_logic(Omg_reg3,">",OmgMax) then
          Omg_reg3 = OmgMax;
        end;
        OmgCC(Vrated_ind-1+n,1) = Omg_reg3;
      end;
    end;
  end;
  //============= end of Variable Speed Fixed Pitch case ============%

  //The optimal rotor speeds have now been calculated
  CC = [V,OmgCC,PitCC];  //Define the Combined Case analysis for WT_Perf      
end;
//==================== end of Variable Speed cases ========================%

//=========================================================================%
//Builds the WT_Perf input file to calculate the Power Curve (P vs. V)
//Build_Input(flag,RElm,Twist,Chord,Thickness,SpdSt,SpdEnd,SpdDel,Omg,PitCtrl,CC)
// !! L.207: Unknown function Build_Input not converted, original calling sequence used.
Build_Input(flag,RElm,Twist,Chord,Thickness,SpdSt,SpdEnd,SpdDel,Omg,0,CC)
cd(RootDir+"\Output_Files\"+filename_main);//Change the working directory to where the WT_Perf executable is located
// !! L.209: Matlab function evalc not yet converted, original calling sequence used.
evalc("!Wt_Perf_Jan1 "+filename_main+".wtp");
cd(RootDir);//Change the working directory back to the root directory

//Reads in the Power vs. Velocity curve from WT_Perf output file
fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".oup","rt");
if mtlb_logic(SpdCtrl,"==",1) then //This is a Combined Case Analysis
 // !! L.215: Unknown function Read_WTP_Output not converted, original calling sequence used.
 [PvsV,WTP_Error] = Read_WTP_Output(fid,1,1);
else //This is a Parametric Analysis
 // !! L.217: Unknown function Read_WTP_Output not converted, original calling sequence used.
 [PvsV,WTP_Error] = Read_WTP_Output(fid,1,0);
end;
mclose(fid);
if mtlb_logic(WTP_Error,"==",1) then
  F = mtlb_i(F,1,%inf);
  if mtlb_logic(StructuralOpt,"==",1) then F = mtlb_i(F,2,%inf);end;
  //        fprintf(''WT_Perf Error line221\n'');
  if mtlb_logic(RecordFailures,"==",0) then return;end;
end;

//Check if Power is negative, if so change Power value to 0
indices = mtlb_logic(PvsV,">",0);
PvsV = PvsV .*indices;

//Check to see if the Betz law has been violated, sometimes the GA can
//find solutions which cause WT_Perf to go bonkers...resulting in
//efficiencies greater than Betz limit
eff = (1000 .*PvsV(:,2)) ./((((0.5 .*rho) .*%pi) .*(RotorRad .^2)) .*(V .^3));
// ! L.235: abs(mtlb_logic(eff,">",0.593)) may be replaced by:
// !    --> mtlb_logic(eff,">",0.593) if mtlb_logic(eff,">",0.593) is Real.

if mtlb_any(abs(mtlb_logic(eff,">",0.593))) then
  WTP_Error = 1;
end;

Torque = 0;//Dummy value for Torque
//For the Variable Speed Fixed Pitch case, calculate the Torque values
%v02 = %f;if mtlb_logic(SpdCtrl,"==",1) then %v02 = mtlb_logic(PitCtrl,"==",0);end;
if %v02 then
  Torque = PvsV(:,2) ./((CC(:,2) .*%pi)/30);  //Torque
end;

//Read in BED file and check for cavitation
%v02 = %t;if ~mtlb_logic(Type,"==",2) then %v02 = mtlb_logic(StructuralOpt,"==",1);end;
if %v02 then
  fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".bed","rt");
  // !! L.248: Unknown function Read_BED not converted, original calling sequence used.
  [BED_Error,BED,AoA,AFang,LocVel,Cl,Cd,CpMin,Thrust] = Read_BED(fid,NumCases,NumSeg,1);
  mclose(fid);

  Cav_Error = Cavitation_Check(RElm,LocVel,CpMin,1);
end;

%v02 = %t;if ~mtlb_logic(WTP_Error,"==",1) then %v02 = mtlb_logic(BED_Error,"==",1);end;%v12 = %t;if ~%v02 then %v12 = Cav_Error==1;end;
if %v12 | ShapeError==1 then
  F = mtlb_i(F,1,%inf);  //This number needs to be REALLY large, especially for multi-megawatt rated powers, else the calculated fitness value could actually be larger than this penalty value
  if mtlb_logic(StructuralOpt,"==",1) then F = mtlb_i(F,2,%inf);end;
  //         fprintf(1,''Errors exist line 255\n'');

  //If recoding the failed cases, copy files to special folder
  if mtlb_logic(RecordFailures,"==",1) then
    // ! L.261: mtlb(FailCount) can be replaced by FailCount() or FailCount whether FailCount is an M-file or not.
  
    if isempty(mtlb(FailCount)) then
      FailCount = 1;
    end;
  
    FailPrefix = "";
    if mtlb_logic(WTP_Error,"==",1) then
      FailPrefix = FailPrefix+"o";
    end;
    if mtlb_logic(BED_Error,"==",1) then
      FailPrefix = FailPrefix+"b";
    end;
    if Cav_Error==1 then
      FailPrefix = FailPrefix+"c";
    end;
    if ShapeError==1 then
      FailPrefix = FailPrefix+"s";
    end;
  
    // !! L.279: string output can be different from Matlab num2str output.
    FailName = string(FailCount)+"_"+FailPrefix+"_"+filename_main;
  
    //copy the failed .wtp and .oup files to the failed directory
    // !! L.282: Matlab function copyfile not yet converted, original calling sequence used.
    // L.282: (Warning name conflict: function name changed from copyfile to %copyfile).
    %copyfile(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".wtp",RootDir+"\Output_Files\"+filename_main+"\Failed_Cases\"+FailName+".wtp");
    // !! L.283: Matlab function copyfile not yet converted, original calling sequence used.
    // L.283: (Warning name conflict: function name changed from copyfile to %copyfile).
    %copyfile(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".oup",RootDir+"\Output_Files\"+filename_main+"\Failed_Cases\"+FailName+".oup");
  
    // !! L.285: No Scilab equivalent for: Matlab exist(), mtlb_exist() is called.
    // !! L.285: mtlb_exist() is just a partial emulation of Matlab exist() function.
  
    if exist(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".bed","file")==2 then
      //.bed file exists, so copy it also
      // !! L.287: Matlab function copyfile not yet converted, original calling sequence used.
      // L.287: (Warning name conflict: function name changed from copyfile to %copyfile).
      %copyfile(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".bed",RootDir+"\Output_Files\"+filename_main+"\Failed_Cases\"+FailName+".bed");
    end;
  
    FailCount = FailCount+1;
    return;
  end;

else //there were no errors, all is good!
 F = mtlb_i(F,1,Fitness_Function(PvsV,Torque,pUvars)); //Calculate the fitness value using the fitness function
 if mtlb_logic(StructuralOpt,"==",1) then
   // !! L.297: Unknown function Struct_Analysis not converted, original calling sequence used.
   Output = Struct_Analysis(RElm,Twist,Chord,Thickness,AoA,AFang,LocVel,Cl,Cd,Thrust);
   F = mtlb_i(F,2,mtlb_e(Output,"TotalMass"));
 end;
end;

// % % % %     %Print GA data to text file
// % % % %     %need to make a string that tells fprintf to print the right length format
// % % % %     s = [];
// % % % %     for n = 1:NumVars
// % % % %         a = [37    51    46    52   102    92   116]; %ASCII representation of: %3.4f\t
// % % % %         s = [s a];
// % % % %     end
// % % % %     %s(end) = 110; %change last character to """"n"""", (next line statment = \n)
// % % % %     format = char(s);
// % % % %     fid = fopen([RootDir ''\Output_Files\'' filename_main ''\'' filename_main ''_Detailed_GA_Output.dat''],''At'');
// % % % %     fprintf(fid,[''%3.4f\t'' format ''\n''],F(1),x'');    
// % % % %     fclose(fid);
//=========================================================================%


//=========================================================================%
//The following code executes ONLY after the Genetic Algorithm has finished%
//=========================================================================%    
//the length of x was changed to """"trick"""" main.m into running these lines of
//code, not very elegant...but it works for now
if mtlb_logic(max(size(x)),">",NumVars) then

  //the last value in x was a dummy value, get rid of the last dummy value in the vector x
  NumSheet = mtlb_e(x,$);
  x = mtlb_i(x,$,[]);

  //======================== Variable Pitch cases ========================%
  if mtlb_logic(PitCtrl,"~=",0) then
    //Hold the rotor speed constant above Vrated
    ii(mtlb_imp(1,NumCases),1) = %f;  //false is a logical = 0, not a string or variable
    for n = mtlb_imp(1,NumCases)
      if mtlb_logic(PvsV(n,2),">=",Prated) then
        ii = mtlb_i(ii,n,%t);  //true is a logical = 1, not a string or variable
        Vrated_ind = n;
        Vrated = V(n);
        V_reg3 = V(n:$);
      
        if mtlb_logic(SpdCtrl,"==",1) then
          OmgCC(Vrated_ind:$,1) = mtlb_e(OmgCC,Vrated_ind);
          CC(:,2) = OmgCC;  //Add the new rotor speeds into the Combined Case
        end;
      
        break
      else
        Vrated_ind = NumCases;
        Vrated = V($);  //Vrated was never acheived, set Vrated as the cut-out speed
        V_reg3 = Vrated;
      end;
    end;
  
    //Regulate the pitch for flow speeds >= Vrated     
    //Build a Parametric Analysis input file with Varying Pitch and Constant Rotor Speed
    if mtlb_logic(SpdCtrl,"==",0) then //Fixed Speed
     Omg = mtlb_e(x,NumVars);
     CC(:,1) = V;
     CC(:,2) = Omg;
    else //Variable Speed   
     Omg = CC($,2); //use the constant rotor speed from the Combined Case 
    end;
    // !! L.361: Unknown function Build_Input not converted, original calling sequence used.
    Build_Input(1,RElm,Twist,Chord,Thickness,Vrated,SpdEnd,SpdDel,Omg,PitCtrl,[])
    cd(RootDir+"\Output_Files\"+filename_main);  //Change the working directory to where the WT_Perf executable is located
    // !! L.363: Matlab function evalc not yet converted, original calling sequence used.
    evalc("!Wt_Perf_Jan1 "+filename_main+".wtp");
    cd(RootDir);  //Change the working directory back to the root directory
  
    fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".oup","rt");
    // !! L.367: Unknown function Read_WTP_Output not converted, original calling sequence used.
    PwrPit = Read_WTP_Output(fid,401,0);  //Should WT_Perf output errors be checked for here?  How would you correct them at this point?
    mclose(fid);
  
    // ! L.370: abs(mtlb_logic(PvsV(:,2),">",Prated)) may be replaced by:
    // !    --> mtlb_logic(PvsV(:,2),">",Prated) if mtlb_logic(PvsV(:,2),">",Prated) is Real.
  
    if mtlb_any(abs(mtlb_logic(PvsV(:,2),">",Prated))) then
      //Find the pitch angles which result in constant power output above Vrated
      PwrDiff = abs(mtlb_s(Prated,PwrPit(:,2:$)));
      PitCases = round(mtlb_a(mtlb_s(SpdEnd,Vrated)/SpdDel,1));
      for n = mtlb_imp(1,PitCases)
        %v0_1 = PwrDiff(n,:);  [min_PwrDiff,index] = min(%v0_1,firstnonsingleton(%v0_1));
        if max(size(index))>1 then
          index = index(1);
        end;
        if mtlb_logic(PitCtrl,"==",1) then
          PitCC(mtlb_s(mtlb_a(Vrated_ind,n),1),1) = mtlb_a(0,0.1*mtlb_s(index,1));  //PitSt + PitDel(index - 1)
        elseif mtlb_logic(PitCtrl,"==",2) then
          PitCC(mtlb_s(mtlb_a(Vrated_ind,n),1),1) = mtlb_s(0,0.1*mtlb_s(index,1));  //PitSt - PitDel(index - 1)
        end;
      end;
    
      //Now make sure the pitch angles are monotonically increasing or decreasing
      fixPitCC = PitCC;
      badcase = [];
      for n = mtlb_imp(2,NumCases)
        if (abs(fixPitCC(n))-abs(fixPitCC(n-1)))<0 then //works for both the pitch to feather and pitch to stall cases
         fixPitCC = mtlb_i(fixPitCC,n,fixPitCC(n-1));
         //mark the index location
         badcase = [badcase;n];
        end;
      end;
      if isempty(badcase) then //no interpolation needed
      else //interpolate the bad cases
       fixV = V;
       fixV = mtlb_i(fixV,badcase,[]);
       fixPitCC = mtlb_i(fixPitCC,badcase,[]);
       PitCC = interp1(fixV,fixPitCC,V,"pchip","extrap");
      end;
    
      CC(:,3) = PitCC;  //Add the new pitch angles into the Combined Case
    else
      CC(:,3) = 0;  //Prated was never exceeded so there is no need to pitch the blades
    end;
  
    //Now build the final Combined Case Analysis WT_Perf input file to calculate the power curve now including the optimal pitch angles
    // !! L.410: Unknown function Build_Input not converted, original calling sequence used.
    Build_Input(2,RElm,Twist,Chord,Thickness,SpdSt,SpdEnd,SpdDel,0,0,CC)
    //Edit the input file to output the BED file
    fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".wtp","r+t");
    for i = mtlb_imp(1,mtlb_a(mtlb_a(mtlb_a(mtlb_a(31,NumSeg),6),NumSeg),3)) %v1_1 = mgetl(fid,1); if meof()~=0 then %v1_1 = -1;end; %v1_1;end;  //move the cursor through the file to the line about writing the BED file
    mseek(0,fid,"cur");  // L.414: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid,"True ");  //Changes the paramter to write the BED file to True
    mclose(fid);
    cd(RootDir+"\Output_Files\"+filename_main);  //Change the working directory to where the WT_Perf executable is located
    // !! L.417: Matlab function evalc not yet converted, original calling sequence used.
    evalc("!Wt_Perf_Jan1 "+filename_main+".wtp");
    cd(RootDir);  //Change the working directory back to the root directory
  
  
    //Need to calculate the AEP for the Variable-Pitch case, since it was not calculated previously
    //Reads in the Power vs. Velocity curve from WT_Perf output file
    fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\"+filename_main+".oup","rt");
    // !! L.424: Unknown function Read_WTP_Output not converted, original calling sequence used.
    PvsV = Read_WTP_Output(fid,1,1);  //This is a Combined Case Analysis output file
    mclose(fid);
  
  end;
  //===================== end of Variable Pitch cases ===================%


  //Send data to Post_Process, which writes the HARP_Opt output files
  if mtlb_logic(StructuralOpt,"==",0) then //use dummy values
   Output = [];
  end;

  // !! L.436: Unknown function Post_Process not converted, original calling sequence used.
  Post_Process(x,RElm,Twist,Chord,Thickness,NumCases,CC,pUvars,NumSheet,Output);

end;
//End of Function
endfunction
