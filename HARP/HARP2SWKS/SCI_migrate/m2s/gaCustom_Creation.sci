function [Population] = gaCustom_Creation(GenomeLength,unused1,options)

// Output variables initialisation (not found in input variables)
Population=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This function creates the initial population for the genetic algorithm.
//Basically, this file needs to create a random distribution of the genetic
//algorithm variables (control points for the Bezier curve fits), but these 
//variables MUST obey the upper and lower bounds, and linear inequality 
//constraints defined in HARP_Opt.m
// 
//Inputs:    GenomeLength: the number of variables that the genetic algorithm is using. 16 if fixed speed and 15 if variable speed turbine.
//           unused1: apparently an unused argument defined by MATLAB
//           options: this is the options structure created for the genetic algorithm, this options struction is defined in HARP_Opt.m
// 
//Outputs:   Population: a matrix, containing all of the 15 or 16 dimensional vectors of the initial population


global("SpdCtrl","RotorRad","ThickMethod","Thickness_values","NumSeg","TwistLB","TwistUB","ChordLB","ChordUB","ThickUB","ThickLB","OmgMin","OmgMax");

//=========================================================================%
//Temporary variables for debugging
// LB = [0 0 0 0 0 0 0 0 0 0 0.2 0.4 0.8 20];
// UB = [10 10 10 10 10 10 10 10 10 10 0.4 0.6 0.9 40];
//=========================================================================%
totalPopulation = mtlb_sum(mtlb_e(options,"PopulationSize"));
//initPopProvided = size(options.InitialPopulation,1);
initPopProvided = 0;
//individualsToCreate = totalPopulation - initPopProvided
individualsToCreate = totalPopulation;

// Initialize Population to be created
Population = zeros(totalPopulation,GenomeLength);
// Use initial population provided already
if initPopProvided>0 then
  Population(1:initPopProvided,:) = mtlb_e(options,"InitialPopulation");
end;

//Preallocate
Twist = zeros(1,5);
Chord = zeros(1,5);
if mtlb_logic(ThickMethod,"==",1) then
  Thick = zeros(1,max(size(Thickness_values))-1);
elseif mtlb_logic(ThickMethod,"==",2) then
  Thick = zeros(1,max(size(Thickness_values)));
end;
X = zeros(1,GenomeLength);


// L.46: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(1,"\nCreating initial population of feasible individuals ... if the\n");
// L.47: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(1,"number of feasible individuals found does not increase after a\n");
// L.48: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(1,"short period of time, you may have set unfeasible input values\n");
// L.49: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(1,"and will need to cancel (CTRL-C) HARP_Opt and input new values:\n\n");

// % % % % fprintf(1,''Feasible Individuals   Total Individuals\n'');
// % % % % fprintf(1,''  Found (need %g)         Created\n'',individualsToCreate);
// % % % % fprintf(1,''------------------------------------------\n'');
// % % % % fprintf(1,''                               '');


// !! L.57: Matlab function waitbar not yet converted, original calling sequence used.
// L.57: (Warning name conflict: function name changed from waitbar to %waitbar).
H = %waitbar(0,"","Name","Creating initial population...","CreateCancelBtn","setappdata(gcbf,''canceling'',1)");
// !! L.58: Matlab function setappdata not yet converted, original calling sequence used.
setappdata(H,"canceling",0)

C = 1;
D = 1;
while mtlb_logic(D,"<=",individualsToCreate)

  //Twist and Chord
  Twist(1,1) = mtlb_a(rand()*mtlb_s(mtlb_e(TwistUB,1),mtlb_e(TwistLB,1)),mtlb_e(TwistLB,1));
  Chord(1,1) = mtlb_a(rand()*mtlb_s(mtlb_e(ChordUB,1),mtlb_e(ChordLB,1)),mtlb_e(ChordLB,1));
  for k = 2:5
    if mtlb_logic(Twist(1,k-1),">",mtlb_e(TwistUB,k)) then
      ubT = mtlb_e(TwistUB,k);
    else
      ubT = Twist(1,k-1);
    end;
    if mtlb_logic(Chord(1,k-1),">",mtlb_e(ChordUB,k)) then
      ubC = mtlb_e(ChordUB,k);
    else
      ubC = Chord(1,k-1);
    end;
    lbT = mtlb_e(TwistLB,k);
    lbC = mtlb_e(ChordLB,k);
    Twist(1,k) = mtlb_a(rand()*mtlb_s(ubT,lbT),lbT);
    Chord(1,k) = mtlb_a(rand()*mtlb_s(ubC,lbC),lbC);
  end;

  //Percent Thickness radial stations
  //if bounds are entered for the Percent Thickness
  if max(size(Thickness_values))>1 then
  
    %v04 = %f;  if isempty(ThickLB) then %v04 = isempty(ThickUB);end;
    if %v04 then
      ThickLB = mtlb_i(ThickLB,1:max(size(Thickness_values)),%nan);
      ThickUB = ThickLB;
    end;
  
    if mtlb_logic(ThickMethod,"==",1) then
      NumAFvars = max(size(Thickness_values))-1;
    elseif mtlb_logic(ThickMethod,"==",2) then
      NumAFvars = max(size(Thickness_values));
    end;
  
    Thick(1,1:NumAFvars) = %nan;  //preallocating
    //first take care of the ones with BOTH lower AND upper bounds
    for n = 1:NumAFvars
      %v14 = %f;  if ~isnan(mtlb_e(ThickLB,n)) then %v14 = ~isnan(mtlb_e(ThickUB,n));end;
      if %v14 then
        Thick = mtlb_i(Thick,n,mtlb_a(rand()*mtlb_s(mtlb_e(ThickUB,n),mtlb_e(ThickLB,n)),mtlb_e(ThickLB,n)));
      end;
    end;
  
    for n = NumAFvars:-1:1
    
      // only upper bound exists
      %v24 = %f;  if isnan(mtlb_e(ThickLB,n)) then %v24 = abs(mtlb_e(ThickUB,n))<%inf;end;
      if %v24 then
      
        randUB(1,n) = matrix(mtlb_e(ThickUB,n),1,-1);
      
        //find a lower bound
        for N = n:-1:1
          %v4_3 = %f;  if isnan(mtlb_e(ThickLB,N)) then %v4_3 = abs(mtlb_e(ThickUB,N))<%inf;end;
          if %v35 then
            randLB(1,n) = matrix(mtlb_e(ThickLB,N),1,-1);
            break
          elseif abs(Thick(N))<%inf then
            randLB(1,n) = Thick(N);
            break
          elseif %v4_3 & N~=n then
            randLB(1,n) = matrix(mtlb_e(ThickUB,N),1,-1);
            break
          elseif N==1 then
            randLB(1,n) = 0;
          else
            continue
          end;
        end;
      end;
    
      // only lower bound exists
      %v54 = %f;  if abs(mtlb_e(ThickLB,n))<%inf then %v54 = isnan(mtlb_e(ThickUB,n));end;
      if %v54 then
      
        randLB = mtlb_i(randLB,n,mtlb_e(ThickLB,n));
      
        //find an upper bound
        for N = n:NumAFvars
          %v65 = %f;  if abs(mtlb_e(ThickUB,N))<%inf then %v65 = isnan(Thick(N));end;
          if %v65 then
            randUB = mtlb_i(randUB,n,mtlb_e(ThickUB,N));
            break
          elseif abs(Thick(N))<%inf then
            randUB = mtlb_i(randUB,n,Thick(N));
            break
          elseif N==NumAFvars then
            randUB = mtlb_i(randUB,n,1);
          else
            continue
          end;
        end;
      end;
    
      // no bounds exist
      %v74 = %f;  if isnan(mtlb_e(ThickLB,n)) then %v74 = isnan(mtlb_e(ThickUB,n));end;
      if %v74 then
        //find a lower bound
        for N = n:-1:1
          %v9_3 = %f;  if isnan(mtlb_e(ThickLB,N)) then %v9_3 = abs(mtlb_e(ThickUB,N))<%inf;end;
          if %v85 then
            randLB = mtlb_i(randLB,n,mtlb_e(ThickLB,N));
            break
          elseif abs(Thick(N))<%inf then
            randLB = mtlb_i(randLB,n,Thick(N));
            break
          elseif %v9_3 & N~=n then
            randLB = mtlb_i(randLB,n,mtlb_e(ThickUB,N));
            break
          elseif N==1 then
            randLB = mtlb_i(randLB,n,0);
          else
            continue
          end;
        end;
      
        //find an upper bound
        for N = n:NumAFvars
          %v105 = %f;  if abs(mtlb_e(ThickUB,N))<%inf then %v105 = isnan(Thick(N));end;
          if %v105 then
            randUB = mtlb_i(randUB,n,mtlb_e(ThickUB,N));
            break
          elseif abs(Thick(N))<%inf then
            randUB = mtlb_i(randUB,n,Thick(N));
            break
          elseif N==NumAFvars then
            randUB = mtlb_i(randUB,n,1);
          else
            continue
          end;
        end;
      end;
    
      // both bounds exist
      %v114 = %f;  if abs(mtlb_e(ThickLB,n))<%inf then %v114 = abs(mtlb_e(ThickUB,n))<%inf;end;
      if %v114 then
        randLB = mtlb_i(randLB,n,mtlb_e(ThickLB,n));
        randUB = mtlb_i(randUB,n,mtlb_e(ThickUB,n));
      end;
    
      if isnan(Thick(n)) then
        Thick = mtlb_i(Thick,n,mtlb_a(rand()*(randUB(n)-randLB(n)),randLB(n)));
      end;
    end;
  
    X(1,1:10+NumAFvars) = [Twist(1,:),Chord(1,:),Thick(1,:)];
  else
    X(1,1:10) = [Twist(1,:),Chord(1,:)];
  end;

  if mtlb_logic(SpdCtrl,"==",0) then
    X(1,GenomeLength) = mtlb_a(rand()*mtlb_s(OmgMax,OmgMin),OmgMin);
  end;

  //evaluate the randomly created individual to see if it produces any kind of error
  F = Main(X);
  if F(1)~=%inf then
    Population(D,:) = X;
    D = D+1;
  end;

  // % % % %     fprintf(1,''\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'');
  // % % % %     fprintf(1,''       %-10.0f\t\t\t%-10.0f'',D-1,C);

  // Check for Cancel button press
  // !! L.221: Matlab function getappdata not yet converted, original calling sequence used.

  if getappdata(H,"canceling") then
    mtlb_delete(H);  // DELETE the waitbar; don''t try to CLOSE it.
    mclose("all");
    error("Program Terminated: User aborted program during creation of initial population.");
  end;
  // Report current estimate in the waitbar''s message field
  // !! L.227: Matlab function sprintf not yet converted, original calling sequence used.
  // L.227: (Warning name conflict: function name changed from sprintf to %sprintf).
  // !! L.227: Matlab function waitbar not yet converted, original calling sequence used.
  // L.227: (Warning name conflict: function name changed from waitbar to %waitbar).
  %waitbar((D-1)/individualsToCreate,H,%sprintf("Feasible Found = %g (Need %g), Total Created = %g",D-1,individualsToCreate,C))

  C = C+1;  //update counter for next iteration
end;
mtlb_delete(H);// DELETE the waitbar; don''t try to CLOSE it.
//fprintf(1,'' Done. \n'');

//=========================================================================%
//Plot the initial population
r = zeros(NumSeg,size(Population,1));
tw = zeros(NumSeg,size(Population,1));
c = zeros(NumSeg,size(Population,1));
pt = zeros(NumSeg,size(Population,1));
dt = zeros(NumSeg,size(Population,1));
for k = 1:size(Population,1)
  [ShapeError,RElm,TWIST,CHORD,PERCENT_THICKNESS,DIMENSIONAL_THICKNESS] = Define_Blade_Shape(Population(k,:));
  r(:,k) = RElm;
  tw(:,k) = TWIST;
  c(:,k) = CHORD;
  pt(:,k) = PERCENT_THICKNESS;
  dt(:,k) = DIMENSIONAL_THICKNESS;
end;
// !! L.249: Matlab function figure not yet converted, original calling sequence used.
Fig1 = figure(1);
// !! L.250: Matlab function set not yet converted, original calling sequence used.
// L.250: (Warning name conflict: function name changed from set to %set).
%set(Fig1,"color","white","name","Initial Population","numberTitle","off");

if mtlb_logic(SpdCtrl,"==",1) then
  subplot(4,1,1);  plot(r,tw,"-b","LineWidth",1,"LineSmoothing","off");  xlabel("Blade Radius (m)");  ylabel("Pre-twist (deg)");  // !! L.253: Matlab function xlim not yet converted, original calling sequence used.
  xlim([0,RotorRad]);  %v0_1 = gca();  %v0_1.box = "on";
  subplot(4,1,2);  plot(r,c,"-r","LineWidth",1,"LineSmoothing","off");  xlabel("Blade Radius (m)");  ylabel("Chord (m)");  // !! L.254: Matlab function xlim not yet converted, original calling sequence used.
  xlim([0,RotorRad]);  %v1_1 = gca();  %v1_1.box = "on";
  subplot(4,1,3);  plot(r,pt,"-","Color",[0,0.5,0],"LineWidth",1,"LineSmoothing","off");  xlabel("Blade Radius (m)");  ylabel("Thickness (%)");  // !! L.255: Matlab function xlim not yet converted, original calling sequence used.
  xlim([0,RotorRad]);  %v2_1 = gca();  %v2_1.box = "on";
  subplot(4,1,4);  plot(r,dt,"-k","LineWidth",1,"LineSmoothing","off");  xlabel("Blade Radius (m)");  ylabel("Thickness (m)");  // !! L.256: Matlab function xlim not yet converted, original calling sequence used.
  xlim([0,RotorRad]);  %v3_1 = gca();  %v3_1.box = "on";
elseif mtlb_logic(SpdCtrl,"==",0) then //Fixed Speed
 subplot(4,3,[1,2]); plot(r,tw,"-b","LineWidth",1,"LineSmoothing","off"); xlabel("Blade Radius (m)"); ylabel("Pre-twist (deg)"); // !! L.258: Matlab function xlim not yet converted, original calling sequence used.
 xlim([0,RotorRad]); %v4_2 = gca(); %v4_2.box = "on";
 subplot(4,3,[4,5]); plot(r,c,"-r","LineWidth",1,"LineSmoothing","off"); xlabel("Blade Radius (m)"); ylabel("Chord (m)"); // !! L.259: Matlab function xlim not yet converted, original calling sequence used.
 xlim([0,RotorRad]); %v5_2 = gca(); %v5_2.box = "on";
 subplot(4,3,[7,8]); plot(r,pt,"-","Color",[0,0.5,0],"LineWidth",1,"LineSmoothing","off"); xlabel("Blade Radius (m)"); ylabel("Thickness (%)"); // !! L.260: Matlab function xlim not yet converted, original calling sequence used.
 xlim([0,RotorRad]); %v6_2 = gca(); %v6_2.box = "on";
 subplot(4,3,[10,11]); plot(r,dt,"-k","LineWidth",1,"LineSmoothing","off"); xlabel("Blade Radius (m)"); ylabel("Thickness (m)"); // !! L.261: Matlab function xlim not yet converted, original calling sequence used.
 xlim([0,RotorRad]); %v7_2 = gca(); %v7_2.box = "on";
 subplot(4,3,[3,12]); plot(1:1:size(Population,1),Population(:,GenomeLength),"ok","LineSmoothing","off"); xlabel("Individual"); ylabel("Rotor Speed (rpm)"); // !! L.262: Matlab function xlim not yet converted, original calling sequence used.
 xlim([0,size(Population,1)]); %v8_2 = gca(); %v8_2.box = "on";
end;
endfunction
