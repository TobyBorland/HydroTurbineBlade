function [] = Interp_Thickness_Coefs(AFfile1,AFfile2,ThickVal1,ThickVal2,delT)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This function requires that the two airfoil files have the same number of
//tables for Reynolds numbers, the Reynolds number values do not have to be
//the same, just the same number of tables.  The Reynolds number of airfoil
//3 will be the weighted average of airfoils 1 and 2.

global("RootDir","filename_main","family")

//========================== TUNING PARAMETERS ============================%
stepsize = 0.001;//largest possible stepsize for AoA
CutOff = 1.5;//region 2.5 starts at CutOff*Stall_AoA3
Blend_L_r1 = 0.5;//start blending """"forward solution"""" into the """"backward solution"""" at 100*(Blend_L_r1)% of the way from AoA = 0 to Stall_AoA3.
Blend_L_r25 = 2;//Rate at which Lift blends into the AF_Prep solution from AoA = CutOff3 to 45. """"weight towards AF_Prep solution"""" = x^(1/Blend_L_r25),0<=x<=1.
Blend_D_r2 = 1;//Rate at which Drag blends into the AF_Prep solution from AoA = Stall_AoA3 to 45. """"weight towards AF_Prep solution"""" = x^(1/Blend_D_r2),0<=x<=1.
//=========================================================================%

//initialize the progress bar
// !! L.18: Matlab function waitbar not yet converted, original calling sequence used.
// L.18: (Warning name conflict: function name changed from waitbar to %waitbar).
H = %waitbar(0,"","Name","Creating interpolated airfoil files...");
// !! L.19: Matlab function get not yet converted, original calling sequence used.
// L.19: (Warning name conflict: function name changed from get to %get).
hchild = %get(H,"children");
// !! L.20: Matlab function get not yet converted, original calling sequence used.
// L.20: (Warning name conflict: function name changed from get to %get).
htitle = %get(hchild,"title");
// !! L.21: Matlab function set not yet converted, original calling sequence used.
// L.21: (Warning name conflict: function name changed from set to %set).
%set(htitle,"Interpreter","None","FontSize",8);

path1 = RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\"+AFfile1;
path2 = RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\"+AFfile2;

ThickVals_i = mtlb_imp(mtlb_s(ThickVal1,delT),-delT,mtlb_a(ThickVal2,delT))';
if isempty(ThickVals_i) then
  return;  //Stop this function early, there is no need to interpolate
end;
w = mtlb_s(ThickVals_i,ThickVal1)/mtlb_s(ThickVal2,ThickVal1);//weight values to use in averaging the coeficients

fid1 = mtlb_fopen(path1,"rt");
// !! L.33: Matlab function textscan not yet converted, original calling sequence used.
NumTables = cell2mat(textscan(fid1,"%f","HeaderLines",3));//number of tables for Reynolds number in the file
mclose(fid1);

fid1 = mtlb_fopen(path1,"rt");
fid2 = mtlb_fopen(path2,"rt");

for k = 1:max(size(ThickVals_i))

  if ThickVals_i(k)>=10 then
    Thick_Suffix = "_0"+msprintf("%3.0f",10*ThickVals_i(k));
  else
    Thick_Suffix = "_00"+msprintf("%3.0f",10*ThickVals_i(k));
  end;
  fid3 = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\"+family+Thick_Suffix+".dat","w+t");
  //Write this header to every new file
  // L.48: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid3,"AeroDyn airfoil file.  Compatible with AeroDyn v13.0 .\n");
  // L.49: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid3,"This file was generated automatically by HARP_Opt.\n");
  // L.50: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid3,"This data represents a %g percent thick airfoil, and was interpolated between the airfoil data in the files %s and %s.\n",ThickVals_i(k),AFfile1,AFfile2);
  // L.51: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid3,"  %g           Number of airfoil tables in this file\n",NumTables);

  mseek(0,fid1);  mseek(0,fid2);  //move the cursor back to the beginning of the files
  for jj = 1:4 %v0 = mgetl(fid1,1); if meof()~=0 then %v0 = -1;end; %v0; %v1 = mgetl(fid2,1); if meof()~=0 then %v1 = -1;end; %v1;end;  //move the cursor to the Reynolds number line of airfoil files 1 and 2

  for n = mtlb_imp(1,NumTables)
  
    //True airfoil 1
    // !! L.59: Matlab function textscan not yet converted, original calling sequence used.
    Re1 = cell2mat(textscan(fid1,"%f"));  %v2 = mgetl(fid1,1);  if meof()~=0 then %v2 = -1;end;  %v2;
    // !! L.60: Matlab function textscan not yet converted, original calling sequence used.
    Stall_AoA1 = cell2mat(textscan(fid1,"%f"));  %v3 = mgetl(fid1,1);  if meof()~=0 then %v3 = -1;end;  %v3;
    for jj = 1:6 %v4 = mgetl(fid1,1); if meof()~=0 then %v4 = -1;end; %v4;end;  //move the cursor to the first line of coefficients
    // !! L.62: Matlab function textscan not yet converted, original calling sequence used.
    AF1 = cell2mat(textscan(fid1,"%f %f %f %f %f","CollectOutput",1));  //reads in the whole table of coefficients
    // ! L.63: abs(isnan(AF1(:,5))) may be replaced by:
    // !    --> isnan(AF1(:,5)) if isnan(AF1(:,5)) is Real.
  
    if mtlb_any(abs(isnan(AF1(:,5)))) then AF1(:,5) = [];end;  //Delete 5th column if the user did not input this coefficient
    // ! L.64: abs(isnan(AF1(:,4))) may be replaced by:
    // !    --> isnan(AF1(:,4)) if isnan(AF1(:,4)) is Real.
  
    if mtlb_any(abs(isnan(AF1(:,4)))) then AF1(:,4) = [];end;  //Delete 4th column if the user did not input this coefficient
    a1 = AF1(:,1);  //Angle of Attack (deg)
    L1 = AF1(:,2);  //Lift polar
    D1 = AF1(:,3);  //Drag polar
    CLmax1 = mtlb_e(L1,FindInd(a1,Stall_AoA1));  //CLmax value
    CL_Alfa01 = interp1(a1,L1,0);  //CL at AoA=0
    CD_Alfa01 = interp1(a1,D1,0);  //CD at A0A = 0
    CutOff1 = CutOff*Stall_AoA1;
  
    //True airfoil 2
    // !! L.74: Matlab function textscan not yet converted, original calling sequence used.
    Re2 = cell2mat(textscan(fid2,"%f"));  %v5 = mgetl(fid2,1);  if meof()~=0 then %v5 = -1;end;  %v5;
    // !! L.75: Matlab function textscan not yet converted, original calling sequence used.
    Stall_AoA2 = cell2mat(textscan(fid2,"%f"));  %v6 = mgetl(fid2,1);  if meof()~=0 then %v6 = -1;end;  %v6;
    for jj = 1:6 %v7 = mgetl(fid2,1); if meof()~=0 then %v7 = -1;end; %v7;end;  //move the cursor to the first line of coefficients
    // !! L.77: Matlab function textscan not yet converted, original calling sequence used.
    AF2 = cell2mat(textscan(fid2,"%f %f %f %f %f","CollectOutput",1));  //reads in the whole table of coefficients
    // ! L.78: abs(isnan(AF2(:,5))) may be replaced by:
    // !    --> isnan(AF2(:,5)) if isnan(AF2(:,5)) is Real.
  
    if mtlb_any(abs(isnan(AF2(:,5)))) then AF2(:,5) = [];end;  //Delete 5th column if the user did not input this coefficient
    // ! L.79: abs(isnan(AF2(:,4))) may be replaced by:
    // !    --> isnan(AF2(:,4)) if isnan(AF2(:,4)) is Real.
  
    if mtlb_any(abs(isnan(AF2(:,4)))) then AF2(:,4) = [];end;  //Delete 4th column if the user did not input this coefficient
    a2 = AF2(:,1);  //Angle of Attack (deg)
    L2 = AF2(:,2);  //Lift polar
    D2 = AF2(:,3);  //Drag polar
    CLmax2 = mtlb_e(L2,FindInd(a2,Stall_AoA2));  //CLmax value
    CL_Alfa02 = interp1(a2,L2,0);  //CL at AoA=0
    CD_Alfa02 = interp1(a2,D2,0);  //CD at A0A = 0
    CutOff2 = CutOff*Stall_AoA2;
  
    //Interpolated airfoil 3, 3 is between 1 and 2
    Re3 = mtlb_a((1-w(k))*Re1,w(k)*Re2);
    Stall_AoA3 = mtlb_a((1-w(k))*Stall_AoA1,w(k)*Stall_AoA2);
    CLmax3 = mtlb_a((1-w(k))*CLmax1,w(k)*CLmax2);
    CL_Alfa03 = mtlb_a((1-w(k))*CL_Alfa01,w(k)*CL_Alfa02);
    CutOff3 = CutOff*Stall_AoA3;
    CD_Alfa03 = mtlb_a((1-w(k))*CD_Alfa01,w(k)*CD_Alfa02);
  
    //Write header for interpolated airfoil file
    // L.97: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"%6.3f         Table ID parameter (Reynolds number in millions)\n",Re3);
    // L.98: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"%6.3f         Stall angle (deg)\n",Stall_AoA3);
    // L.99: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"  0           Zero lift angle of attack (deg)\n");
    // L.100: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"  0           Cn slope for zero lift (dimensionless)\n");
    // L.101: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"  0           Cn at stall value for positive angle of attack\n");
    // L.102: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"  0           Cn at stall value for negative angle of attack\n");
    // L.103: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"  0           Angle of attack for minimum CD (deg)\n");
    // L.104: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"  0           Minimum CD value\n");
  
    if mtlb_logic(ThickVal1,"==",100) then //we are blending between a cirlce and airfoil, use linear interpolation
     //AirfoilPrep linear interpolation algorithm
     //Entire range, AoA = -180 to 180, using AoA values from original data in airfoil 2
     AoA3 = a2; //angles of attack from airfoil 2
     Lift3 = mtlb_a(w(k) .*L2,(1-w(k)) .*interp1(a1,L1,AoA3));
     Drag3 = mtlb_a(w(k) .*D2,(1-w(k)) .*interp1(a1,D1,AoA3));
    
     if size(AF1,2)==5 then //4 coefficients are present
      //Use AirfoilPrep algorithm to interpolate the 3rd and 4th coefficients (moment AND pressure coefficients)
      C3 = mtlb_a(w(k) .*AF2(:,4),(1-w(k)) .*interp1(AF1(:,1),AF1(:,4),AoA3));
      C4 = mtlb_a(w(k) .*AF2(:,5),(1-w(k)) .*interp1(AF1(:,1),AF1(:,5),AoA3));
      AF3 = [AoA3,Lift3,Drag3,C3,C4];
      // L.118: No simple equivalent, so mtlb_fprintf() is called.
      mtlb_fprintf(fid3,"%5.3f    %6.4f    %6.4f    %6.4f    %6.4f\n",AF3');
     elseif size(AF1,2)==4 then //only 3 coefficients are present
      //Use AirfoilPrep algorithm to interpolate the 3rd and 4th coefficients (moment OR pressure coefficients)
      C3 = mtlb_a(w(k) .*AF2(:,4),(1-w(k)) .*interp1(AF1(:,1),AF1(:,4),AoA3));
      AF3 = [AoA3,Lift3,Drag3,C3];
      // L.123: No simple equivalent, so mtlb_fprintf() is called.
      mtlb_fprintf(fid3,"%5.3f    %6.4f    %6.4f    %6.4f\n",AF3');
     else //only 2 coefficients (lift and drag) are present
      AF3 = [AoA3,Lift3,Drag3];
      // L.126: No simple equivalent, so mtlb_fprintf() is called.
      mtlb_fprintf(fid3,"%5.3f    %6.4f    %6.4f\n",AF3');
     end;
    
    else //We are blending between two airfoils, use the special interpolation algorithm
    
     //Divide into the 4 regions of angle-of-attack:
     //region 1: AoA = 0 to Stall AoA  (uses special interpolation)
     //region 2: AoA = Stall AoA to 45 (uses special interpolation)
     //region 3: AoA = 45 to 180       (uses linear interpolation)
     //region 4: AoA = -180 to 0       (uses linear interpolation)   
     //Regions 1 and 2 need to be same number of points but can have different stepsizes
     //divide into N points that guarantee a spacing = stepsize or smaller
     %v8_2 = floor(mtlb_a([Stall_AoA1;mtlb_s(CutOff1,Stall_AoA1);Stall_AoA2;mtlb_s(CutOff2,Stall_AoA2)]/stepsize,1)); N = mtlb_max(%v8_2,firstnonsingleton(%v8_2));
    
     //Angle of Attack vectors
     a1_r1 = mtlb_t(mtlb_linspace(0,Stall_AoA1,N));
     a1_r2 = mtlb_t(mtlb_linspace(Stall_AoA1,45,N));
     a1_r12 = [a1_r1;mtlb_e(a1_r2,2:$)]; //Caution: stepsize changes in this vector
     a2_r1 = mtlb_t(mtlb_linspace(0,Stall_AoA2,N));
     a2_r2 = mtlb_t(mtlb_linspace(Stall_AoA2,45,N));
     a2_r12 = [a2_r1;mtlb_e(a2_r2,2:$)]; //Caution: stepsize changes in this vector
     a3_r1 = linspace(0,Stall_AoA3,N)';
     a3_r2 = linspace(Stall_AoA3,45,N)';
     c = FindInd(a3_r2,CutOff3);
     a3_r25 = a3_r2(mtlb_imp(c,$));
    
     //Lift vectors
     L1_r1 = interp1(a1,L1,a1_r1);
     L1_r2 = interp1(a1,L1,a1_r2);
     L1_r12 = [L1_r1;L1_r2(2:$)];
     L2_r1 = interp1(a2,L2,a2_r1);
     L2_r2 = interp1(a2,L2,a2_r2);
     L2_r12 = [L2_r1;L2_r2(2:$)];
    
     //Drag vectors
     D1_r1 = interp1(a1,D1,a1_r1);
     D1_r2 = interp1(a1,D1,a1_r2);
     D1_r12 = [D1_r1;D1_r2(2:$)];
     D2_r1 = interp1(a2,D2,a2_r1);
     D2_r2 = interp1(a2,D2,a2_r2);
     D2_r12 = [D2_r1;D2_r2(2:$)];
    
     //take 1st derivatives of Regions 1 and 2
     //Lift Slopes
     dLda1_r1 = FiniteDiff(a1_r1,L1_r1);
     dLda1_r2 = FiniteDiff(a1_r2,L1_r2);
     dLda2_r1 = FiniteDiff(a2_r1,L2_r1);
     dLda2_r2 = FiniteDiff(a2_r2,L2_r2);
     dLda3_r1 = (1-w(k))*dLda1_r1+w(k)*dLda2_r1;
     dLda3_r2 = (1-w(k))*dLda1_r2+w(k)*dLda2_r2;
     //Drag Slopes
     dDda1_r1 = FiniteDiff(a1_r1,D1_r1);
     dDda1_r2 = FiniteDiff(a1_r2,D1_r2);
     dDda2_r1 = FiniteDiff(a2_r1,D2_r1);
     dDda2_r2 = FiniteDiff(a2_r2,D2_r2);
     dDda3_r1 = (1-w(k))*dDda1_r1+w(k)*dDda2_r1;
     dDda3_r2 = (1-w(k))*dDda1_r2+w(k)*dDda2_r2;
    
     //Calculate the Lift coefficients by solving an initial value problem
     //Region 1:
     //moving forwards, from AoA=0 to Stall AoA
     L3_r1_f = Euler(a3_r1,dLda3_r1,CL_Alfa03);
     //moving backwards, from Stall AoA to AoA=0
     %v9_2 = Euler(a3_r1($:-1:1,:),dLda3_r1($:-1:1,:),CLmax3); L3_r1_b = %v9_2($:-1:1,:);
     //now blend the forward & backward solutions w/ weighted averages
     W1 = (linspace(0,1,N)'-1) ./(1-Blend_L_r1)+1; W1 = mtlb_i(W1,mtlb_logic(W1,"<",0),0); W1 = mtlb_i(W1,isnan(W1),0);
     L3_r1 = (1-W1) .*L3_r1_f+W1 .*L3_r1_b;
     //Region 2:
     L3_r2 = Euler(a3_r2,dLda3_r2,CLmax3);
     //Region 2.5:
     L3_r25 = L3_r2(mtlb_imp(c,$));
    
     //Calculate the Drag coefficients by solving an initial value problem
     D3_r1 = Euler(a3_r1,dDda3_r1,CD_Alfa03);
     D3_r2 = Euler(a3_r2,dDda3_r2,D3_r1($));
    
     //     %AirfoilPrep linear interpolation algorithm
     //     %Entire range, AoA = -180 to 180, using AoA values from original data airfoil A
     //     L3_afp = (1-w(k)).*L1 + w(k).*interp1(a2,L2,a1);
     //     D3_afp = (1-w(k)).*D1 + w(k).*interp1(a2,D2,a1);
     //regions 1 and 2 combined
     L3_r12_afp = mtlb_a((1-w(k)) .*L1_r12,w(k) .*interp1(a2_r12,L2_r12,a1_r12));
     D3_r12_afp = mtlb_a((1-w(k)) .*D1_r12,w(k) .*interp1(a2_r12,D2_r12,a1_r12));
     //region 2
     D3_r2_afp = interp1(a1_r12,D3_r12_afp,a3_r2);
     //region 2.5
     L3_r25_afp = interp1(a1_r12,L3_r12_afp,a3_r25);
    
     //Blend Lift with the AirfoilPrep solution in Region 2.5
     W25 = (linspace(0,1,max(size(a3_r25)))'-1+1) .^(1/Blend_L_r25);
     L3_r25_blend = mtlb_a((1-W25) .*L3_r25,W25 .*L3_r25_afp);
     //Now combine regions 2 and 2.5 into a single vector
     L3_r2_blended = [L3_r2(mtlb_imp(1,c));L3_r25_blend(2:$)];
     //Blend Drag with the AirfoilPrep solution in Region 2
     W2 = (linspace(0,1,max(size(a3_r2)))'-1+1) .^(1/Blend_D_r2);
     W2 = mtlb_i(W2,mtlb_logic(W2,"<",0),0); W2 = mtlb_i(W2,isnan(W2),0);
     D3_r2_blended = mtlb_a((1-W2) .*D3_r2,W2 .*D3_r2_afp);
    
     //Calculate Lift & Drag in Regions 3 and 4 using the AirfoilPrep algorithm
     //region 3
     %v10_2 = mtlb_e(a1,mtlb_logic(a1,">",45)); a1_r3 = mtlb_e(a1,mtlb_imp(FindInd(a1,mtlb_min(%v10_2,firstnonsingleton(%v10_2))),$));
     L1_r3 = interp1(a1,L1,a1_r3);
     D1_r3 = interp1(a1,D1,a1_r3);
     L3_r3 = mtlb_a((1-w(k)) .*L1_r3,w(k) .*interp1(a2,L2,a1_r3));
     D3_r3 = mtlb_a((1-w(k)) .*D1_r3,w(k) .*interp1(a2,D2,a1_r3));
     //region 4
     %v11_2 = mtlb_e(a1,mtlb_logic(a1,"<",0)); a1_r4 = mtlb_e(a1,mtlb_imp(1,FindInd(a1,mtlb_max(%v11_2,firstnonsingleton(%v11_2)))));
     L1_r4 = interp1(a1,L1,a1_r4);
     D1_r4 = interp1(a1,D1,a1_r4);
     L3_r4 = mtlb_a((1-w(k)) .*L1_r4,w(k) .*interp1(a2,L2,a1_r4));
     D3_r4 = mtlb_a((1-w(k)) .*D1_r4,w(k) .*interp1(a2,D2,a1_r4));
    
     //Put Lift & Drag together into a single region again
     h1 = hcosspace(0,Stall_AoA3,50,1);
     h2 = hcosspace(Stall_AoA3,45,50,0);
     AoA3 = [a1_r4;h1;h2(2:$);a1_r3];
    
     Lift3 = interp1([a1_r4;a3_r1;a3_r2(2:$);a1_r3],[L3_r4;L3_r1;L3_r2_blended(2:$);L3_r3],AoA3);
    
     Drag3 = interp1([a1_r4;a3_r1;a3_r2(2:$);a1_r3],[D3_r4;D3_r1;D3_r2_blended(2:$);D3_r3],AoA3);
    
     if size(AF1,2)==5 then //4 coefficients are present
      //Use AirfoilPrep algorithm to interpolate the 3rd and 4th coefficients (moment AND pressure coefficients)
      C3_afp = mtlb_a((1-w(k)) .*AF1(:,4),w(k) .*interp1(AF2(:,1),AF2(:,4),AF1(:,1)));
      C3 = interp1(AF1(:,1),C3_afp,AoA3);
      C4_afp = mtlb_a((1-w(k)) .*AF1(:,5),w(k) .*interp1(AF2(:,1),AF2(:,5),AF1(:,1)));
      C4 = interp1(AF1(:,1),C4_afp,AoA3);
      AF3 = [AoA3,Lift3,Drag3,C3,C4];
      // L.254: No simple equivalent, so mtlb_fprintf() is called.
      mtlb_fprintf(fid3,"%5.3f    %6.4f    %6.4f    %6.4f    %6.4f\n",AF3');
     elseif size(AF1,2)==4 then //only 3 coefficients are present
      //Use AirfoilPrep algorithm to interpolate the 3rd and 4th coefficients (moment OR pressure coefficients)
      C3_afp = mtlb_a((1-w(k)) .*AF1(:,4),w(k) .*interp1(AF2(:,1),AF2(:,4),AF1(:,1)));
      C3 = interp1(AF1(:,1),C3_afp,AoA3);
      AF3 = [AoA3,Lift3,Drag3,C3];
      // L.260: No simple equivalent, so mtlb_fprintf() is called.
      mtlb_fprintf(fid3,"%5.3f    %6.4f    %6.4f    %6.4f\n",AF3');
     else //only 2 coefficients (lift and drag) are present
      AF3 = [AoA3,Lift3,Drag3];
      // L.263: No simple equivalent, so mtlb_fprintf() is called.
      mtlb_fprintf(fid3,"%5.3f    %6.4f    %6.4f\n",AF3');
     end;
    
    end;
  
  
    // L.269: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3,"EOT %g\n",n);
    %v12 = mgetl(fid1,1);  if meof()~=0 then %v12 = -1;end;  %v12;  //Skip through the line that reads ''EOT''
    %v13 = mgetl(fid2,1);  if meof()~=0 then %v13 = -1;end;  %v13;  //Skip through the line that reads ''EOT''
  end;

  mclose(fid3);

  // Report current progress in the waitbar''s message field
  // !! L.277: Matlab function waitbar not yet converted, original calling sequence used.
  // L.277: (Warning name conflict: function name changed from waitbar to %waitbar).
  %waitbar(k/max(size(ThickVals_i)),H,"Interpolating between "+AFfile1+" and "+AFfile2)


end;
mtlb_delete(H);// DELETE the waitbar; don''t try to CLOSE it.

//Now close the files
mclose(fid1);
mclose(fid2);


//============================= PLOTS =====================================%
// figure(''Name'',''Lift: Comparison to AirfoilPrep'',''color'',[1 1 1]);
// hold on;
// plot(a1,L1,''o-b'',a2,L2,''o-r'',...
//     [a1_r4;a3_r1;a3_r2(2:end);a1_r3],[L3_r4;L3_r1;L3_r2(2:end);L3_r3],''-k'',...
//      AoA3,Lift3,'':k'',a1,L3_afp,''-g'')
// legend(''1'',''2'',''3 special'',''3 special w/ blending'',''3 AF_Prep'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''Lift Coefficient'',''FontSize'',32);
// axis([-180 180 -0.6 2.0]);
// hold off;
// 
// figure(''Name'',''Drag: Comparison to AirfoilPrep'',''color'',[1 1 1]);
// hold on;
// plot(a1,D1,''o-b'',a2,D2,''o-r'',...
//     [a1_r4;a3_r1;a3_r2(2:end);a1_r3],[D3_r4;D3_r1;D3_r2(2:end);D3_r3],''-k'',...
//      AoA3,Drag3,'':k'',a1,D3_afp,''-g'')
// legend(''1'',''2'',''3 special'',''3 special w/ blending'',''3 AF_Prep'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''Drag Coefficient'',''FontSize'',32);
// axis([0 45 0 0.7]);
// hold off;
// 
// if size(AF1,2)>3; %3rd coefficient is present
// figure(''Name'',''3rd Coefficient'',''color'',[1 1 1]);
// hold on;
// plot(a1,AF1(:,4),''o-b'',a2,AF2(:,4),''o-r'',...
//      a1,C3_afp,''-g'')
// legend(''1'',''2'',''3 AF_Prep'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''3rd Coefficient'',''FontSize'',32);
// hold off;
// end
// 
// if size(AF1,2)>4; %4th coefficient is present
// figure(''Name'',''4th Coefficient'',''color'',[1 1 1]);
// hold on;
// plot(a1,AF1(:,5),''o-b'',a2,AF2(:,5),''o-r'',...
//      a1,C4_afp,''-g'')
// legend(''1'',''2'',''3 AF_Prep'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''4th Coefficient'',''FontSize'',32);
// hold off;
// end
// 
// figure(''Name'',''Lift: comparison to original data'',''color'',[1 1 1]);
// hold on;
// plot(a1,L1,''o-b'',a2,L2,''o-r'',...
//      AoA3,Lift3,''x-k'')
// legend(''1: original'',''2: original'',''3: final interpolated'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''Lift Coefficient'',''FontSize'',32);
// axis([-180 180 -0.6 1.5]);
// hold off;
// 
// figure(''Name'',''Drag: comparison to original data'',''color'',[1 1 1]);
// hold on;
// plot(a1,D1,''o-b'',a2,D2,''o-r'',...
//      AoA3,Drag3,''x-k'')
// legend(''1: original'',''2: original'',''3: final interpolated'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''Drag Coefficient'',''FontSize'',32);
// axis([0 45 0 0.7]);
// hold off;
// 
// if size(AF1,2)>3; %3rd coefficient is present
// figure(''Name'',''3rd Coefficient:  comparison to original data'',''color'',[1 1 1]);
// hold on;
// plot(a1,AF1(:,4),''o-b'',a2,AF2(:,4),''o-r'',...
//      AoA3,C3,''x-k'')
// legend(''1'',''2'',''3 AF_Prep'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''3rd Coefficient'',''FontSize'',32);
// hold off;
// end
// 
// if size(AF1,2)>4; %4th coefficient is present
// figure(''Name'',''4th Coefficient:  comparison to original data'',''color'',[1 1 1]);
// hold on;
// plot(a1,AF1(:,5),''o-b'',a2,AF2(:,5),''o-r'',...
//      AoA,C4,''x-k'')
// legend(''1'',''2'',''3 AF_Prep'',2);
// xlabel(''Angle of Attack (deg)'',''FontSize'',32);
// ylabel(''4th Coefficient'',''FontSize'',32);
// hold off;
// end
endfunction
