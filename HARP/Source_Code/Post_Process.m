function Post_Process(x,RElm,Twist,Chord,Thickness,NumCases,CC,pUvars,NumSheet,Output)
%This function performs various tasks after the optimization code has found
%the best solution:
%       -Builds & Evaluates the final WT_Perf files which correspond to the
%        best solution found by the genetic algorithm
%       -Edits the headers of the WT_Perf files to show the time/date they
%        were created
%       -Prints the output data to the Excel output file

 
%Inputs:    x: vector, the best individual found by the genetic algorithm
%           RELm: vector, contains the radius values of the blade elements
%           Twist: vector, contains the Bezier curve Twist values of the blade elements, corresponding to the best solution found by the GA
%           Chord: vector, contains the Bezier curve Chord values of the blade elements, corresponding to the best solution found by the GA
%           ThicknessStep: vector, contains the step function % Thickness values of the blade elements, corresponding to the best solution found by the GA
%           NumCases: scalar, the number of flow speeds for the analysis
%           CC: matrix, contains the flow speeds, rotor speeds, and blade pitches for a Combined Case Analysis
%           AEP: scalar, the annual energy production value, corresponding to the best solution found by the GA
% Output variables from Structural_Fitness.m:
                        % Output.TotalMass = BladeMass;
                        % Output.MassUnitLen = MassDist./LocSpan;
                        % Output.ShellThickness = ShellThickness;
                        % Output.Iuu = Iuu_sh;
                        % Output.Ivv = Ivv_sh;
                        % Output.StrainLE = StrainLE;
                        % Output.StrainTE = StrainTE;
                        % Output.StrainUpper = StrainUpper;
                        % Output.StrainLower = StrainLower;
                        % Output.StressLE = StressLE;
                        % Output.StressTE = StressTE;
                        % Output.StressUpper = StressUpper;
                        % Output.StressLower = StressLower;
                        
%Outputs:   none, executing this file just performs the tasks listed above

global RootDir filename_main RotorRad NumSeg Prated SpdSt SpdEnd SpdDel SpdCtrl PitCtrl...
       PopSize ProbDist U_mean Weib_Umean Weib_k Weib_c user_pU_interp NumVars...
       StructuralOpt ParetoSet Pareto_filename RecordFailures;    

persistent ExcelData

V = (SpdSt:SpdDel:SpdEnd)';


%% Build the final WT_Perf Input_Files, then evaluate, then read in the output

%Fixed-Speed Fixed-Pitch case
if SpdCtrl == 0 && PitCtrl == 0;
    %For Passive Control Turbine (FSFP)
    %Edit the input .wtp file to output Power,Cp,Flap,Torque,Thrust, and BED
    fid = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '.wtp'],'r+t');
    for i = 1:(31 + NumSeg + 6 + NumSeg + 3); fgetl(fid); end;
    fseek(fid,0,'cof'); count = fprintf(fid,'True '); fseek(fid, -count, 'cof');
    for i = 1:10; fgetl(fid); end;
    fseek(fid,0,'cof'); count = fprintf(fid,'True '); fseek(fid, -count, 'cof');
    fgetl(fid);
    fseek(fid,0,'cof'); count = fprintf(fid,'True '); fseek(fid, -count, 'cof');
    fgetl(fid);
    fseek(fid,0,'cof'); count = fprintf(fid,'True '); fseek(fid, -count, 'cof');
    fgetl(fid); 
    fseek(fid,0,'cof'); count = fprintf(fid,'True '); fseek(fid, -count, 'cof');
    fgetl(fid);
    fseek(fid,0,'cof'); count = fprintf(fid,'True '); fseek(fid, -count, 'cof');
    fclose(fid);
    
    cd([RootDir '\Output_Files\' filename_main]);   %Change the working directory to where the WT_Perf executable is located
    evalc(['!Wt_Perf_Jan1 ' filename_main '.wtp']);
    cd(RootDir);                                    %Change the working directory back to the root directory
    
    %Need to define a dummy Combined Case, only to write to Excel
    OmgCC(1:NumCases,1) = x(end);
    PitCC(1:NumCases,1) = 0;
    CC = [V OmgCC PitCC];

else %for any other combinations of variable-speed or variable-pitch
    %Edit the input .wtp file to output the BED file
    fid = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '.wtp'],'r+t');
    for i = 1:(31 + NumSeg + 6 + NumSeg + 3); fgetl(fid); end;
    fseek(fid,0,'cof'); fprintf(fid,'True ');
    fclose(fid);
    
    cd([RootDir '\Output_Files\' filename_main]);   %Change the working directory to where the WT_Perf executable is located
    evalc(['!Wt_Perf_Jan1 ' filename_main '.wtp']);
    cd(RootDir);                                    %Change the working directory back to the root directory
end

%Rename the WT_Perf files to correspond to the Pateto frontier
if StructuralOpt == 1;
movefile([RootDir '\Output_Files\' filename_main '\' filename_main '.wtp'],...
         [RootDir '\Output_Files\' filename_main '\' Pareto_filename '.wtp']);
movefile([RootDir '\Output_Files\' filename_main '\' filename_main '.oup'],...
         [RootDir '\Output_Files\' filename_main '\' Pareto_filename '.oup']);
movefile([RootDir '\Output_Files\' filename_main '\' filename_main '.bed'],...
         [RootDir '\Output_Files\' filename_main '\' Pareto_filename '.bed']);
end

%Read-in the output data
fid = fopen([RootDir '\Output_Files\' filename_main '\' Pareto_filename '.oup'],'rt');
if SpdCtrl == 0 && PitCtrl == 0;%This was a Parametric Analysis output file
    Power  = textscan(fid,'%f %f','HeaderLines',10,'delimiter','\t','CollectOutput',0);
    Cp     = textscan(fid,'%f %f','HeaderLines',5,'delimiter','\t','CollectOutput',0);
    Torque = textscan(fid,'%f %f','HeaderLines',5,'delimiter','\t','CollectOutput',0);
    Flap   = textscan(fid,'%f %f','HeaderLines',5,'delimiter','\t','CollectOutput',0);
    Thrust = textscan(fid,'%f %f','HeaderLines',5,'delimiter','\t','CollectOutput',0);
    Power  = Power{2};
    Cp     = Cp{2};
    Torque = Torque{2};
    Flap   = Flap{2};
    Thrust = Thrust{2};
else %This was a Combined Case output file
    A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f','HeaderLines',8);
    Power  = A{5};
    Cp     = A{9};
    Torque = A{6};
    Flap   = A{8};
    Thrust = A{7};    
end
fclose(fid);


%% Edit the input & output file headers to show the program version and date
fid = fopen([RootDir '\Output_Files\' filename_main '\' Pareto_filename '.wtp'],'r+t');
date = datestr(now, 'mmmm dd, yyyy HH:MM AM');
fgetl(fid); 
fseek(fid,0,'cof');
fprintf(fid,'This file (Optimal Solution) was generated automatically by HARP_Opt on %s \n',date);
fprintf(fid,'This line is for user comments.                                   ');
fclose(fid);

fid_oup = fopen([RootDir '\Output_Files\' filename_main '\' Pareto_filename '.oup'],'r+t');
for n = 1:3; fgetl(fid_oup); end;
fseek(fid_oup,0,'cof');
fprintf(fid_oup,'This file (Optimal Solution) was generated automatically by HARP_Opt                                 ');
fclose(fid_oup);

fid_bed = fopen([RootDir '\Output_Files\' filename_main '\' Pareto_filename '.bed'],'r+t');
for n = 1:3; fgetl(fid_bed); end;
fseek(fid_bed,0,'cof');
fprintf(fid_bed,'This file (Optimal Solution) was generated automatically by HARP_Opt                                 ');
fclose(fid_bed);

%% Now Print the Output Data to an Excel File

if SpdCtrl == 0 && PitCtrl == 0
    CC(1:NumCases,2) = x(end);
    CC(1:NumCases,3) = 0;
end

%Calculate the Annual Energy Production (AEP)
    %First check if Power is negative, if so change Power value to 0
    indices = Power>0;
    PWR = Power.*indices;
AEP = get_AEP([V PWR],pUvars);

%Calculate the Capacity Factor
CF = 100*AEP./(Prated.*8760); %Capacity factor in percent

if ProbDist == 0;
    AEP = 'N/A';
    AEPstr = [', AEP was not calculated. Capacity Factor CF = '];
    FlowProb(1:NumCases,1) = 0;
elseif ProbDist == 1;
    AEPstr = [' kW-hr/yr, based on a Rayleigh distribution with a Umean = ' num2str(U_mean) 'm/s. Capacity Factor CF = '];
    FlowProb = (pi./2).*(V./U_mean.^2).*exp(-(pi./4).*(V./U_mean).^2);
elseif ProbDist == 2;
    AEPstr = [' kW-hr/yr, based on a Weibull distribution with a Umean = ' num2str(Weib_Umean) 'm/s, k = ' num2str(Weib_k) ', and c = ' num2str(Weib_c) '. Capacity Factor CF = '];
    FlowProb = (Weib_k./Weib_c).*((V./Weib_c).^(Weib_k-1)).*exp(-(V./Weib_c).^Weib_k);
else
    AEPstr = [' kW-hr/yr, based on a user defined flow probability. Capacity Factor CF = '];
    FlowProb = user_pU_interp;
end

if ~isempty(Output)
BladeMass = Output.TotalMass;
MassUnitLen = Output.MassUnitLen;
ShellThickness = Output.ShellThickness; %convert to (mm)
Iuu = Output.Iuu;
Ivv = Output.Ivv;
StrainLE = 1e6*Output.StrainLE; %convert to microstrain
StrainTE = 1e6*Output.StrainTE;
StrainUpper = 1e6*Output.StrainUpper;
StrainLower = 1e6*Output.StrainLower;
StressLE = Output.StressLE;
StressTE = Output.StressTE;
StressUpper = Output.StressUpper;
StressLower = Output.StressLower;
Mnorm = Output.Mnorm;
Mtan = Output.Mtan;
else
    %Need to give dummy values
BladeMass = 0;
MassUnitLen = zeros(NumSeg,1);
ShellThickness = zeros(NumSeg,1);
Iuu = zeros(NumSeg,1);
Ivv = zeros(NumSeg,1);
StrainLE = zeros(NumSeg,1);
StrainTE = zeros(NumSeg,1);
StrainUpper = zeros(NumSeg,1);
StrainLower = zeros(NumSeg,1);
StressLE = zeros(NumSeg,1);
StressTE = zeros(NumSeg,1);
StressUpper = zeros(NumSeg,1);
StressLower  = zeros(NumSeg,1);
Mnorm = zeros(NumSeg,1);
Mtan = zeros(NumSeg,1);
end

M1 = [RElm./RotorRad RElm Twist Chord Thickness Chord.*Thickness./100 ...
      ShellThickness MassUnitLen Iuu Ivv StrainUpper StrainLower StrainLE StrainTE StressUpper StressLower StressLE StressTE Mnorm Mtan];
M2 = [CC(:,1) CC(:,2) CC(:,3) Power Cp Flap Torque Thrust FlowProb];
xls_filename = [RootDir '\Output_Files\' filename_main '\' filename_main '_Output.xls'];

HeaderLines = 6;
ExcelData{NumSeg+HeaderLines,50,NumSheet} = [];
ExcelData{1,23,NumSheet} = ['All data & figures correspond to the filename "' Pareto_filename '"'];
ExcelData{2,23,NumSheet} = ['AEP =  ' num2str(AEP,'%3.0f') AEPstr num2str(CF,'%4.1f') '%'];
ExcelData{3,23,NumSheet} = ['Total Blade Mass =  ' num2str(BladeMass,'%2.2f') ' kg.'];
ExcelData{4,23,NumSheet} = 'Individual x=';
ExcelData(4,24:(24+NumVars-1),NumSheet) = num2cell(x);
ExcelData(5,23:52,NumSheet) = {'r/R','Radius','Pre-Twist','Chord','% Thick','Thickness',...
                               'Shell Thickness','Mass Dist.','Iuu','Ivv','Up. Strain','Low. Strain','LE Strain','TE Strain','Up. Stress','Low. Stress','LE Stress','TE Stress',...
                               'Mnorm','Mtan',' ','Flow Spd','Rotor Spd','Blade Pitch','Power','Power Coef','Root Flap','Torque','Thrust','Flow Probability'};
ExcelData(6,23:52,NumSheet) = {'(-)','(m)','(deg)','(m)','(t/c)','(m)',...
                               '(mm)','(kg/m)','(m^4)','(m^4)','(microstrain)','(microstrain)','(microstrain)','(microstrain)','(MPa)','(MPa)','(MPa)','(MPa)',...
                               'kN-m','kN-m',' ','(m/s)','(rpm)','(deg)','(kW)','(-)','(kN-m)','(kN-m)','(kN)','(-)'};
ExcelData(7:(NumSeg+HeaderLines),23:42,NumSheet) = num2cell(M1);
ExcelData(7:(NumCases+HeaderLines),44:52,NumSheet) = num2cell(M2);

if NumSheet == size(ParetoSet,1);
    Excel = actxserver('excel.application');
    %Excel.Visible = 1;
    ExcelWorkbook = Excel.Workbooks.Open(xls_filename);
    
    for n = 1:NumSheet
    xlswrite2007(xls_filename,ExcelData(:,:,n),num2str(n));
    end
          
    Sheets = ExcelWorkbook.Sheets; % Get the list of sheets in the workbook
    Sheets.Item('1').Select %select sheet 1 as the active sheet
    ExcelWorkbook.Save % save the changes
    ExcelWorkbook.Close(false)  % Close Excel workbook.
    Excel.Quit % close excel
    delete(Excel); 
end
%% copy the airfoil files if recording failed cases
if RecordFailures == 1
   copyfile([RootDir '\Output_Files\' filename_main '\Airfoil_Data'],[RootDir '\Output_Files\' filename_main '\Failed_Cases\Airfoil_Data']); 
end

%% Read in the Detailed_GA_Output.dat file and organize the data
% % % % % fid1 = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '_Detailed_GA_Output.dat'],'rt');
% % % % % for n = 1:2; %moves cursor through the header of the file
% % % % %     fgetl(fid1);
% % % % % end
% % % % % n=1; current_line = fgetl(fid1);
% % % % % while current_line ~= -1
% % % % % Output_temp(n,:) = str2num(current_line);
% % % % % n = n+1;
% % % % % current_line = fgetl(fid1);
% % % % % end
% % % % % fclose(fid1);
% % % % % 
% % % % % fid2 = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '_Summary_GA_Output.dat'],'Wt');
% % % % % s = [];
% % % % % for n = 1:(NumVars - 10 - (SpdCtrl==0))
% % % % %     a = abs(['Thick' num2str(n) '	']);
% % % % %     s = [s a];
% % % % % end
% % % % % format2 = char(s);
% % % % % fprintf(fid2,['F:    AEP(kW-hr/yr):	Twist1	Twist2	Twist3	Twist4	Twist5	Chord1	Chord2	Chord3	Chord4	Chord5	'	format2 '(RPM or TSR)\n']);
% % % % % fprintf(fid2,'============================================================================================================================================================\n');     
% % % % % 
% % % % % NumGens = (size(Output_temp,1)-1)/PopSize;
% % % % % for n = 0:(NumGens-1)
% % % % %     
% % % % %     Gen = Output_temp((n*PopSize+1):((n+1)*PopSize),1:end);
% % % % %     ii = find(Gen(:,1)==min(Gen(:,1)));
% % % % %     if size(ii,1) > 1; ii = ii(1); end
% % % % %     Oup(n+1,1) = Gen(ii,1);
% % % % %     %Oup(n+1,2) = Gen(ii,2);
% % % % %     %Oup(n+1,3:(2+NumVars)) = Gen(ii,3:end);
% % % % %     Oup(n+1,2:(1+NumVars)) = Gen(ii,2:end);
% % % % %     s = [];
% % % % %     for jj = 1:size(Oup,2)
% % % % %         a = [37    51    46    52   102    92   116]; %ASCII representation of: %3.4f\t
% % % % %         s = [s a];
% % % % %     end
% % % % %     s(end) = 110; %change last character to "n", (next line statment = \n)
% % % % %     format = char(s);
% % % % %     fprintf(fid2,format,Oup(n+1,1),Oup(n+1,2:end)');    
% % % % %      
% % % % % end
% % % % % fclose(fid2);
fclose('all');
%% Generate a loads document in a .txt file 

