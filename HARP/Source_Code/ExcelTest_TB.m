function Optimize_Callback(hObject, eventdata, handles)
% hObject    handle to Optimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% clc;
% clear functions %reinitalizes the persistent variables
% 
% fprintf(1,'User inputs defined. Initializing...\n')
% 
% %Calls the "Save" function which writes all the WT_Perf template files
% Save_Callback(hObject, eventdata, handles);

global SpdSt SpdEnd SpdDel OmgMin OmgMax Prated rho Patm Pv HubHt WatDepth Type...
       NumSeg HubRad Thickness_values ThickMethod RotorRad NumBlade family RootDir filename_main NumAFPoints...
       Correct_3D SpdCtrl PitCtrl CavSF OptMethod StructuralOpt ProbDist U_mean Weib_Umean Weib_k Weib_c user_pU_interp...
       ShaftTilt PreCone NumVars NumGens PopSize EliteCount ParetoFraction CrossFrc KinVisc CircleRoot...
       minRootChord maxRootChord RootTranSt RootTranEnd Normalized_AF_Coordinates...
       TwistLB TwistUB ChordLB ChordUB ThickLB ThickUB LB UB...
       Einput allowableStrain matDensity DecreasingST STmin STdel SFstruct...
       ParetoSet Pareto_filename OpenFile_Error RadialSpacing RecordFailures; %NumAFPoints added

% %Read in the input variables from the GUI
% SpdCtrl = str2num(handles.SpdCtrl);
% PitCtrl = str2num(handles.PitCtrl);
% NumBlade = str2num(get(handles.NumBlade,'String'));
% NumSeg = str2num(get(handles.NumSeg,'String'));
% RotorRad = 0.5*str2num(get(handles.RotorDia,'String'));
% HubRad = 0.5*str2num(get(handles.HubDia,'String'));
% HubHt = str2num(get(handles.HubHt,'String'));
% Prated = str2num(get(handles.RatedPwr,'String'));
% ShaftTilt = str2num(get(handles.ShaftTilt,'String'));
% PreCone = str2num(get(handles.PreCone,'String'));
% SpdSt = str2num(get(handles.SpdSt,'String'));
% SpdEnd = str2num(get(handles.SpdEnd,'String'));
% SpdDel = str2num(get(handles.SpdDel,'String'));
% OmgMin = str2num(get(handles.OmgMin,'String'));
% OmgMax = str2num(get(handles.OmgMax,'String'));
% rho = str2num(get(handles.Rho,'String'));
% KinVisc = str2num(get(handles.KinVisc,'String'));
% CC = get(handles.CheckCavit,'Value');
% Type = CC+1;
% Pv = str2num(get(handles.VapPres,'String'));
% Patm = str2num(get(handles.Patm,'String'));
% WatDepth = str2num(get(handles.WatDepth,'String'));
% CavSF = str2num(get(handles.CavSF,'String'));
% Correct_3D = str2num(handles.dragmodel);
% ThickMethod = str2num(handles.ThickMethod);
% CircleRoot = get(handles.CircularRoot,'Value');
% minRootChord = str2num(get(handles.minRtChord,'String'));
% maxRootChord = str2num(get(handles.maxRtChord,'String'));
% RootTranSt = str2num(get(handles.RtTranSt,'String'));
% RootTranEnd = str2num(get(handles.RtTranEnd,'String'));
% 
% StructuralOpt = get(handles.StructuralOpt,'Value');
% RecordFailures = get(handles.RecFail,'Value');
% 
% %blade element spacing
% RadialSpacing = str2num(handles.Elm_spacing);
% 
% Einput = str2num(get(handles.Einput,'String')); %GPa bulk material elasticity
% allowableStrain = str2num(get(handles.allowableStrain,'String'))./1e6; %user input in microstrain, convert to strain
% matDensity = str2num(get(handles.matDensity,'String')); %(kg/m^3) bulk material density
% STmin = str2num(get(handles.STmin,'String'))/1000; %minimum shell thickness (m)
% STdel = str2num(get(handles.STdel,'String'))/1000; %shell thickness increment (m)
% SFstruct = str2num(get(handles.SFstruct,'String')); %safety factor multiplied to bending moments
% DecreasingST = str2num(handles.DecreasingST);
% 
% OptMethod = str2num(handles.OptMethod);
% ProbDist = str2num(handles.ProbDist);
% opened_file = get(handles.FlowDist_Filename,'String');    
%     if strcmp(opened_file,'no file selected')
%         OpenFile_Error = 1;
%     else
%         OpenFile_Error = 0;
%     end
% 
% U_mean = str2num(get(handles.U_box,'String'));
% Weib_Umean = str2num(get(handles.Weib_U_box,'String'));
% Weib_k = str2num(get(handles.k_box,'String'));
% Weib_c = str2num(get(handles.c_box,'String'));
%     
% NumGens = str2num(get(handles.NumGens,'String'));
% PopSize = str2num(get(handles.PopSize,'String'));
% EliteCount = str2num(get(handles.EliteCount,'String'));
% ParetoFraction = str2num(get(handles.ParetoFraction,'String'));
% CrossFrc = str2num(get(handles.CrossFrc,'String'));
% GATol = str2num(get(handles.GATol,'String'));
% TwistLB = [str2num(get(handles.TwLB1,'String')) str2num(get(handles.TwLB2,'String')) str2num(get(handles.TwLB3,'String')) str2num(get(handles.TwLB4,'String')) str2num(get(handles.TwLB5,'String'))];
% TwistUB = [str2num(get(handles.TwUB1,'String')) str2num(get(handles.TwUB2,'String')) str2num(get(handles.TwUB3,'String')) str2num(get(handles.TwUB4,'String')) str2num(get(handles.TwUB5,'String'))];
% ChordLB = [str2num(get(handles.CLB1,'String')) str2num(get(handles.CLB2,'String')) str2num(get(handles.CLB3,'String')) str2num(get(handles.CLB4,'String')) str2num(get(handles.CLB5,'String'))];
% ChordUB = [str2num(get(handles.CUB1,'String')) str2num(get(handles.CUB2,'String')) str2num(get(handles.CUB3,'String')) str2num(get(handles.CUB4,'String')) str2num(get(handles.CUB5,'String'))];
% ThickLB = str2num(get(handles.ThickLB,'String'));
% ThickUB = str2num(get(handles.ThickUB,'String'));

SpdCtrl = 1;
PitCtrl = 0;
NumBlade = 3;
NumSeg = 30;
RotorRad = 0.5;
HubRad = 0.12;
HubHt = 6;
Prated = 1.0;
ShaftTilt = 0.0;
PreCone = 0.0;
SpdSt = 0.2;
SpdEnd = 3.5;
SpdDel = 0.1;
OmgMin = 4;
OmgMax = 90;
rho = 1025.0;
KinVisc = 0.0000010300;
CC = 1;
Type = CC+1;
Pv = 2500.0;
Patm = 101325.0;
WatDepth = 7.0;
CavSF = 1.0;
Correct_3D = 0;
ThickMethod = 2;
CircleRoot = 0;
minRootChord = 0.8;
maxRootChord = 1.5;
RootTranSt = 1.6667;
RootTranEnd = 0.25;

StructuralOpt = 1;
RecordFailures = 0;

%blade element spacing
RadialSpacing = 0;

Einput = 27.6; %GPa bulk material elasticity
allowableStrain = 3000/1e6; %user input in microstrain, convert to strain
matDensity = 1800; %(kg/m^3) bulk material density
STmin = 1/1000; %minimum shell thickness (m)
STdel = 0.2/1000; %shell thickness increment (m)
SFstruct = 1.2; %safety factor multiplied to bending moments
DecreasingST = 1;

OptMethod = 0;
ProbDist = 0;

U_mean = 6.03;
Weib_Umean = 6.03;
Weib_k = 1.91;
Weib_c = 6.8;
    
NumGens = 60;
PopSize = 80;
EliteCount = 1;
ParetoFraction = 0.5;
CrossFrc = 0.25;
GATol = 1.0e-03;
TwistLB = [-10 -10 -10 -10 -10];
TwistUB = [40 40 40 40 40];
ChordLB = [0.25 0.2 0.1 0.1 0.1];
ChordUB = [0.4 0.4 0.4 0.4 0.4];

LB = [-10.0 -10.0 -10.0 -10.0 -10.0 0.25]; %gets passed into ga function
UB = [40.0 40.0 40.0 40.0 40.0 0.4000]; %gets passed into ga function
         
family = 'BW3';
Thickness_values = 5;  %Thickness values must be in descending order
filename_main = '1W5_BW3_1m_GRP3';
    
%Create a new output directory for the user's filename, and copy the WT_Perf executable to this directory
RootDir = pwd; %pwd is a function which returns the present working directory

% %Check if output directory already exists, if yes, then delete and recreate
% if exist([RootDir '\Output_Files\' filename_main],'dir') == 7;
%    rmdir([RootDir '\Output_Files\' filename_main],'s');
%    mkdir([RootDir '\Output_Files\' filename_main]);
% else
%    mkdir([RootDir '\Output_Files\' filename_main]);
% end

% %save a screenshot of the GUI to output directory
% saveas(gcf,[RootDir '\Output_Files\' filename_main '\' filename_main '_Input.bmp'],'bmp');
% 
% mkdir([RootDir '\Output_Files\' filename_main '\Airfoil_Data']);
% copyfile([RootDir '\Input_Files\Templates\WT_Perf_Jan1.exe'],[RootDir '\Output_Files\' filename_main '\WT_Perf_Jan1.exe']);
% 
%Copy the template for the Excel output file to the new output directory
copyfile([RootDir '\Input_Files\Templates\Output_Template.xls'],[RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xls']);
% 
% %If a Stall-Delay model is being used, need to create a directory for the 3D Airfoil Data as well
% if Correct_3D ~= 0
% mkdir([RootDir '\Output_Files\' filename_main '\Airfoil_Data\3D_Airfoil_Data'])
% end
% 
% %If we are recording the failed cases, make new folder
% if RecordFailures == 1
%    mkdir([RootDir '\Output_Files\' filename_main '\Failed_Cases'])
%    copyfile([RootDir '\Input_Files\Templates\WT_Perf_Jan1.exe'],[RootDir '\Output_Files\' filename_main '\Failed_Cases\WT_Perf_Jan1.exe']);
% end



%Establish the Linear Inequality Constraints Ax<=b. Establish constraints
%such that the twist, chord, and %thickness distributions are monotonically decreasing
A = [-1 1 0 0 0 0 0 0 0 0;   %c1>=c2
      0 -1 1 0 0 0 0 0 0 0;  %c2>=c3
      0 0 -1 1 0 0 0 0 0 0;  %c3>=c4
      0 0 0 -1 1 0 0 0 0 0;  %c4>=c5
      0 0 0 0 0 -1 1 0 0 0;  %tw1>=tw2
      0 0 0 0 0 0 -1 1 0 0;  %tw2>=tw3
      0 0 0 0 0 0 0 -1 1 0;  %tw3>=tw4
      0 0 0 0 0 0 0 0 -1 1]; %tw4>=tw5

if length(Thickness_values) == 1; %only 1 airfoil profile exists
    NumVars = 10;
elseif length(Thickness_values) > 1 && ThickMethod == 1 %Piecewise Constant
    NumVars = 10 + length(Thickness_values)-1; %multiple airfoil profiles exist
    A(:,NumVars) = 0;
    nn=0;
    for n = 9:6+length(Thickness_values)
    A(n,11+nn) = 1;
    A(n,11+nn+1) = -1;
    nn=nn+1;
    end  
elseif length(Thickness_values) > 1 && ThickMethod == 2 %Piecewise Linear
    NumVars = 10 + length(Thickness_values); %multiple airfoil profiles exist
    nn=0;
    for n = 9:7+length(Thickness_values)
    A(n,11+nn) = 1;
    A(n,11+nn+1) = -1;
    nn=nn+1;
    end 
end 
b = zeros(size(A,1),1); %Linear Inequality Constraints Ax<=b

if SpdCtrl == 0   %Fixed Speed: rotor speed becomes the last variable
   NumVars = NumVars + 1; 
   LB = [LB OmgMin]; %add the min/max rotors speeds into the bounds
   UB = [UB OmgMax];
   A(:,NumVars) = 0;
end

    %=====================================================================%
% %     NumSeedIndiv = xlsread('Input_Files\Initial_Population.xls',1,'B6');
% %     SeedInitPop = get(handles.SeedInitPop,'Value');
% %     
% %     %Read in Initial Populations from the "Initial Population.xls" file
% %     %NOTE: can only seed individuals with the same number of variables as
% %     %the user specified in the current HARP_Opt session
% %     
% %     if SeedInitPop == 1 && NumSeedIndiv ~= 0;
% %         InitPop = xlsread('Input_Files\Initial_Population.xls',1,['A10:P' num2str(10 + NumSeedIndiv - 1)]);
% %         if NumVars == 16;       %Rotor Speed is the 16th variable       
% %         EmptyCells = find(isnan(InitPop));
% %         OmgAvg = (OmgMin + OmgMax)/2;
% %         InitPop(EmptyCells) = OmgAvg;  
% %             if size(InitPop,2) < 16; %the user has neglected to seed a 16th variable
% %             InitPop(:,16) = OmgAvg;  %If the inital population does not contain a value for Rotor Spd., use this average value
% %             end
% %         else                       %Only need 15 variables
% %         InitPop = InitPop(:,1:15); %get rid of the 16th variables
% %         end
% %     else
         InitPop = [];
% %     end
    %=====================================================================%

%Configure the options for the Genetic Algorithm

%'PopInitRange',[LB;UB], ---
%{@gaplotbestindiv,@gaplotbestf,@gaplotdistance,
GAoptions = gaoptimset('InitialPopulation',InitPop,'CreationFcn',{@gaCustom_Creation},'Generations',NumGens,'PopulationSize',PopSize,'EliteCount',EliteCount,...
                     'CrossoverFraction',CrossFrc,'TolFun',GATol,'CrossoverFcn',{@crossoverintermediate},...
                     'SelectionFcn',{@selectionstochunif},'FitnessScalingFcn',{@fitscalingrank},'MutationFcn',{@mutationadaptfeasible},...
                     'PlotFcns',{@gaCustomPlot},'OutputFcns',@gaoutputfcn,'Display','diagnose');
disp('GAoptions:');disp(GAoptions);

MOGAoptions = gaoptimset('InitialPopulation',InitPop,'CreationFcn',{@gaCustom_Creation},'Generations',NumGens,'PopulationSize',PopSize,...
                     'CrossoverFraction',CrossFrc,'TolFun',GATol,'CrossoverFcn',{@crossoverintermediate},...
                     'MutationFcn',{@mutationadaptfeasible},'ParetoFraction',ParetoFraction,'PlotFcn',@gaCustomPlot,...
                     'OutputFcns',@gaoutputfcn,'Display','diagnose');
disp('MOGAoptions:');disp(MOGAoptions);
                 
%Need to check for user input errors before starting the Genetic Algorithm
UserInputError = Error_Check; %Calls the Error_Check.m function
if UserInputError == 0;
    
    %=====================================================================%
%     %Read in variables and perform tasks related to AEP
%        if ProbDist == 3;
%            %User Defined Distribution: read in data from "Flow_Distribution.dat"
%            fid = fopen(get(handles.Hidden_FlowDist_Filename,'String'),'rt');
%            user_data = textscan(fid,'%f %f','HeaderLines',13);
%            fclose(fid);
%            U_custom = user_data{1};
%            pU_custom = user_data{2};
%            %Now interpolate the custom p(U) distribution
%            V = (SpdSt:SpdDel:SpdEnd)';
%            user_pU_interp = interp1(U_custom,pU_custom,V,'pchip');
%            %if the spline caused any negative values change them to zero
%            user_pU_interp(user_pU_interp<0) = 0;
%        else
%            user_pU_interp = [];    
%        end  

    %=====================================================================%
 
    
%     %If the airfoil files can be found (i.e. there were no user input errors,
%     %then copy the 2D airfoil files that will be used to the new ouput directory
%     if CircleRoot == 1;
%         ThickVals = [100 Thickness_values];
%     else
%         ThickVals = Thickness_values;
%     end    
%     
%     for n = 1:length(ThickVals)
%         if ThickVals(n) >= 99.95;
%         Thick_Suffix = '_1000';
%         elseif ThickVals(n) >= 9.95 && ThickVals(n) < 99.95;
%         Thick_Suffix = ['_0' num2str(10*ThickVals(n),'%3.0f')];
%         else  
%         Thick_Suffix = ['_00' num2str(10*ThickVals(n),'%3.0f')];
%         end
% 
% %         AFfile(n,:) = [family Thick_Suffix '.dat'];
%         AFfile(n,:) = [family Thick_Suffix];
%         copyfile([RootDir '\Input_Files\Airfoil_Data\' AFfile(n,:) '.dat'],[RootDir '\Output_Files\' filename_main '\Airfoil_Data\' AFfile(n,:) '.dat']);
%     end

%     %We need to create new airfoil files which are interpolated between the available airfoil percent thicknesses
%     thick_stepsize = 0.1; %airfoil thickness will be interpolated using this interval stepsize
%     if CircleRoot == 0 && ThickMethod == 1
%         %do nothing
%     else
%         for n = 1:length(ThickVals)-1;
%         Interp_Thickness_Coefs([AFfile(n,:) '.dat'],[AFfile(n+1,:) '.dat'],ThickVals(n),ThickVals(n+1),thick_stepsize); 
%         end    
%     end
    
%     if StructuralOpt == 1;
%        mkdir([RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates']);
%           
%         if length(ThickVals) > 1
% 
%            for n = 1:length(ThickVals)-1;    
%                 [AFcoordinates] = Interp_Thickness_Profile([AFfile(n,:) '.prof'],[AFfile(n+1,:) '.prof'],ThickVals(n),ThickVals(n+1),thick_stepsize);
%                 if n == 1
%                 Normalized_AF_Coordinates = AFcoordinates;
%                 else
%                 Normalized_AF_Coordinates = [Normalized_AF_Coordinates;AFcoordinates(2:end,:)];
%                 end
%            end
%        elseif length(ThickVals) == 1;
%             %only one airfoil being used
%             in_path1 = fopen([RootDir '\Input_Files\Airfoil_Data\' AFfile '.prof'],'rt');
%             A1 = cell2mat(textscan(in_path1,'%f %f','HeaderLines',4,'CollectOutput',1));
%             fclose(in_path1);
%             x1 = A1(:,1); %x coordinates
%             y1 = A1(:,2); %y coordinates
%             NumAFPoints = 40; %Number of points for new file, actual # of points will be NumAFPoints -1
%             x_upper = hcosspace(0,1,NumAFPoints/2,3);
%             x_lower = flipud(x_upper(1:(end-1)));
%             Xnew = [x_upper;x_lower];
%             a1 = find(A1(:,1)==1);
%             y_upper1 = interp1(x1(1:a1),y1(1:a1),x_upper,'pchip');
%             y_lower1 = interp1(x1(a1:end),y1(a1:end),x_lower,'pchip');
%             Ynew1 = [y_upper1;y_lower1];
%             Normalized_AF_Coordinates = {[Xnew Ynew1], ThickVals};
%             %Now write the newly interpolated profile to the output directory
%             out_path1 = fopen([RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates\' AFfile '.prof'],'Wt');
%             %Print the headers and data for the newly interpolated files
%             fprintf(out_path1,'This file was generated automatically by HARP_Opt.\n');
%             fprintf(out_path1,'These airfoil coordinates were interpolated from the file %s to have %g points, using cosine spacing.\n',[AFfile '.prof'],NumAFPoints-1);
%             fprintf(out_path1,'X\tY\n\n');
%             fprintf(out_path1,'%3.6f\t%3.6f\n',[Xnew Ynew1]');
%             fclose(out_path1);
%         end
%           coords_filename = [RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates\Normalized_AF_Coordinates'];
%           save(coords_filename,'Normalized_AF_Coordinates')
%     end
    
     if StructuralOpt == 1;
%   
%     [ParetoSet,Fvals,exitFlag,Output] = gamultiobj(@Main,NumVars,A,b,[],[],LB,UB,MOGAoptions);
% 
%     %sort the Pareto set
%     [sortedFvals sortIndex] = sortrows(Fvals);
%     sortedParetoSet = zeros(size(ParetoSet));
%     for n = 1:length(sortIndex)
%         sortedParetoSet(n,:) = ParetoSet((sortIndex(n)),:);
%     end
%     
%     Fvals = sortedFvals;
%     ParetoSet = sortedParetoSet;%pause on
%     
%     NaNcols = isnan(Fvals(:,1));
%     Fvals(NaNcols,:) = [];
%     ParetoSet(NaNcols,:) = [];
%     
%     %disp('ParetoSet');
%     if size(ParetoSet,1) > 1
%         for n = 1:size(ParetoSet,1)
%              disp('ParetoSet');disp(n);
%             disp(ParetoSet(n,:));
%         end
%     end
%     %pause on;
%     %keyboard;
%=========================================================================%
    
   ParetoSet = ...
   [34.6969   29.4882   17.9729    8.4780    2.9690    0.2720    0.2250    0.1592    0.1166    0.1019
    34.6964   32.0363   18.6594    6.4919    3.1431    0.2716    0.2249    0.1590    0.1166    0.1013
    34.6964   32.0363   18.6594    6.4919    3.1431    0.2716    0.2249    0.1590    0.1166    0.1013
    32.6781   30.9619   17.9872    5.9116    2.3195    0.2717    0.2249    0.1588    0.1163    0.1013
    33.9249   32.6522   18.6207    4.7558    2.5051    0.2717    0.2248    0.1587    0.1160    0.1013
    33.9249   32.6522   18.6207    4.7558    2.5051    0.2717    0.2248    0.1587    0.1160    0.1013
    33.9249   32.6522   18.6207    4.7558    2.5051    0.2717    0.2248    0.1587    0.1160    0.1013
    31.8266   31.7211   17.2397    4.1655    2.4562    0.2717    0.2247    0.1581    0.1157    0.1012
    32.0593   27.0749   16.5108    4.2979    0.5544    0.2717    0.2243    0.1582    0.1157    0.1012
    32.0593   27.0749   16.5108    4.2979    0.5544    0.2717    0.2243    0.1582    0.1157    0.1012
    30.7040   27.0359   16.2544    4.1316    0.6311    0.2717    0.2243    0.1582    0.1156    0.1012
    30.2443   27.1600   15.8707    4.0953    0.7476    0.2717    0.2242    0.1581    0.1150    0.1012
    30.2443   27.1600   15.8707    4.0953    0.7476    0.2717    0.2242    0.1581    0.1150    0.1012
    28.0167   25.4638   16.2510    3.8305    0.6489    0.2714    0.2239    0.1544    0.1140    0.1013
    27.3676   24.7951   16.1383    3.7847    0.4123    0.2714    0.2239    0.1542    0.1137    0.1013
    27.2890   24.6931   15.6412    3.5687    0.5206    0.2714    0.2238    0.1539    0.1137    0.1013];
%     8.3042    -5.9165   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.9165   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.0542    -5.8540   -8.3431   -9.1857   -9.5448    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.9165   -8.3431   -9.4982   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -6.1665   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.6665   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.9165   -8.3431   -9.4982   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     7.8042    -6.1665   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.4165   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.0542    -5.8540   -8.3431   -9.1857   -9.5448    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.6665   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.6665   -8.3431   -9.7482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.5542    -5.4165   -8.3431   -9.2482   -9.5448    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.5542    -5.6665   -8.3431   -8.7482   -9.2948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.9165   -8.3431   -9.3232   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.9165   -7.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     7.5542    -5.4165   -8.3431   -9.2482   -9.5448    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -5.6665   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.0542    -5.8540   -8.3431   -9.1857   -9.5448    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.5542    -5.6665   -8.3431   -8.7482   -9.2948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.3042    -6.1665   -8.3431   -9.2482   -9.7948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.5542    -4.6665   -8.3431   -8.7482   -9.2948    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.4260    -5.6665   -8.3431   -8.9893   -9.4707    0.2707    0.2195    0.1330    0.1019    0.1010
%     8.0542    -5.8540   -8.3431   -9.1857   -9.5448    0.2707    0.2195    0.1330    0.1019    0.1010];


    %dumpfile = [RootDir '\Output_Files\' filename_main '\' filename_main  '_dump.txt'];
    %save dumpfile cell2mat(ParetoSet) -ASCII
        
    %Start writing the Pareto data to Excel file
    Excel = actxserver('excel.application');
    Excel.Visible = 1;
    Excel.DisplayAlerts = 0;
    ExcelFile = [RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xls'];
    % MS EXCEL 2007
    %ExcelFileTemplate = [RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xltx'];
    % MS EXCEL 2003
    %ExcelFileTemplate = [RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xlt'];
    %open the excel file, full path needs to be mentioned or else excel will pick it from most recently opened files.
    ExcelWorkbook = Excel.Workbooks.Open(ExcelFile); %This file has a sheet named "1" which has some formating
    wksheet = ExcelWorkbook.Worksheets.Item('1'); % Choose desired sheet
    
    
    
    %Make a copy of the Excel template sheet for the multiple Pareto solutions
    if size(ParetoSet,1) > 1
        for n = 2:size(ParetoSet,1)
            if n > 26
                break;
                disp();
            end
            %make a copy of the template sheet % PROBLEM W/ EXCEL 2000,
            %2007.. does this fail on >19 sets?
            % copy failure recorded: see support microsoft.com/kb/210684
            wksheet.Copy(wksheet); %this will create a sheet called "1 (2)" and places it before "1" (not sure if this will be consistent)
            %wksheet.Add Type=ExcelFileTemplate
            disp(' ..generating Excel workbook sheet');
            newSheet=ExcelWorkbook.Worksheets.Item('1 (2)'); %get a handle to this copied sheet
            newSheet.Name=num2str(n); %rename it with a new name
            if mod(n,10)==0
                % copy method still fails at over 26 instances or so
                ExcelWorkbook.Save 
                ExcelWorkbook.Close(false)  
                ExcelWorkbook = Excel.Workbooks.Open(ExcelFile); 
                wksheet = ExcelWorkbook.Worksheets.Item('1');
            end
        end
        
        % Get the list of sheets in the workbook
        Sheets = ExcelWorkbook.Sheets;
        % Moves Sheet "1" to before sheet "2"
        Sheets.Item('1').Move(Sheets.Item('2'))          
    end
    
    ExcelWorkbook.Save % save the changes
    ExcelWorkbook.Close(false)  % Close Excel workbook.
    Excel.Quit % close excel
    delete(Excel); 
    %system('taskkill /F /IM EXCEL.EXE');
    
    %Ok, now we have copied all the sheets we needed, now let's go through
    %the sheets one by one and write the data
    for n = 1:size(ParetoSet,1)
        Pareto_filename = [filename_main '_' num2str(n)]; %this is a global variable which gets changed, Post_Process.m uses this variable
        Main([ParetoSet(n,:)';n]);
    end
    
    else
    %Initialze the Genetic Algorithm and return the best individual found
    [x_final] = ga(@Main,NumVars,A,b,[],[],LB,UB,[],GAoptions);
    
    Pareto_filename = filename_main; %this is a global variable which gets changed, Post_Process.m uses this variable
    ParetoSet = x_final;
    Main([x_final 1]);
    end
    


    %Append the best individual to the "Initial Population.xls" file
%     xlswrite('Input_Files\Initial_Population.xls',NumSeedIndiv+1,'B6:B6');
%     if NumVars == 15;
%     xlswrite('Input_Files\Initial_Population.xls',x_final',['A' num2str(10 + NumSeedIndiv) ':O' num2str(10 + NumSeedIndiv)]);
%     else
%     xlswrite('Input_Files\Initial_Population.xls',x_final',['A' num2str(10 + NumSeedIndiv) ':P' num2str(10 + NumSeedIndiv)]);
%     end

    disp('...Optimization Complete.'); %The entire program has finished running!
    if StructuralOpt == 1
    fprintf('The number of points on the Pareto front was: %d\n', size(ParetoSet,1));
    end
else
    disp('User input error. Program terminated.');
    fclose('all');
end
%=========================================================================%
