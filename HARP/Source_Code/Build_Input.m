function Build_Input(flag,RElm,Twist,Chord,Thickness,SpdSt,SpdEnd,SpdDel,Omg,PitCtrl,CC)
%This function creates the WT_Perf input file for various configurations.
%flag = 0 for Cp vs. TSR, 
%flag = 1 for P vs. V
%flag = 2 for Combined Case analysis (active control turbines)

%Inputs:    flag: if flag = 0, build a Parametric Analysis to calculate the Cp vs. TSR curve
%                 if flag = 1, build a Parametric Analysis to calculate the Power vs. Velocity curve
%                 if flag = 2, build a Combined Case Analysis to calculate the Power vs. Velocity curve
%           RELm: vector, contains the radius values of the blade elements
%           Twist: vector, contains the Bezier curve Twist values of the blade elements
%           Chord: vector, contains the Bezier curve Chord values of the blade elements
%           Thickness: vector, contains the step function % Thickness values of the blade elements
%           SpdSt: scalar, the beginning flow speed
%           SpdEnd: scalar, the end flow speed
%           SpdDel: scalar, the flow speed increment
%           Omg: scalar, the rotor speed
%           PitCtrl: scalar, this is a flag which tells this function to build a Parametric Analysis that varies the blade pitch
%           CC: matrix, contains the flow speeds, rotor speeds, and blade pitches for a Combined Case Analysis
%
%Outputs:   none, after this function executes, the WT_Perf input file will
%           appear in the Output_Files\filename_main\ folder


global RootDir filename_main Correct_3D Type family SpdCtrl RotorRad StructuralOpt;

NumCases = round((SpdEnd - SpdSt)/SpdDel + 1);         
r_over_R = RElm/RotorRad;

if Type == 1 && StructuralOpt == 0
    Print_BED = 'False'; %Wind turbine AND no stuctural optimization
else
    Print_BED = 'True '; %Hydro turbine AND/OR structural optimization
end

if flag == 0;
    Print_BED = 'False'; %dont need to print the BED (dont check for cavitation at this point)
    Input_TSR = 'True ';
    Output_Pwr = 'False';
    Output_Cp = 'True ';
else
    Input_TSR = 'False';
    Output_Pwr = 'True ';
    Output_Cp = 'False';
end

if PitCtrl == 1 %Variable Pitch cases, builds a Parametric Analysis to determine the optimal pitch angles
   PitSt  = 0;
   PitEnd = 40;
   PitDel = 0.1;
   ParCol = 2;
   ParTab = 1;
elseif PitCtrl == 2 %pitch to stall turbine
   PitSt  = 0;
   PitEnd = -40;
   PitDel = -0.1;
   ParCol = 2;
   ParTab = 1; 
else
   PitSt  = 0;
   PitEnd = 0;
   PitDel = 0;
   ParCol = 1;
   ParTab = 2;
end


%=========================================================================%
%%Writes the blade geometry into the Blade Geometry.wtp file

fid = fopen('Input_Files\Templates\Blade_Geometry.wtp', 'Wt'); %opens for writing, discarding any existing text
for n = 1:length(RElm)
    fprintf(fid, '%2.4f\t%2.3f\t%2.4f\t  %d\t%s\n',RElm(n),Twist(n),Chord(n),n,'True');
end
fclose(fid);
%=========================================================================%


%=========================================================================%
%%Writes the name of the airfoil files into the Airfoils.wtp file

fid = fopen('Input_Files\Templates\Airfoils.wtp', 'Wt'); %opens for writing, discarding any existing text
fprintf(fid, '%d                                        NumAF:               Number of airfoil files.\n', length(r_over_R));

  for n = 1:length(r_over_R);
       
      if Thickness(n) >= 99.95;
        Thick_Suffix = '_1000';
      elseif Thickness(n) >= 9.95 && Thickness(n) < 99.95;
        Thick_Suffix = ['_0' num2str(10*Thickness(n),'%3.0f')];
      else
        Thick_Suffix = ['_00' num2str(10*Thickness(n),'%2.0f')];
      end;
      
      %Use the Stall Delay Models ONLY FOR THE FIXED SPEED CASES
      if ((Correct_3D == 1) || (Correct_3D == 2)) && SpdCtrl == 0;  
        filename_2D =  ['Output_Files\' filename_main '\Airfoil_Data\' family Thick_Suffix '.dat'];
        Vavg = (SpdSt + SpdEnd)/2; %Use the average flow speed for the Stall-Delay models
        
        Build_3D_Airfoil_File(filename_2D,r_over_R(n),Chord(n),Thickness(n),Vavg,Omg)
        
        
        %All radial stations inboard of r/R <= 0.35 use the same 3D corrections as the r/R = 0.35 station. This is because the stall delay models get really weird at inboard stations.
        if round(100*r_over_R(n)) == 100;
        fprintf(fid, '"Airfoil_Data\\3D_Airfoil_Data\\%s_1000%s.dat"\n', family, Thick_Suffix);
        elseif round(100*r_over_R(n)) <= 35;
        fprintf(fid, '"Airfoil_Data\\3D_Airfoil_Data\\%s_0350%s.dat"\n', family, Thick_Suffix);
        else
        fprintf(fid, '"Airfoil_Data\\3D_Airfoil_Data\\%s_0%s%s.dat"\n', family, num2str(1000*r_over_R(n),'%3.0f'), Thick_Suffix);
        end
    
      else %Don't use Stall Delay models, just use the 2D airfoil files
        fprintf(fid, '"Airfoil_Data\\%s%s.dat"\n', family, Thick_Suffix); 
      end 
  
  end
  fclose(fid);
%=========================================================================%

%=========================================================================%
%Constructs a complete WT_Perf input file by merging these files in order:
% "Input_Files\Templates\Input_Configuration.wtp"
% "Input_Files\Templates\Blade_Geometry.wtp"
% "Input_Files\Templates\Aerodynamic_Inputs.wtp"
% "Input_Files\Templates\Airfoils.wtp"
% "Input_Files\Templates\Output_Configuration.wtp

%Start Merging all files together to the file "Merged Inputs.wtp"
fid(1) = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '.wtp'],'Wt'); %discards any existing text in this file
fclose(fid(1));
fid(1) = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '.wtp'],'At'); %opens file for appending data
fid(2) = fopen([RootDir '\Input_Files\Templates\Input_Configuration.wtp'], 'rt');
fid(3) = fopen([RootDir '\Input_Files\Templates\Blade_Geometry.wtp'], 'rt');
fid(4) = fopen([RootDir '\Input_Files\Templates\Aerodynamic_Inputs.wtp'],'rt');
fid(5) = fopen([RootDir '\Input_Files\Templates\Airfoils.wtp'],'rt');

%Configures the "Output Configuration.wtp" file accordingly
if flag == 0 || flag == 1;	%This is a Parametric Analysis WT_Perf file
    fid(6) = fopen([RootDir '\Input_Files\Templates\Output_Configuration.wtp'],'r+t');
    
    for i = 1:3; fgetl(fid(6)); end;
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%s',Print_BED);
    fseek(fid(6), -count, 'cof');
    fgetl(fid(6));
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%s',Input_TSR);
    fseek(fid(6), -count, 'cof');
    for i = 1:7; fgetl(fid(6)); end;
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%g',ParCol);
    fseek(fid(6), -count, 'cof');
    fgetl(fid(6));
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%g',ParTab);
    fseek(fid(6), -count, 'cof');
    fgetl(fid(6));    
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%s',Output_Pwr);
    fseek(fid(6), -count, 'cof');
    fgetl(fid(6));
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%s',Output_Cp);
    fseek(fid(6), -count, 'cof');
    for i = 1:4; fgetl(fid(6)); end;
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%3.1f, %3.1f, %3.1f        ',PitSt,PitEnd,PitDel);
    fseek(fid(6), -count, 'cof');
    fgetl(fid(6));
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%3.2f, %3.2f, %3.2f        ',Omg,Omg,0);
    fseek(fid(6), -count, 'cof');
    fgetl(fid(6));
    fseek(fid(6),0,'cof');
    fprintf(fid(6),'%3.3f, %3.3f, %3.3f                ', SpdSt, SpdEnd, SpdDel);  

    fclose(fid(6));
    fid(6) = fopen([RootDir '\Input_Files\Templates\Output_Configuration.wtp'],'rt'); %open again for reading only now
    
else    %This is for a Combined-Case WT_Perf file
    fid(6) = fopen([RootDir '\Input_Files\Templates\CC_Output_Configuration.wtp'],'r+t');  
    fid(7) = fopen([RootDir '\Input_Files\Templates\Combined_Case.wtp'],'Wt');
    
    for i = 1:3; fgetl(fid(6)); end;
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%s',Print_BED);
    fseek(fid(6), -count, 'cof');
    fgetl(fid(6));
    fseek(fid(6),0,'cof');
    count = fprintf(fid(6),'%s',Input_TSR);
    fseek(fid(6), -count, 'cof');
    for i = 1:3; fgetl(fid(6)); end;
    fseek(fid(6),0,'cof');
    fprintf(fid(6),'%g               ',NumCases);
    
    fprintf(fid(7),'%3.3f\t%3.2f\t%3.1f\n',CC');
    fclose(fid(6));
    fclose(fid(7));
    fid(6) = fopen([RootDir '\Input_Files\Templates\CC_Output_Configuration.wtp'],'rt'); %open again for reading only now
    fid(7) = fopen([RootDir '\Input_Files\Templates\Combined_Case.wtp'],'rt'); %open again for reading only now
    
end

%Merges all files together to the file "filename_main".wtp
for n = 2:length(fid);
    current_line = 0; %initial dummy value
    while current_line ~= -1;
    current_line = fgetl(fid(n));
    if current_line == -1; %stops writing when the end of file is reached
        break
    else
    fprintf(fid(1),'%s\n', current_line);
    end
    end
end

fclose('all');
%=========================================================================%


%End of Function
