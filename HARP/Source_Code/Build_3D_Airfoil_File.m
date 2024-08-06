function Build_3D_Airfoil_File(filename_2D,r_over_R,Chord,Thickness,Vavg,Omg)
%This function builds a stall-delay airfoil file.
%Inputs:    filename_2D: name of the 2D airfoil data file, for which the stall-delay equations will be applied to
%           r_over_R: scalar value, the non-dimensional radius
%           Chord: scalar value, the local chord
%           Thickness: scalar value, the local % thickness value
%           Vavg: scalar value, the average of the SpdSt and SpdEnd flow speeds
%           Omg: scalar value, the rotor speed
%
%Outputs:   none, after this function executes, the new airfoil files will
%           appear in the Output_Files\filename_main\Airfoil_Data\3D_Airfoil_Data directory

global RootDir filename_main Type family;

%% Determine what the 3D filename will be and open the file for writing
if Thickness >= 99.95;
    Thick_Suffix = '_1000';
elseif Thickness >= 9.95 && Thickness < 99.95;
    Thick_Suffix = ['_0' num2str(10*Thickness,'%3.0f')];
else
    Thick_Suffix = ['_00' num2str(10*Thickness,'%2.0f')];
end

if round(100*r_over_R) == 100;
        filename_3D = [RootDir '\Output_Files\' filename_main '\Airfoil_Data\3D_Airfoil_Data\' family '_1000' Thick_Suffix '.dat'];
    elseif round(100*r_over_R) <= 35;
        filename_3D = [RootDir '\Output_Files\' filename_main '\Airfoil_Data\3D_Airfoil_Data\' family '_0350' Thick_Suffix '.dat'];
    else
        filename_3D = [RootDir '\Output_Files\' filename_main '\Airfoil_Data\3D_Airfoil_Data\' family '_0' num2str(1000*r_over_R,'%3.0f') Thick_Suffix '.dat'];
end

fid3D = fopen(filename_3D,'Wt'); %first discard any previous contents
fid3D = fopen(filename_3D,'At'); %Opens file for appending data to end of file

%% Determine the number of tables in the 2D airfoil file
fid2D = fopen(filename_2D,'rt');
NumTables = cell2mat(textscan(fid2D,'%f','HeaderLines',3));

%% Write the opening header of the 2D file to the 3D file
frewind(fid2D); %rewind cursor to beginning of 2D file
for n = 1:4; fprintf(fid3D,'%s\n',fgetl(fid2D)); end

%% Read in the Aerodynamic Coefficients from the 2D airfoil file for all Reynolds numbers
for n = 1:NumTables
    
    for j = 1:8; fprintf(fid3D,'%s\n',fgetl(fid2D)); end %writes the header of the 2D table to the new 3D airfoil file 
        
    Aero_Coefficients = cell2mat(textscan(fid2D,'%f %f %f %f %f','CollectOutput',1)); %reads in the 2D aerodynamic coefficients
    if any(isnan(Aero_Coefficients(:,5))); Aero_Coefficients(:,5)=[]; end; %Delete 5th column if the user did not input this coefficient
    if any(isnan(Aero_Coefficients(:,4))); Aero_Coefficients(:,4)=[]; end; %Delete 4th column if the user did not input this coefficient
    fgetl(fid2D); %Skip through the line that reads 'EOT'

    %Apply 3D Corrections to 2D airfoil data, then append 3D data to 3D airfoil file
    if size(Aero_Coefficients,2) == 5;
      [Alpha CL_3D CD_3D C3_2D C4_2D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
      fprintf(fid3D,'%5.3f    %6.4f    %6.4f    %6.4f    %6.4f\n',[Alpha CL_3D CD_3D C3_2D C4_2D]');
      fprintf(fid3D,'EOT %s\n',num2str(n));
    elseif size(Aero_Coefficients,2) == 4;
      [Alpha CL_3D CD_3D C3_2D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
      fprintf(fid3D,'%5.3f    %6.4f    %6.4f    %6.4f\n',[Alpha CL_3D CD_3D C3_2D]'); 
      fprintf(fid3D,'EOT %s\n',num2str(n));
    else
      [Alpha CL_3D CD_3D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
      fprintf(fid3D,'%5.3f    %6.4f    %6.4f\n',[Alpha CL_3D CD_3D]'); 
      fprintf(fid3D,'EOT %s\n',num2str(n));  
    end

%     if Type == 2; %if Water turbine, expect the minimum pressure coefficients in the airfoil file
%       [Alpha CL_3D CD_3D CP_2D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
%       fprintf(fid3D,'%3.2f    %3.4f    %3.4f    %3.4f\n',[Alpha CL_3D CD_3D CP_2D]');
%       fprintf(fid3D,'EOT %s\n',num2str(n));
%     else
%       [Alpha CL_3D CD_3D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
%       fprintf(fid3D,'%3.2f    %3.4f    %3.4f\n',[Alpha CL_3D CD_3D]'); 
%       fprintf(fid3D,'EOT %s\n',num2str(n));
%     end

end
fclose(fid2D);
fclose(fid3D);
%End of Function