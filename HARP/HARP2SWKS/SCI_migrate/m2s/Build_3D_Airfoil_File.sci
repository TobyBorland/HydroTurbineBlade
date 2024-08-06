function [] = Build_3D_Airfoil_File(filename_2D,r_over_R,Chord,Thickness,Vavg,Omg)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This function builds a stall-delay airfoil file.
//Inputs:    filename_2D: name of the 2D airfoil data file, for which the stall-delay equations will be applied to
//           r_over_R: scalar value, the non-dimensional radius
//           Chord: scalar value, the local chord
//           Thickness: scalar value, the local % thickness value
//           Vavg: scalar value, the average of the SpdSt and SpdEnd flow speeds
//           Omg: scalar value, the rotor speed
// 
//Outputs:   none, after this function executes, the new airfoil files will
//           appear in the Output_Files\filename_main\Airfoil_Data\3D_Airfoil_Data directory

global("RootDir","filename_main","Type","family");

//% Determine what the 3D filename will be and open the file for writing
%v0_2 = %f;if mtlb_logic(Thickness,">=",9.95) then %v0_2 = mtlb_logic(Thickness,"<",99.95);end;
if mtlb_logic(Thickness,">=",99.95) then
  Thick_Suffix = "_1000";
elseif %v0_2 then
  Thick_Suffix = "_0"+msprintf("%3.0f",10*Thickness);
else
  Thick_Suffix = "_00"+msprintf("%2.0f",10*Thickness);
end;

if mtlb_logic(round(100*r_over_R),"==",100) then
  filename_3D = RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\3D_Airfoil_Data\"+family+"_1000"+Thick_Suffix+".dat";
elseif mtlb_logic(round(100*r_over_R),"<=",35) then
  filename_3D = RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\3D_Airfoil_Data\"+family+"_0350"+Thick_Suffix+".dat";
else
  filename_3D = RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\3D_Airfoil_Data\"+family+"_0"+msprintf("%3.0f",1000*r_over_R)+Thick_Suffix+".dat";
end;

fid3D = mtlb_fopen(filename_3D,"Wt");//first discard any previous contents
fid3D = mtlb_fopen(filename_3D,"At");//Opens file for appending data to end of file

//% Determine the number of tables in the 2D airfoil file
fid2D = mtlb_fopen(filename_2D,"rt");
// !! L.37: Matlab function textscan not yet converted, original calling sequence used.
NumTables = cell2mat(textscan(fid2D,"%f","HeaderLines",3));

//% Write the opening header of the 2D file to the 3D file
mseek(0,fid2D);//rewind cursor to beginning of 2D file
for n = 1:4 %v0 = mgetl(fid2D,1); if meof()~=0 then %v0 = -1;end; // L.41: No simple equivalent, so mtlb_fprintf() is called.
 mtlb_fprintf(fid3D,"%s\n",%v0);end;

//% Read in the Aerodynamic Coefficients from the 2D airfoil file for all Reynolds numbers
for n = mtlb_imp(1,NumTables)

  for j = 1:8 %v0 = mgetl(fid2D,1); if meof()~=0 then %v0 = -1;end; // L.46: No simple equivalent, so mtlb_fprintf() is called.
   mtlb_fprintf(fid3D,"%s\n",%v0);end;  //writes the header of the 2D table to the new 3D airfoil file 

  // !! L.48: Matlab function textscan not yet converted, original calling sequence used.
  Aero_Coefficients = cell2mat(textscan(fid2D,"%f %f %f %f %f","CollectOutput",1));  //reads in the 2D aerodynamic coefficients
  // ! L.49: abs(isnan(Aero_Coefficients(:,5))) may be replaced by:
  // !    --> isnan(Aero_Coefficients(:,5)) if isnan(Aero_Coefficients(:,5)) is Real.

  if mtlb_any(abs(isnan(Aero_Coefficients(:,5)))) then Aero_Coefficients(:,5) = [];end;  //Delete 5th column if the user did not input this coefficient
  // ! L.50: abs(isnan(Aero_Coefficients(:,4))) may be replaced by:
  // !    --> isnan(Aero_Coefficients(:,4)) if isnan(Aero_Coefficients(:,4)) is Real.

  if mtlb_any(abs(isnan(Aero_Coefficients(:,4)))) then Aero_Coefficients(:,4) = [];end;  //Delete 4th column if the user did not input this coefficient
  %v1 = mgetl(fid2D,1);  if meof()~=0 then %v1 = -1;end;  %v1;  //Skip through the line that reads ''EOT''

  //Apply 3D Corrections to 2D airfoil data, then append 3D data to 3D airfoil file
  if size(Aero_Coefficients,2)==5 then
    [Alpha,CL_3D,CD_3D,C3_2D,C4_2D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
    // L.56: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3D,"%5.3f    %6.4f    %6.4f    %6.4f    %6.4f\n",mtlb_t([Alpha,CL_3D,CD_3D,C3_2D,C4_2D]));
    // !! L.57: string output can be different from Matlab num2str output.
    // L.57: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3D,"EOT %s\n",string(n));
  elseif size(Aero_Coefficients,2)==4 then
    [Alpha,CL_3D,CD_3D,C3_2D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
    // L.60: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3D,"%5.3f    %6.4f    %6.4f    %6.4f\n",mtlb_t([Alpha,CL_3D,CD_3D,C3_2D]));
    // !! L.61: string output can be different from Matlab num2str output.
    // L.61: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3D,"EOT %s\n",string(n));
  else
    [Alpha,CL_3D,CD_3D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
    // L.64: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3D,"%5.3f    %6.4f    %6.4f\n",mtlb_t([Alpha,CL_3D,CD_3D]));
    // !! L.65: string output can be different from Matlab num2str output.
    // L.65: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf(fid3D,"EOT %s\n",string(n));
  end;

  //     if Type == 2; %if Water turbine, expect the minimum pressure coefficients in the airfoil file
  //       [Alpha CL_3D CD_3D CP_2D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
  //       fprintf(fid3D,''%3.2f    %3.4f    %3.4f    %3.4f\n'',[Alpha CL_3D CD_3D CP_2D]'');
  //       fprintf(fid3D,''EOT %s\n'',num2str(n));
  //     else
  //       [Alpha CL_3D CD_3D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg);
  //       fprintf(fid3D,''%3.2f    %3.4f    %3.4f\n'',[Alpha CL_3D CD_3D]''); 
  //       fprintf(fid3D,''EOT %s\n'',num2str(n));
  //     end

end;
mclose(fid2D);
mclose(fid3D);
//End of Function
endfunction
