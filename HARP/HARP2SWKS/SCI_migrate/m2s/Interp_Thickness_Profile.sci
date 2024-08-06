function [AFcoordinates] = Interp_Thickness_Profile(AFfile1,AFfile2,ThickVal1,ThickVal2,delT)

// Output variables initialisation (not found in input variables)
AFcoordinates=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


global("RootDir","filename_main","family","NumAFPoints")

//initialize the progress bar
// !! L.6: Matlab function waitbar not yet converted, original calling sequence used.
// L.6: (Warning name conflict: function name changed from waitbar to %waitbar).
H = %waitbar(0,"","Name","Interpolating airfoil profiles...");
// !! L.7: Matlab function get not yet converted, original calling sequence used.
// L.7: (Warning name conflict: function name changed from get to %get).
hchild = %get(H,"children");
// !! L.8: Matlab function get not yet converted, original calling sequence used.
// L.8: (Warning name conflict: function name changed from get to %get).
htitle = %get(hchild,"title");
// !! L.9: Matlab function set not yet converted, original calling sequence used.
// L.9: (Warning name conflict: function name changed from set to %set).
%set(htitle,"Interpreter","None","FontSize",8);

//First need to open both airfoil coordinate files and interpolate the
//coordinates to have the same x-values
in_path1 = mtlb_fopen(RootDir+"\Input_Files\Airfoil_Data\"+AFfile1,"rt");
in_path2 = mtlb_fopen(RootDir+"\Input_Files\Airfoil_Data\"+AFfile2,"rt");
// !! L.15: Matlab function textscan not yet converted, original calling sequence used.
A1 = cell2mat(textscan(in_path1,"%f %f","HeaderLines",4,"CollectOutput",1));
// !! L.16: Matlab function textscan not yet converted, original calling sequence used.
A2 = cell2mat(textscan(in_path2,"%f %f","HeaderLines",4,"CollectOutput",1));
mclose(in_path1);
mclose(in_path2);

x1 = A1(:,1);
y1 = A1(:,2);
x2 = A2(:,1);
y2 = A2(:,2);

//Use full cosine spacing for the new x values
//NumAFPoints MUST BE AN EVEN NUMBER!
NumAFPoints = 40;//Number of points for new file, actual # of points will be NumAFPoints -1
x_upper = hcosspace(0,1,NumAFPoints/2,3);
%v0 = x_upper(1:$-1);x_lower = %v0($:-1:1,:);
Xnew = [x_upper;x_lower];

// ! L.32: abs(A1(:,1)==1) may be replaced by:
// !    --> A1(:,1)==1 if A1(:,1)==1 is Real.
a1 = find(abs(A1(:,1)==1))';
// ! L.33: abs(A2(:,1)==1) may be replaced by:
// !    --> A2(:,1)==1 if A2(:,1)==1 is Real.
a2 = find(abs(A2(:,1)==1))';

y_upper1 = interp1(mtlb_e(x1,1:a1),mtlb_e(y1,1:a1),x_upper,"pchip");
y_lower1 = interp1(mtlb_e(x1,a1:$),mtlb_e(y1,a1:$),x_lower,"pchip");
Ynew1 = [y_upper1;y_lower1];
AF1 = [Xnew,Ynew1];

y_upper2 = interp1(mtlb_e(x2,1:a2),mtlb_e(y2,1:a2),x_upper,"pchip");
y_lower2 = interp1(mtlb_e(x2,a2:$),mtlb_e(y2,a2:$),x_lower,"pchip");
Ynew2 = [y_upper2;y_lower2];
AF2 = [Xnew,Ynew2];

// subplot(2,1,1);
// plot(x1,y1,''o-b'',Xnew,Ynew1,''o-r'');
// axis equal; grid on; legend(''Original'',''Interpolated'');
// subplot(2,1,2);
// plot(x2,y2,''o-b'',Xnew,Ynew2,''o-r'');
// axis equal; grid on; legend(''Original'',''Interpolated'');

//Now write the newly interpolated profiles to the output directory
out_path1 = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\Coordinates\"+AFfile1,"Wt");
out_path2 = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\Coordinates\"+AFfile2,"Wt");
//Print the headers and data for the newly interpolated files
// L.56: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path1,"This file was generated automatically by HARP_Opt.\n");
// L.57: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path1,"These airfoil coordinates were interpolated from the file %s to have %g points, using cosine spacing.\n",AFfile1,NumAFPoints-1);
// L.58: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path1,"X\tY\n\n");
// L.59: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path1,"%3.6f\t%3.6f\n",[Xnew,Ynew1]');
// L.60: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path2,"This file was generated automatically by HARP_Opt.\n");
// L.61: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path2,"These airfoil coordinates were interpolated from the file %s to have %g points, using cosine spacing.\n",AFfile2,NumAFPoints-1);
// L.62: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path2,"X\tY\n\n");
// L.63: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf(out_path2,"%3.6f\t%3.6f\n",[Xnew,Ynew2]');
mclose(out_path1);
mclose(out_path2);

//=========================================================================%
ThickVals = mtlb_imp(mtlb_s(ThickVal1,delT),-delT,mtlb_a(ThickVal2,delT))';
if isempty(ThickVals) then
  return;  //Stop this function early, there is no need to interpolate
end;
ThickVals = [ThickVal1;ThickVals;ThickVal2];

AFcoordinates = cell();AFcoordinates(1,1) = makecell([1,1],[Xnew,Ynew1]);
AFcoordinates(max(size(ThickVals)),1) = makecell([1,1],[Xnew,Ynew2]);
AFcoordinates(1,2) = makecell([1,1],ThickVal1);
AFcoordinates(max(size(ThickVals)),2) = makecell([1,1],ThickVal2);
for n = 2:max(size(ThickVals))-1
  maxThick = 101;
  w = 1;
  while mtlb_logic(maxThick,">",mtlb_e(ThickVals,n))
    Yinterp = mtlb_a((1-w) .*Ynew2,w .*Ynew1);  //weight starts off towards y1, y1 is the thicker airfoil
    %v0_1 = mtlb_s(Yinterp,Yinterp($:-1:1,:));  maxThick = max(%v0_1,firstnonsingleton(%v0_1))*100;
    w = w-0.00005;
  end;
  MT(n,1) = maxThick;
  AFcoordinates(n,1) = makecell([1,1],[Xnew,Yinterp]);
  AFcoordinates(n,2) = makecell([1,1],mtlb_e(ThickVals,n));

  //Write the normalized airfoil coordinates each to a text file

  %v1_2 = %f;  if mtlb_logic(mtlb_e(ThickVals,n),">=",9.95) then %v1_2 = mtlb_logic(mtlb_e(ThickVals,n),"<",99.95);end;
  if mtlb_logic(mtlb_e(ThickVals,n),">=",99.95) then
    Thick_Suffix = "_1000";
  elseif %v1_2 then
    Thick_Suffix = "_0"+msprintf("%3.0f",10*mtlb_e(ThickVals,n));
  else
    Thick_Suffix = "_00"+msprintf("%3.0f",10*mtlb_e(ThickVals,n));
  end;
  AFfile(n,1:length(family+Thick_Suffix+".prof")) = family+Thick_Suffix+".prof";

  fid = mtlb_fopen(RootDir+"\Output_Files\"+filename_main+"\Airfoil_Data\Coordinates\"+part(AFfile(n),":"),"Wt");
  // L.102: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid,"This file was generated automatically by HARP_Opt.\n");
  // L.103: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid,"These airfoil coordinates represent a %g percent thick airfoil. These coordinates were interpolated between the data in the files %s and %s. %g points with cosine spacing.\n",mtlb_e(ThickVals,n),AFfile1,AFfile2,NumAFPoints-1);
  // L.104: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid,"X\tY\n\n");
  // L.105: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf(fid,"%3.6f\t%3.6f\n",[Xnew,Yinterp]');
  mclose(fid);

  // Report current progress in the waitbar''s message field
  // !! L.109: Matlab function waitbar not yet converted, original calling sequence used.
  // L.109: (Warning name conflict: function name changed from waitbar to %waitbar).
  %waitbar(n/(max(size(ThickVals))-1),H,"Interpolating between "+AFfile1+" and "+AFfile2)
end;
mtlb_delete(H);// DELETE the waitbar; don''t try to CLOSE it.

// plot(Xnew,Ynew1,''o-b'');
// for n = 1:size(AFcoordinates,1);
// hold on; 
// plot(AFcoordinates(n,1).entries(:,1),AFcoordinates(n,1).entries(:,2),''x-r''); axis equal;
// end
// plot(Xnew,Ynew2,''o-b'');
// hold off
endfunction
