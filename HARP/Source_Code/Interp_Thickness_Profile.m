function [AFcoordinates] = Interp_Thickness_Profile(AFfile1,AFfile2,ThickVal1,ThickVal2,delT)

global RootDir filename_main family NumAFPoints

%initialize the progress bar
H = waitbar(0,'','Name','Interpolating airfoil profiles...');
hchild = get(H,'children');
htitle = get(hchild,'title');
set(htitle,'Interpreter','None','FontSize',8);

%First need to open both airfoil coordinate files and interpolate the
%coordinates to have the same x-values
in_path1 = fopen([RootDir '\Input_Files\Airfoil_Data\' AFfile1],'rt');
in_path2 = fopen([RootDir '\Input_Files\Airfoil_Data\' AFfile2],'rt');
A1 = cell2mat(textscan(in_path1,'%f %f','HeaderLines',4,'CollectOutput',1));
A2 = cell2mat(textscan(in_path2,'%f %f','HeaderLines',4,'CollectOutput',1));
fclose(in_path1);
fclose(in_path2);

x1 = A1(:,1);
y1 = A1(:,2);
x2 = A2(:,1);
y2 = A2(:,2);

%Use full cosine spacing for the new x values
%NumAFPoints MUST BE AN EVEN NUMBER!
NumAFPoints = 40; %Number of points for new file, actual # of points will be NumAFPoints -1
x_upper = hcosspace(0,1,NumAFPoints/2,3);
x_lower = flipud(x_upper(1:(end-1)));
Xnew = [x_upper;x_lower];

a1 = find(A1(:,1)==1);
a2 = find(A2(:,1)==1);

y_upper1 = interp1(x1(1:a1),y1(1:a1),x_upper,'pchip');
y_lower1 = interp1(x1(a1:end),y1(a1:end),x_lower,'pchip');
Ynew1 = [y_upper1;y_lower1];
AF1 = [Xnew Ynew1];

y_upper2 = interp1(x2(1:a2),y2(1:a2),x_upper,'pchip');
y_lower2 = interp1(x2(a2:end),y2(a2:end),x_lower,'pchip');
Ynew2 = [y_upper2;y_lower2];
AF2 = [Xnew Ynew2];

% subplot(2,1,1);
% plot(x1,y1,'o-b',Xnew,Ynew1,'o-r');
% axis equal; grid on; legend('Original','Interpolated');
% subplot(2,1,2);
% plot(x2,y2,'o-b',Xnew,Ynew2,'o-r');
% axis equal; grid on; legend('Original','Interpolated');

%Now write the newly interpolated profiles to the output directory
out_path1 = fopen([RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates\' AFfile1],'Wt');
out_path2 = fopen([RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates\' AFfile2],'Wt');
%Print the headers and data for the newly interpolated files
fprintf(out_path1,'This file was generated automatically by HARP_Opt.\n');
fprintf(out_path1,'These airfoil coordinates were interpolated from the file %s to have %g points, using cosine spacing.\n',AFfile1,NumAFPoints-1);
fprintf(out_path1,'X\tY\n\n');
fprintf(out_path1,'%3.6f\t%3.6f\n',[Xnew Ynew1]');
fprintf(out_path2,'This file was generated automatically by HARP_Opt.\n');
fprintf(out_path2,'These airfoil coordinates were interpolated from the file %s to have %g points, using cosine spacing.\n',AFfile2,NumAFPoints-1);
fprintf(out_path2,'X\tY\n\n');
fprintf(out_path2,'%3.6f\t%3.6f\n',[Xnew Ynew2]');
fclose(out_path1);
fclose(out_path2);

%=========================================================================%
ThickVals = ((ThickVal1-delT):-delT:(ThickVal2+delT))';
if isempty(ThickVals)
    return; %Stop this function early, there is no need to interpolate
end
ThickVals = [ThickVal1;ThickVals;ThickVal2];

AFcoordinates(1,1) = {[Xnew,Ynew1]};
AFcoordinates(length(ThickVals),1) = {[Xnew,Ynew2]};
AFcoordinates(1,2) = {[ThickVal1]};
AFcoordinates(length(ThickVals),2) = {[ThickVal2]};
for n = 2:length(ThickVals)-1
    maxThick = 101;
    w = 1;
    while maxThick > ThickVals(n);
        Yinterp = (1-w).*Ynew2 + w.*Ynew1; %weight starts off towards y1, y1 is the thicker airfoil
        maxThick = max(Yinterp-flipud(Yinterp))*100;
        w = w - 0.00005;
    end
    MT(n,1) = maxThick;
    AFcoordinates(n,1) = {[Xnew,Yinterp]};
    AFcoordinates(n,2) = {[ThickVals(n)]};
   
    %Write the normalized airfoil coordinates each to a text file
    
    if ThickVals(n) >= 99.95;
        Thick_Suffix = '_1000';
        elseif ThickVals(n) >= 9.95 && ThickVals(n) < 99.95;
        Thick_Suffix = ['_0' num2str(10*ThickVals(n),'%3.0f')];
        else  
        Thick_Suffix = ['_00' num2str(10*ThickVals(n),'%3.0f')];
    end
    AFfile(n,:) = [family Thick_Suffix '.prof'];

fid = fopen([RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates\' AFfile(n,:)],'Wt');
fprintf(fid,'This file was generated automatically by HARP_Opt.\n');
fprintf(fid,'These airfoil coordinates represent a %g percent thick airfoil. These coordinates were interpolated between the data in the files %s and %s. %g points with cosine spacing.\n',ThickVals(n),AFfile1,AFfile2,NumAFPoints-1);
fprintf(fid,'X\tY\n\n');
fprintf(fid,'%3.6f\t%3.6f\n',[Xnew Yinterp]');
fclose(fid);

% Report current progress in the waitbar's message field
waitbar(n/(length(ThickVals)-1),H,['Interpolating between ' AFfile1 ' and ' AFfile2])
end
delete(H) % DELETE the waitbar; don't try to CLOSE it.

% plot(Xnew,Ynew1,'o-b');
% for n = 1:size(AFcoordinates,1);
% hold on; 
% plot(AFcoordinates{n,1}(:,1),AFcoordinates{n,1}(:,2),'x-r'); axis equal;
% end
% plot(Xnew,Ynew2,'o-b');
% hold off
