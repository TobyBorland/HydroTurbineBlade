function coords = interpAFprofile(afp_file1,...
                                  afp_file2,...
                                  thickVal1,...
                                  thickVal2,...
                                  FOIL_FAM, ...
                                  N_AF, ...
                                  outputDir)
                                       
shortName1 = strcat(FOIL_FAM,'_',num2str(10.*thickVal1,'%04.0f'),'.prof');
shortName2 = strcat(FOIL_FAM,'_',num2str(10.*thickVal2,'%04.0f'),'.prof');
fileStr = fullfile(outputDir,'Airfoil_Data',FOIL_FAM);

fprintf(1,['Interpolating: ' shortName1 ' to ' shortName2 '...\n'])

fid1 = fopen(afp_file1,'r');
fid2 = fopen(afp_file2,'r');                                    
xy1  = cell2mat(textscan(fid1,'%f %f','HeaderLines',4));
xy2  = cell2mat(textscan(fid2,'%f %f','HeaderLines',4));
fclose(fid1);
fclose(fid2);

x1 = xy1(:,1);
y1 = xy1(:,2);
x2 = xy2(:,1);
y2 = xy2(:,2);

[unused a1]  = min(abs(x1 - 1)); % find the index of the trailing edge
[unused a2]  = min(abs(x2 - 1)); % find the index of the trailing edge

%interpolate airfoils to have same x-coordinates
xUp  = cosspace(0,1,N_AF/2+1,'both');
xLow = flipud( xUp(1:(N_AF/2-1)) );
xNew = [xUp; xLow(1:end-1)];

yUp1  = interp1(x1(1:a1),y1(1:a1),xUp,'pchip');
yLow1 = interp1(x1(a1:end),y1(a1:end),xLow,'pchip');
yNew1 = [yUp1; yLow1(1:end-1)];

yUp2  = interp1(x2(1:a2),y2(1:a2),xUp,'pchip');
yLow2 = interp1(x2(a2:end),y2(a2:end),xLow,'pchip');
yNew2 = [yUp2; yLow2(1:end-1)];

%write the interpolated profiles for the "true" airfoils to the output directory
fid1 = fopen(fullfile(outputDir,'Airfoil_Data',shortName1),'w');
fid2 = fopen(fullfile(outputDir,'Airfoil_Data',shortName2),'w');
fprintf(fid1,'%g    N_af_nodes: number of airfoil nodes, counted clockwise starting at LE.\n',N_AF-1);
fprintf(fid1,'This file was generated automatically by HARP_Opt.\n');
fprintf(fid1,'Coordinates were interpolated from file %s to have %g points.\n',shortName1,N_AF-1);
fprintf(fid1,'Xnode\tYnode\n');
fprintf(fid1,'%3.6f\t%3.6f\n',[xNew yNew1]');
fprintf(fid2,'%g    N_af_nodes: number of airfoil nodes, counted clockwise starting at LE.\n',N_AF-1);
fprintf(fid2,'This file was generated automatically by HARP_Opt.\n');
fprintf(fid2,'Coordinates were interpolated from file %s to have %g points.\n',shortName2,N_AF-1);
fprintf(fid2,'Xnode\tYnode\n');
fprintf(fid2,'%3.6f\t%3.6f\n',[xNew yNew2]');
fclose(fid1);
fclose(fid2);

%now interpolate the profiles using weighted averaging
delT = 0.1;
thickVals = ((thickVal1-delT):-delT:(thickVal2+delT))';
if isempty(thickVals)
    %Stop this function early, there is no need to interpolate
    %because the thickness values are close enough already
    return;
end
thickVals = [thickVal1;thickVals;thickVal2];

%store coordinates and correspsonding thickness value in a cell array
coords(1,1)                = {[xNew,yNew1]};
coords(1,2)                = {thickVal1};
coords(numel(thickVals),1) = {[xNew,yNew2]};
coords(numel(thickVals),2) = {thickVal2};

for n = 2:length(thickVals)-1
    maxThick = 101;
    w = 1;
    while maxThick > thickVals(n);
        yInterp  = (1-w).*yNew2 + w.*yNew1;
        maxThick = max(yInterp-flipud(yInterp))*100;
        w        = w - 0.001;
    end
    coords(n,1) = {[xNew,yInterp]};
    coords(n,2) = {thickVals(n)};

    %Write the interpolated airfoil coordinates to a new text file
    afp_file = [fileStr '_' num2str(10.*thickVals(n),'%04.0f') '.prof'];

    fid = fopen(afp_file,'w');
    fprintf(fid,'%g    N_af_nodes: number of airfoil nodes, counted clockwise starting at LE.\n',N_AF-1);
    fprintf(fid,'This file was generated automatically by HARP_Opt.\n');
    fprintf(fid,'Coordinates represent a %g percent thick airfoil, interpolated between data in files %s and %s.\n',thickVals(n),shortName1,shortName2);
    fprintf(fid,'X\tY\n');
    fprintf(fid,'%3.6f\t%3.6f\n',[xNew yInterp]');
    fclose(fid);

   
end