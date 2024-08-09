function R = nRad(R)
    //disp('R1');    disp(R)
    if R >= 2*%pi then
      R = R - 2*%pi;
    elseif R <= -2*%pi then
      R = R + 2*%pi;
    end
    //disp('R2');    disp(R)
endfunction

//function A3 = unduplicate(A3)
//    if size(A3,2) ~=3 then
//        error('unduplicate() requires tuple input');
//    end
//    if size(A3,1) ~=1 then
//        return;
//    end
//    for i = 2: size(A3,2)
//        if and((A3(i,1)-A3(i+1,1)<=%eps),(A3(i,2)-A3(i+1,2)<=%eps),(A3(i,3)-A3(i+1,3)<=%eps))  then
//            
//        end
//    end
//endfunction

function iA = nearest(A,n)
    // A: vector
    // n: real value
    // iA: nearest vector value
    a = min(abs(A(:)-n));
    if length(a)>1 then
        a=a(1); // somewhat arbitrary
    end
    iA = find(A(:)==(a+n));
    if iA==[] then
       iA = find(A(:)==abs(a-n));
    end
endfunction

function p_area = area_calc(x,y)
    // Get the number of vertices
    n = length(x);
    
    // Initialize the area
    p_area = 0;
    
    for i = 1 : n-1
        p_area = p_area + (x(i) + x(i+1))*(y(i) - y(i+1));
    end
    p_area = abs(p_area)/2;
endfunction

function [cx,cy] = chordcentreline(ax,ay)
    // uses airfoil plotting convention of wperf .prof
    cO = floor(length(ax)/2);
    cx = zeros(cO,1);
    cy = cx;disp('chordcentrelinedbg');
    for i = 1:cO
        mX = (ax(i) - ax($ + 1 - i))/2;
        mY = (ay(i) - ay($ + 1 - i))/2;
        cx(i) = ax(i) - mX;
        cy(i) = ay(i) - mY;
    end
    cy = interp1(cx, cy, ax,'spline');
    cx = ax;
endfunction

function [Tx,Ty,CLSx,CLSy,C3x,C3y,CLSth] = pTransform(Px,Py,a,b,hR,bR)
    // find intersection x1, y1 of log spiral R = a*exp(b*theta)
    // with circle of radius hR
    th1 = -log((a/hR)*sqrt(b^2+1))/b;
    chR = a*exp(b*th1);
    ix = chR*cos(th1);
    iy = chR*sin(th1);
    //plot(ix, iy,'y.');
    
    // centre of circle ccx, ccy
    th2 = th1 + atan(1/b) - %pi/2;
    ccx = hR*cos(th2) + ix;
    ccy = hR*sin(th2) + iy;
    //plot(ccx, ccy,'y.');
    
    // log spiral rotation correction angle thmin
    thmin = log(bR/a)/b;
    ixr = chR*cos(th1-thmin);
    iyr = chR*sin(th1-thmin);
    //plot(ixr, iyr,'y.');
    CLSth = th1-thmin;
    
    ccR = sqrt(ccx^2 + ccy^2);// != hubR!
    ccA = atan(ccy/ccx)+%pi;//_SLOPPISSIMO
    // log spiral rotation correction
    ccxr = ccR*cos(ccA-thmin);
    ccyr = ccR*sin(ccA-thmin);
    //plot(ccxr, ccyr,'y.');

    CLSx = ixr-ccxr;
    CLSy = iyr-ccyr;
    // determine jx, jy, rotated 2*%pi/3 
    // from circle spiral intersection about hub
    th3 = %pi + th2 - (2*%pi/3);
    jx = hR*cos(th3) + ccx;
    jy = hR*sin(th3) + ccy;
    //plot(jx, jy,'y.');
    // log spiral rotation correction
    jxr = hR*cos(th3-thmin)+ccxr;
    jyr = hR*sin(th3-thmin)+ccyr;
    //plot(jxr, jyr,'y.');
    
    C3x = jxr-ccxr;
    C3y = jyr-ccyr;
    //plot(C3x, C3y,'y.');
    
    Px = Px+ccxr;
    Pth = log(Px./a)./b - thmin;
    Px2 = sqrt((Px.*sin(Pth)).^2 + Px.^2);
    Pth2 = (log(Px2./a)./b)-thmin;
    Tx = Px-ccxr;
    Ty = Px2.*sin(Pth2)-ccyr;
endfunction
        
function [Tx,Ty] = pTransform2(Px,Py,a,b,hR,bR)
    // create log curve shifted to marry hub radius
    th1 = log(Px./A)./B;
    thmin = log(bR/a)/b;
   
    th2 = -log((a/hR)*sqrt(b^2+1))/b;

    ix = A*exp(B*th2)*cos(th2);
    iy = A*exp(B*th2)*sin(th2);
    
    th3 = th2 + atan(1/B) - %pi/2;
    ccx = hubR*cos(th3) + ix;
    ccy = hubR*sin(th3) + iy;
    
    ccR = sqrt(ccx^2 + ccy^2);
    ccA = atan(ccy/ccx)+%pi;
    
    ccxr = ccR*cos(ccA-thmin);
    ccyr = ccR*sin(ccA-thmin);
    
    Tx = Px.*cos(th1-thmin)-ccxr;
    Ty = Px.*sin(th1-thmin)-ccyr;
endfunction
        
function [rOverR, r, preTwist, chord, percT, t, pitchAxis] = interpBladeData(fname, numSects, interpType)
    // This program takes a tabluated blade design
    // and interpolates between it to create more airfoil sections
    
    // import data
    // this code assumes that the data is formatted in columns of
    // r/R, r, preTwist, chord, % thickness, thickness, pitch axis location
    
    [LU_xls,HO_SST,SheetN,SheetP] = xls_open(fname);
    // LU_xls: a number, the logical unit on the Excel stream.
    // HO_SST: vector of all character strings which appear in the Excel sheets.
    // SheetN: a vector of strings, the sheet names.
    // SheetP: a vector of numbers, the position of the beginning of sheets in the Excel stream.
    
    [data, iText] = xls_read(LU_xls, SheetP(1));
    // LU_xls: the logical unit on the Excel stream returned by xls_open.
    //SheetP: the position of the beginning of the sheet in the Excel stream.
    //        This position is one of those returned by xls_open.
    // data: numbers matrix, the numerical data found in the sheet. 
    //        Cells without numerical data are represented by NaN values.
    // iText: matrix of indices with the same size as Value. 
    //        The 0 indices indicates that no string exists in the corresponding Excel cell. 
    //        A positive index i points to the string SST(i) where SST is given by xls_open.
    mclose(LU_xls);
    
    // data = importdata(dataFile);
    offsetR = 7;
    offsetC = 22;
    offsetRend = 1;
    while ~isnan(data(offsetR + offsetRend, 1 + offsetC))
        offsetRend = offsetRend + 1;
    end
    offsetRend = offsetRend - 1;
    harpOpt = struct('r', data(offsetR:(offsetRend + offsetR), 2 + offsetC),...
                     'rOverR', data(offsetR:(offsetRend + offsetR), 1 + offsetC),...
                     'preTwist', data(offsetR:(offsetRend + offsetR), 3 + offsetC),...
                     'chord', data(offsetR:(offsetRend + offsetR), 4 + offsetC),...
                     'percT', data(offsetR:(offsetRend + offsetR), 5 + offsetC),...
                     't', data(offsetR:(offsetRend + offsetR), 6 + offsetC),...
                     'pitchAxis', data(offsetR:(offsetRend + offsetR), 7 + offsetC),...
                     'length', data(offsetR:(offsetRend + offsetR), 7 + offsetC));
                     
    radius = harpOpt.r(1)/harpOpt.rOverR(1); 
    harpOpt.length = radius - harpOpt.r(1);  
    //r = linspace(harpOpt.r(1), harpOpt.r($), numSects);//________________________
    r = linspace(harpOpt.r(1), radius, numSects);// practical hub radius to wingtip
    //r = linspace(0, radius, numSects); // origin to wingtip
    rOverR = r/radius;
    //pause
    
    // pitchAxis is Shell Thickness on compiled & uncompiled HARP_Opt
    
    preTwist = interp1(harpOpt.r,harpOpt.preTwist,r,interpType);
    chord = interp1(harpOpt.r,harpOpt.chord,r,interpType);
    percT = interp1(harpOpt.r,harpOpt.percT,r,interpType);
    t = interp1(harpOpt.r,harpOpt.t,r,interpType);
    
    
    pitchAxis = interp1(harpOpt.r,harpOpt.pitchAxis,r,interpType);
    lengthB = harpOpt.length;
endfunction


function createSWKswb(fname,A,rangeStart,rangeStop)
    // Write the data in a Macro that can be read by SolidWorks
    // will only read 14 airfoil sections -> split into several swb files & convert to swp
    // SWKS imports curve as metres, this script set in millimetres.
    
    if iscell(A) then
        if size(A(rangeStart).entries,1)*size(A(rangeStart).entries,2)*(rangeStop - rangeStart) > 3444 then
            error("specified range too great")
        end
    else
        if size(A,2)*(rangeStop - rangeStart) > 3444 then // assuming tuples type columns
            error("specified range too great")
        end
    end
    fid = mopen(fname,'w');
    mfprintf(fid,'%s\n','Dim swApp As Object');
    mfprintf(fid,'%s\n','Dim Part As Object');
    mfprintf(fid,'%s\n','Dim SelMgr As Object');
    mfprintf(fid,'%s\n','Dim boolstatus As Boolean');
    mfprintf(fid,'%s\n','Dim longstatus As Long, longwarnings As Long');
    mfprintf(fid,'%s\n','Dim Feature As Object');
    mfprintf(fid,'%s\n','Sub main()');
    mfprintf(fid,'\n%s\n','Set swApp = Application.SldWorks');
    
    mfprintf(fid,'\n%s\n','Set Part = swApp.ActiveDoc');
    mfprintf(fid,'%s\n','Set SelMgr = Part.SelectionManager');
    mfprintf(fid,'%s\n','swApp.ActiveDoc.ActiveView.FrameState = 1');
    
    if iscell(A) then
        for i = rangeStart:rangeStop
            temp = cell2mat(A(i));
            mfprintf(fid,'%s\n','Part.InsertCurveFileBegin');
            for j = 1:max(size(temp))//length(temp)
                mfprintf(fid,'%s','Part.InsertCurveFilePoint ');
                mfprintf(fid,'%f, %f, %f\n',temp(j,1)/1000,temp(j,2)/1000,temp(j,3)/1000);
            end
            mfprintf(fid,'%s\n','Part.InsertCurveFileEnd');
        end
    else
        mfprintf(fid,'%s\n','Part.InsertCurveFileBegin');
        for i = rangeStart:rangeStop
            mfprintf(fid,'%s','Part.InsertCurveFilePoint ');
            if (abs(A(i,1)) - 0 < 2*%eps) then
              mfprintf(fid,'%i,',0);
            else
              mfprintf(fid,'%f,',A(i,1)/1000);
            end
            if (abs(A(i,2)) - 0 < 2*%eps) then
              mfprintf(fid,' %i,',0);
            else
              mfprintf(fid,' %f,',A(i,2)/1000);
            end
            if (abs(A(i,3)) - 0 < 2*%eps) then
              mfprintf(fid,' %i\n',0);
            else
              mfprintf(fid,' %f\n',A(i,3)/1000);
            end
            
            //mfprintf(fid,'%f, %f, %f\n',A(i,1),A(i,2),A(i,3)); // bug in zero format specifier -> 0.000000 not 0
        end
        mfprintf(fid,'%s\n','Part.InsertCurveFileEnd');
    end
    
    mfprintf(fid,'%s\n','End Sub');
    mclose(fid);
endfunction

function SWKSave(sDir,swkstring)
    // generic save process
    A = eval(swkstring);
    if size(A,2)~=3 then
        error('tuple input');
    end
    saveName = fullfile(sDir, swkstring+'.sldcrv');
    csvWrite([A(:,1), A(:,2), A(:,3)],saveName,ascii(32),".","%.8f");
    saveName = fullfile(sDir, swkstring+'.swb');
    createSWKswb(saveName,A,1,size(A,1));
endfunction


// This program interpolates between blade sections generated by Harp_Opt to create a
// "blade skeleton". The xy coords that define the skeleton are saved in text files.
// These data files can be used in a SolidWorks macro to create a fully 3D blade geometry.
// A solid works macro file that should automatically create the blade geometry is also
// created.
//
// Input: Blade shape, airfoil coordinates
// Output: Coordinates the define airfoil skeleton and a Solidworks macro for creating
// the airfoil geometry

//clc; clear all; //xdel(winsid());//close all;

// Plot intermediate airfoils, y=yes, n=no
plotBlade = 'y'; // java memory leak?

// Load airfoil data and define the blade geometry
disp("--XLS OUTPUT FILES MUST BE REFORMATTED TO MS EXCEL 2003 WITHIN LIBREOFFICE--");
disp('COORD FRAME: Blade || Z-axis, Aerofoil chord || X-axis, Aerofoil thickness || Y-axis');


HARP_dir = pwd();
//dataFile =  fullfile(HARP_dir, "\Output_Files\", "\1W5_BW3_1m\", "1W5_BW3_1m_OutputXLS3.xls");
dataFile =  fullfile(HARP_dir, "\Output_Files\", "\1W5_BW3_1m_hub\", "1W5_BW3_1m_hub_OutputXLS3.xls");

disp('LOADING->');
disp(dataFile);

//saveDir =  fullfile(HARP_dir, "\Output_Files\", "\1W5_BW3_1m_GRP3\SWKScoords\");
saveDir =  fullfile(HARP_dir, "\Output_Files\", "\1W5_BW3_1m_hub\SWKScoords\");

subdiv = 30;

//r_R = 0;r = 0;preTwist = 0;chord = 0;perc_t = 0;t = 0;pitchAxis = 0;
[r_R,r,preTwist,chord,perc_t,t,pitchAxis] = interpBladeData(dataFile,subdiv,"spline");

//disp(".. interpolateBladeData() return call");

// metres to mm
t = 1000.*t;
r = 1000.*r;
chord = 1000.*chord;
      
perc_t = round(perc_t.*10)./10; // round the thickness to the nearest 0.1 %
preTwist = preTwist.*%pi/180;


// Root blade chord wider than 2*%pi hub radius / 3
// Find polynomial approximation to chord distribution, 

// hub calculation
// re-run HARP_Opt with modified hub diameter
// 
// chord = minimum hub radius * 2*(3^0.5).. tan(%pi/3)
//chordPoly = polyfit(r', chord', 6, GraphChk = 'Y', pChar = 'r')
//cPc = coeff(chordPoly);
//cPc(2) = cPc(2) - (2*(3^0.5));
//roots(cPc); // something wrong.
//
//plot(r(1:33),(2*(3^0.5)).*r(1:33),'r');//=> 117mm hub radius


minRad = 117;
hubRoot = nearest(r,minRad);
disp(".. hubRoot");disp(hubRoot);

//r = r(hubRoot:$);
//r_R = r_R(hubRoot:$);
//preTwist = preTwist(hubRoot:$);
//chord = chord(hubRoot:$);
//perc_t = perc_t(hubRoot:$);
//t = t(hubRoot:$);
//pitchAxis = pitchAxis(hubRoot:$);


HARP_Airfoil_fname = fullfile(HARP_dir, "\Input_Files\", "\Airfoil_Data\", "BW3_0050.prof");

[Airfoil_coords, Prof_text] = fscanfMat(HARP_Airfoil_fname);
disp(Prof_text(1));

disp(".. check airfoils file (redundant)");
//// thickness calculation on HARP_Opt fails with non-unique X values, i.e. at 0.0 & 1.0
//// interpolate airfoil underside values
//plot(Airfoil_coords(:,1),Airfoil_coords(:,2),'k+');
//plot(Airfoil_coords(41:$,1),Airfoil_coords(41:$,2),'bo');
//upperSurface = Airfoil_coords(1:40,:);
//lowerSurface = Airfoil_coords(41:$,:);
//shiftLowerSurface = [lowerSurface(1:$-1,1) + 0.001; lowerSurface($,1) - 0.001];
//shiftLowerSurface = [shiftLowerSurface interp1(lowerSurface(:,1) , lowerSurface(:,2) ,shiftLowerSurface , "spline")];
//plot(shiftLowerSurface(:,1),shiftLowerSurface(:,2),'g*');
//Airfoil_coords = [upperSurface; shiftLowerSurface];
//

//%Only need the profiles for the current percent thickness values
//perT = cell2mat(Normalized_AF_Coordinates(:,2));
//Scaled_AF_Coordinates = cell(NumSeg,2);
//for n = 1:length(PERCENT_THICKNESS)
//    Ind = FindInd(perT,PERCENT_THICKNESS(n)); 
//    Scaled_AF_Coordinates(n,1) = {CHORD(n).*Normalized_AF_Coordinates{Ind,1}};
//    Scaled_AF_Coordinates(n,2) = {perT(Ind)};
//end


// Find the indices of the airfoil data
// used for matching multiple airfoil profile - not relevant for constant thickness value
// or single airfoil type
disp("single airfoil profile type presumed");
//for i = 1:length(perc_t);
//    ind = find(airfoilCoords_t==perc_t(i));
//    airfoils{i} = airfoilCoords{ind,1};
//end

disp("..twist blade about chord and hub - check ");//pause

// Hydrofoil section inverted for use as turbine
Airfoil_coords = [Airfoil_coords(:,1) -Airfoil_coords(:,2)];

//airfoils = cell(size(Airfoil_coords,1),size(Airfoil_coords,2),subdiv);
airfoils = cell(subdiv);
airfoils2D = zeros((subdiv), 3);// get [min(x) max(x) z] for each aerofoil

// Scale the initial airfoil for each section along the chord length
// data structure of cells containing (X,Y,Z) tuples in columns representing airfoil chord points
for i = 1:subdiv
    airfoils(i).entries = [Airfoil_coords.*chord(i) ones(max(size(Airfoil_coords)), 1).*r(i)];
end

// Transform the airfoil points for the twisted airfoils
for i = 1:(subdiv)
    temp = cell2mat(airfoils(i));
    // no-flip X coords for blades rotating in ANTICLOCKWISE direction from the perspective o fthe advancing water
    airfoils(i).entries = [-temp(:,1)*cos(preTwist(i))-temp(:,2)*sin(preTwist(i)),...
                           -temp(:,1)*sin(preTwist(i))+temp(:,2)*cos(preTwist(i)),...
                           temp(:,3)];
end

// Adjust the airfoils so they are aligned along the pitch axis
// This the leading edge: simple design geometry and mold.
for i = 1:subdiv
    temp = cell2mat(airfoils(i));
    airfoils(i).entries = [temp(:,1) - pitchAxis(i)*chord(i) temp(:,2:3)];
    //airfoils2D(i,:) = [max(airfoils(i).entries(:,1)) min(airfoils(i).entries(:,1)) airfoils(i).entries(1,3)];
end

disp(".. sweep airfoil points by rotation around rotation axis");//pause

// Transform the airfoil points by rotation around rotation axis
// Hub rotates in ANTIclockwise direction (to prevent unscrewing) from reversed prop shaft orientation.
// Swept back blades 
// logarithmic spiral eyeball range determined from spiraltest.sce
A = 1;
B = -3.5; // airfoil profiles would have to be flipped back to front to sweep BACK in CLOCKWISE spiral
bladeR = 500;
hubR = 82.5;
hubIR = 13.75;

// hub height
temp = cell2mat(airfoils(hubRoot));
maxY = max(temp(:,2));
minY = min(temp(:,2));

// only the intersection hubIz, hubIx etc used
[Tz,Tx,hubIz,hubIx,hubIz2,hubIx2,angleI] = pTransform(r(hubRoot:$),zeros(r(hubRoot:$)),A,B,hubR,bladeR); 

// create a root structural guide projected on to curved hub surface
angleI = atan(hubIz/hubIx);
hubHeight = -minY;//maxY-minY;
// the rootFoil shape is projected on to a cylinder "nudge" smaller than the hub to allow easy trimming of the surface
nudge = 3;
// this is a cheap method to introduce curvature at the LE of the rootChord profile, drop data points at the corner and
// rely on SWKS to interpret the gap as a curved spline
slack = 10;
rootWidthRad = (hubR-nudge)*(5*%pi/6-angleI)+2*nudge;
alpha = atan(hubHeight/rootWidthRad);
halfSector = 0.5*sqrt(hubHeight^2 + rootWidthRad^2);
chordR = halfSector/cos(%pi/2-alpha);
xR1 = (hubR-nudge)*angleI + rootWidthRad-2*nudge;
yR1 = chordR-hubHeight;//+maxY;
// coincides with log spiral / hub intersection
xR2 = (hubR-nudge)*angleI-2*nudge;
yR2 = -chordR;//+maxY; // check is origin at top or base of hub

// closed rootChord has same number of points as Airfoil_coords foil, i.e. 82+smidgin
// requires identical value at each end of tuple set to close loop in SWKS
arcR2 = linspace(-2*alpha+%pi/2,%pi/2,-1+max(size(Airfoil_coords))/2+slack);
arcR1 = linspace(-2*alpha+3*%pi/2,3*%pi/2,max(size(Airfoil_coords))/2+slack);
R1 = [chordR.*cos(arcR1)'+xR1 chordR.*sin(arcR1)'+yR1];
R2 = [chordR.*cos(arcR2)'+xR2 chordR.*sin(arcR2)'+yR2];

//rootChord_1 = [hubR.*cos(R1(:,1)./hubR) R1(:,2) hubR.*sin(R1(:,1)./hubR)];
//rootChord_2 = [hubR.*cos(R2(:,1)./hubR) R2(:,2) hubR.*sin(R2(:,1)./hubR)];
//rootChordLE = [hubIx 0 hubIz];
//rootChordTE = [hubR.*sin(-%pi/3) minY hubR.*cos(-%pi/3)];
rootChord_o1 = [(hubR-nudge).*cos(R1(:,1)./(hubR-nudge)) R1(:,2) (hubR-nudge).*sin(R1(:,1)./(hubR-nudge))];
rootChord_o2 = [(hubR-nudge).*cos(R2(:,1)./(hubR-nudge)) R2(:,2) (hubR-nudge).*sin(R2(:,1)./(hubR-nudge))];
//rootChord = [rootChordLE; rootChord_o1; rootChordTE; rootChord_o2; rootChordLE];
rootChord = [rootChord_o1(slack:$,:); rootChord_o2(2:$-slack,:); rootChord_o1(slack,:)];

// place hydrofoil outlines along Z axis
Rz = r;//[1e-10 r(2:$)-r(1)];// r(2:hubRoot-1)];
for i = 1:length(Rz)
    [Zdiff,Xdiff] = pTransform2(Rz(i),0,A,B,hubR,bladeR); 
    Tz(i) = Rz(i)-Zdiff;
    Tx(i) = Xdiff;

    while abs(Rz(i)-Zdiff) > 1e-10
        [Zdiff,Xdiff] = pTransform2(Tz(i),0,A,B,hubR,bladeR); 
        Tz(i) = Tz(i)+(Rz(i)-Zdiff);//disp(Rz(i)-Zdiff);pause
        Tx(i) = Xdiff;
    end
end

// create leading edge curve
LEspline_leader = [1e-10 r(2:$)-r(1)];
[LEsplineZ_leader,LEsplineX_leader] = pTransform2(LEspline_leader,zeros(LEspline_leader),A,B,hubR,bladeR);
LEsplineX_leader = LEsplineX_leader(find((LEsplineZ_leader > hubIz)&(LEsplineZ_leader < r(hubRoot))));
LEsplineZ_leader = LEsplineZ_leader(find((LEsplineZ_leader > hubIz)&(LEsplineZ_leader < r(hubRoot))));

//LEsplineL = [rootChord_o1(slack,:)
//             hubIx 0 hubIz
//             LEsplineX_leader' 0.*ones(LEsplineX_leader') LEsplineZ_leader'
//             Tx(hubRoot) 0 Rz(hubRoot)];
             
slack2 = slack+1;//_________________________________________________________
LEsplineL = [rootChord_o1(slack,:)
             LEsplineX_leader(:,slack2:$)' 0.*ones(LEsplineX_leader(:,slack2:$)') LEsplineZ_leader(:,slack2:$)'
             Tx(hubRoot) 0 Rz(hubRoot)];
             
// do not discard wingtip in createBladeSkeleton8             
LEspline = [Tx(hubRoot:$)' 0.*ones(Tx(hubRoot:$)') Rz(hubRoot:$)']; 

for i = 1:length(Rz)//hubRoot-1
    temp = cell2mat(airfoils(i));
    airfoils(i).entries = [temp(:,1) + Tx(i), temp(:,2), temp(:,3)];
end

// Trailing Edge at consistent TE point
TEspline = zeros(subdiv,3);
topSpline = zeros(subdiv,3);
baseSpline = zeros(subdiv,3);
for i = 1:subdiv
    temp = cell2mat(airfoils(i));
    TEspline(i,:) = [temp(find(temp(:,1)==min(temp(:,1)),1),1) temp(find(temp(:,1)==min(temp(:,1)),1),2) temp(1,3)];
    topSpline(i,:) = [temp(63,1) temp(63,2) temp(63,3)];
    baseSpline(i,:) = [temp(21,1) temp(21,2) temp(21,3)];
end

TEsplineL = [rootChord_o2(1,:); hubR.*sin(-%pi/3) minY hubR.*cos(-%pi/3); TEspline(hubRoot,:)];
//TEsplineL2 = [hubR.*sin(-%pi/3) minY hubR.*cos(-%pi/3); TEspline(hubRoot:$,:)];
TEspline = TEspline(hubRoot:$,:);
topSpline = [rootChord(63,:); topSpline(hubRoot:$,:)];
baseSpline = [rootChord(21,:); baseSpline(hubRoot:$,:)];

//temp = cell2mat(airfoils(hubRoot)); //wingtip
////split hubRoot into upper and lower splines
//midChord = floor(size(temp,1)/2);
//hubRootFoil1 = temp(1:midChord,:);
//hubRootFoil2 = temp(midChord:$,:);

// create hub circles, 2 concentric bolt circles
twopi = linspace(0,2*%pi,3*subdiv);
ZfacingThird = linspace(-%pi/3+%pi/2,%pi/3+%pi/2,subdiv);
ZfacingSixth = linspace(-%pi/12+%pi/2,%pi/12+%pi/2,subdiv);
XfacingHalf = linspace(-%pi/2,%pi/2,subdiv);

hub3RDZTop = [hubR.*cos(ZfacingThird') maxY.*ones(subdiv,1) hubR.*sin(ZfacingThird') ];
hub3RDZBase = [hubR.*cos(ZfacingThird') minY.*ones(subdiv,1) hubR.*sin(ZfacingThird')];
hubInner3RDZTop = [hubIR.*cos(ZfacingThird') maxY.*ones(subdiv,1) hubIR.*sin(ZfacingThird') ];
hubInner3RDZBase = [hubIR.*cos(ZfacingThird') minY.*ones(subdiv,1) hubIR.*sin(ZfacingThird')];

// fractional overlap at wingtip
bladeTop = [bladeR.*cos(ZfacingSixth') (maxY+20).*ones(subdiv,1) bladeR.*sin(ZfacingSixth') ];
//bladeTop = bladeTop(1:$-1,:);
bladeBase = [bladeR.*cos(ZfacingSixth') (minY-20).*ones(subdiv,1) bladeR.*sin(ZfacingSixth')];
//bladeBase = bladeBase(1:$-1,:);
hubTop = [hubR.*cos(twopi') maxY.*ones(3*subdiv,1) hubR.*sin(twopi') ];
hubBase = [hubR.*cos(twopi') minY.*ones(3*subdiv,1) hubR.*sin(twopi')];
hubInnerTop = [hubIR.*cos(twopi') maxY.*ones(3*subdiv,1) hubIR.*sin(twopi') ];
hubInnerBase = [hubIR.*cos(twopi') minY.*ones(3*subdiv,1) hubIR.*sin(twopi')];

hub3RDXTop = [hubR.*cos(XfacingHalf') maxY.*ones(subdiv,1) hubR.*sin(XfacingHalf') ];
hub3RDXBase = [hubR.*cos(XfacingHalf') minY.*ones(subdiv,1) hubR.*sin(XfacingHalf')];
hubInner3RDXTop = [hubIR.*cos(XfacingHalf') maxY.*ones(subdiv,1) hubIR.*sin(XfacingHalf') ];
hubInner3RDXBase = [hubIR.*cos(XfacingHalf') minY.*ones(subdiv,1) hubIR.*sin(XfacingHalf')];

hubSegX1U = [hubIx maxY hubIz; 0 maxY 0];
hubSegX2U = [0 maxY 0; hubIx2 maxY hubIz2];

hubSegX1L = [hubIx minY hubIz; 0 minY 0];
hubSegX2L = [0 minY 0; hubIx2 minY hubIz2];

hubSegZ1U = [hubR.*sin(%pi/3) maxY hubR.*cos(%pi/3)
             hubIR.*sin(%pi/3) maxY hubIR.*cos(%pi/3)];
hubSegZ2U = [hubIR.*sin(-%pi/3) maxY hubIR.*cos(-%pi/3)
             hubR.*sin(-%pi/3) maxY hubR.*cos(-%pi/3)];

hubSegZ1L = [hubR.*sin(%pi/3) minY hubR.*cos(%pi/3)
             hubIR.*sin(%pi/3) minY hubIR.*cos(%pi/3)];
hubSegZ2L = [hubR.*sin(-%pi/3) minY hubR.*cos(-%pi/3)
             hubIR.*sin(-%pi/3) minY hubIR.*cos(-%pi/3)];
                   
boltcircleR = 37.5;
boltcircleTop = [boltcircleR.*cos(twopi') maxY.*ones(3*subdiv,1) boltcircleR.*sin(twopi')];
//boltcircleBase = [boltcircleR.*cos(twopi') minY.*ones(3*subdiv,1) boltcircleR.*sin(twopi')];

boltcircle2R = 66.5;
//boltcircle2Top = [boltcircle2R.*cos(twopi') maxY.*ones(3*subdiv,1) boltcircle2R.*sin(twopi')];
//boltcircle2Base = [boltcircle2R.*cos(twopi') minY.*ones(3*subdiv,1) boltcircle2R.*sin(twopi')];

bolt1U = [boltcircleR*cos(2*%pi/3) maxY boltcircleR*sin(2*%pi/3); (boltcircleR+7)*cos(2*%pi/3) maxY (boltcircleR+7)*sin(2*%pi/3)];
bolt2U = [boltcircle2R*cos(2*%pi/3) maxY boltcircle2R*sin(2*%pi/3); (boltcircle2R+7)*cos(2*%pi/3) maxY (boltcircle2R+7)*sin(2*%pi/3)];
bolt3U = [boltcircleR*cos(%pi/3) maxY boltcircleR*sin(%pi/3); (boltcircleR+7)*cos(%pi/3) maxY (boltcircleR+7)*sin(%pi/3)];
bolt4U = [boltcircle2R*cos(%pi/3) maxY boltcircle2R*sin(%pi/3); (boltcircle2R+7)*cos(%pi/3) maxY (boltcircle2R+7)*sin(%pi/3)];

bolt1L = [boltcircleR*cos(2*%pi/3) minY boltcircleR*sin(2*%pi/3); (boltcircleR+7)*cos(2*%pi/3) minY (boltcircleR+7)*sin(2*%pi/3)];
bolt2L = [boltcircle2R*cos(2*%pi/3) minY boltcircle2R*sin(2*%pi/3); (boltcircle2R+7)*cos(2*%pi/3) minY (boltcircle2R+7)*sin(2*%pi/3)];
bolt3L = [boltcircleR*cos(%pi/3) minY boltcircleR*sin(%pi/3); (boltcircleR+7)*cos(%pi/3) minY (boltcircleR+7)*sin(%pi/3)];
bolt4L = [boltcircle2R*cos(%pi/3) minY boltcircle2R*sin(%pi/3); (boltcircle2R+7)*cos(%pi/3) minY (boltcircle2R+7)*sin(%pi/3)];

// some values/guidelines for construction an enclosing mold-box
LEmoldboxEdge = [LEspline(14,1) minY LEspline(14,3); LEspline(14,1) maxY LEspline(14,3)];

temp = cell2mat(airfoils(subdiv)); //wingtip
wingtipZmin = bladeR.*cos(asin((min(temp(:,1)))./bladeR));
wingtipXmin = bladeR.*sin(asin((min(temp(:,1)))./bladeR));
TEmoldboxEdge1 = [wingtipXmin minY wingtipZmin; wingtipXmin maxY wingtipZmin];
temp = cell2mat(airfoils(hubRoot));
TEmoldboxEdge2 = [min(temp(:,1)) minY temp(find(temp(:,1)==min(temp(:,1)),1),3)
                  min(temp(:,1)) maxY temp(find(temp(:,1)==min(temp(:,1)),1),3)];


disp("..plot the blade in 3d");//pause
// Plot the blade in 3d
if plotBlade == 'y';
    fig=figure; 
    //clf()
    a = newaxes(); 
    a.data_bounds=[-200,-200,-100;150,50,510];

    a.rotation_angles=[90 90];
    a.isoview="on"; // isoview mode
    a.axes_visible="off";
    for i = hubRoot:subdiv
        temp = cell2mat(airfoils(i));
        plot3d(temp(:,1),temp(:,2),temp(:,3),flag=[0,4,0]);
    end
    
    plot3d(LEsplineL(:,1),LEsplineL(:,2),LEsplineL(:,3),flag=[0,4,0])
    plot3d(LEspline(:,1),LEspline(:,2),LEspline(:,3),flag=[0,4,0])
    plot3d(TEsplineL(:,1),TEsplineL(:,2),TEsplineL(:,3),flag=[0,4,0])
    plot3d(TEspline(:,1),TEspline(:,2),TEspline(:,3),flag=[0,4,0])
    plot3d(topSpline(:,1),topSpline(:,2),topSpline(:,3),flag=[0,4,0])
    plot3d(baseSpline(:,1),baseSpline(:,2),baseSpline(:,3),flag=[0,4,0])
       
    plot3d(hub3RDZTop(:,1),hub3RDZTop(:,2),hub3RDZTop(:,3),flag=[0,4,0])
    plot3d(hub3RDZBase(:,1),hub3RDZBase(:,2),hub3RDZBase(:,3),flag=[0,4,0])
    plot3d(hubInner3RDZTop(:,1),hubInner3RDZTop(:,2),hubInner3RDZTop(:,3),flag=[0,4,0])
    plot3d(hubInner3RDZBase(:,1),hubInner3RDZBase(:,2),hubInner3RDZBase(:,3),flag=[0,4,0])
    
    //plot3d(rootChord1(:,1),rootChord1(:,2),rootChord1(:,3),flag=[0,4,0])
    plot3d(rootChord(:,1),rootChord(:,2),rootChord(:,3),flag=[0,4,0])
    //plot3d(rootChord2(:,1),rootChord2(:,2),rootChord2(:,3),flag=[0,4,0])

    plot3d(hubTop(:,1),hubTop(:,2),hubTop(:,3),flag=[0,4,0])
    plot3d(hubBase(:,1),hubBase(:,2),hubBase(:,3),flag=[0,4,0])
    plot3d(bladeTop(:,1),bladeTop(:,2),bladeTop(:,3),flag=[0,4,0])
    plot3d(bladeBase(:,1),bladeBase(:,2),bladeBase(:,3),flag=[0,4,0])
    //plot3d(hubInnerTop(:,1),hubInnerTop(:,2),hubInnerTop(:,3),flag=[0,4,0])
    //plot3d(hubInnerBase(:,1),hubInnerBase(:,2),hubInnerBase(:,3),flag=[0,4,0])

    plot3d(hubSegX1U(:,1),hubSegX1U(:,2),hubSegX1U(:,3),flag=[0,4,0])
    plot3d(hubSegX1L(:,1),hubSegX1L(:,2),hubSegX1L(:,3),flag=[0,4,0])
    plot3d(hubSegX2U(:,1),hubSegX2U(:,2),hubSegX2U(:,3),flag=[0,4,0])
    plot3d(hubSegX2L(:,1),hubSegX2L(:,2),hubSegX2L(:,3),flag=[0,4,0])
    plot3d(hubSegZ1U(:,1),hubSegZ1U(:,2),hubSegZ1U(:,3),flag=[0,4,0])
    plot3d(hubSegZ1L(:,1),hubSegZ1L(:,2),hubSegZ1L(:,3),flag=[0,4,0])
    plot3d(hubSegZ2U(:,1),hubSegZ2U(:,2),hubSegZ2U(:,3),flag=[0,4,0])
    plot3d(hubSegZ2L(:,1),hubSegZ2L(:,2),hubSegZ2L(:,3),flag=[0,4,0])
    //plot3d(boltcircleTop(:,1),boltcircleTop(:,2),boltcircleTop(:,3),flag=[0,4,0])
    //plot3d(boltcircleBase(:,1),boltcircleBase(:,2),boltcircleBase(:,3),flag=[0,4,0])
    //plot3d(boltcircle2Top(:,1),boltcircle2Top(:,2),boltcircle2Top(:,3),flag=[0,4,0])
    //plot3d(boltcircle2Base(:,1),boltcircle2Base(:,2),boltcircle2Base(:,3),flag=[0,4,0])
    plot3d(bolt1U(:,1),bolt1U(:,2),bolt1U(:,3),flag=[0,4,0])
    plot3d(bolt2U(:,1),bolt2U(:,2),bolt2U(:,3),flag=[0,4,0])
    plot3d(bolt3U(:,1),bolt3U(:,2),bolt3U(:,3),flag=[0,4,0])

    plot3d(LEmoldboxEdge(:,1),LEmoldboxEdge(:,2),LEmoldboxEdge(:,3),flag=[0,4,0])
    plot3d(TEmoldboxEdge1(:,1),TEmoldboxEdge1(:,2),TEmoldboxEdge1(:,3),flag=[0,4,0])
    plot3d(TEmoldboxEdge2(:,1),TEmoldboxEdge2(:,2),TEmoldboxEdge2(:,3),flag=[0,4,0])
    //a=gca();
end

disp('SAVING TO ->');
disp(saveDir);
pause

SWKSave(saveDir,'hubTop');
SWKSave(saveDir,'hubBase');
SWKSave(saveDir,'bladeTop');
SWKSave(saveDir,'bladeBase');
SWKSave(saveDir,'hubInnerTop');
SWKSave(saveDir,'hubInnerBase');

SWKSave(saveDir,'hub3RDZTop');
SWKSave(saveDir,'hub3RDZBase');
SWKSave(saveDir,'hubInner3RDZTop');
SWKSave(saveDir,'hubInner3RDZBase');
SWKSave(saveDir,'hub3RDXTop');
SWKSave(saveDir,'hub3RDXBase');
SWKSave(saveDir,'hubInner3RDXTop');
SWKSave(saveDir,'hubInner3RDXBase');

SWKSave(saveDir,'rootChord');
//SWKSave(saveDir,'rootChord1');
//SWKSave(saveDir,'rootChord2');

SWKSave(saveDir,'boltcircleTop'); // visual reference
//SWKSave(saveDir,'boltcircleBase');
//SWKSave(saveDir,'boltcircle2Top');
//SWKSave(saveDir,'boltcircle2Base');

SWKSave(saveDir,'hubSegX1U');
SWKSave(saveDir,'hubSegX1L');
SWKSave(saveDir,'hubSegX2U');
SWKSave(saveDir,'hubSegX2L');

SWKSave(saveDir,'hubSegZ2U');
SWKSave(saveDir,'hubSegZ2L');
SWKSave(saveDir,'hubSegZ1U');
SWKSave(saveDir,'hubSegZ1L');

SWKSave(saveDir,'LEsplineL');
SWKSave(saveDir,'LEspline');
SWKSave(saveDir,'TEsplineL');
//SWKSave(saveDir,'TEsplineL2');
SWKSave(saveDir,'TEspline');
SWKSave(saveDir,'topSpline');
SWKSave(saveDir,'baseSpline');

SWKSave(saveDir,'bolt1U');
SWKSave(saveDir,'bolt2U');
SWKSave(saveDir,'bolt3U');
SWKSave(saveDir,'bolt4U');

SWKSave(saveDir,'bolt1L');
SWKSave(saveDir,'bolt2L');
SWKSave(saveDir,'bolt3L');
SWKSave(saveDir,'bolt4L');

SWKSave(saveDir,'LEmoldboxEdge');
SWKSave(saveDir,'TEmoldboxEdge1');
SWKSave(saveDir,'TEmoldboxEdge2');


// save hub circles, bolt circles
//saveName = fullfile(saveDir, 'hubTop.sldcrv');
//csvWrite([hubTop(:,1), hubTop(:,2), hubTop(:,3)],saveName,ascii(32),".","%.8f");
//saveName = fullfile(saveDir, 'hubTop.swb');
//createSWKswb(saveName,hubTop,1,size(hubTop,1));

//// Export the airfoil points
//for i = 1:(subdiv - hubRoot + 1)//length(chord)
//    // Solidworks rounds to 2 decimal places
//    temp = 1000.*cell2mat(airfoils(i));
//    saveName = fullfile(saveDir, 'airfoil_'+string(i)+'.sldcrv');
//    csvWrite([temp(:,1), temp(:,2), temp(:,3)],saveName,ascii(32),".","%.8f");//ascii(9));
//end
//

//for i = 1:floor((subdiv - hubRoot + 1)/14)
//    swb_pathname = fullfile(saveDir, 'AF_'+string(i)+'.swb');
//    createSWKswb(swb_pathname,airfoils,((i-1)*14)+1,(i*14));
//end
//if modulo((subdiv - hubRoot + 1),14) > 0 then
//    swb_pathname = fullfile(saveDir, 'AF_'+string(i+1)+'.swb');
//    createSWKswb(swb_pathname,airfoils,i*14,(subdiv - hubRoot + 1));    
//end



for i = 1:floor(subdiv/14) 
    swb_pathname = fullfile(saveDir, 'AF_'+string(i)+'.swb');
    createSWKswb(swb_pathname,airfoils,((i-1)*14)+1,(i*14));
end
if modulo(subdiv,14) > 0 then
    swb_pathname = fullfile(saveDir, 'AF_'+string(i+1)+'.swb');
    createSWKswb(swb_pathname,airfoils,i*14,(subdiv));    
end

