// adjust Selig UIUC tables for use in HARP-Opt
clf;
HARP_dir = "C:\Documents and Settings\Foobert\Desktop\Geezer_Hydro\HARP\";
HARP_Airfoil_fname = fullfile(HARP_dir, "\Input_Files\", "\Airfoil_Data\", "new2_BW3_0050.prof");

// "new2_BW3_0050.prof" has smoothed BW3 coords

[Airfoil_coords, Prof_text] = fscanfMat(HARP_Airfoil_fname);
disp(Prof_text(1));

// *.prof format starts ax X,Y = 0.0, 0,0 
// X-values singular for pairs of Y-values
// tuple ordered in traversal from top leading edge to trailing edge, then traversal via lower surface

// Selig ordered from 
// interpolate airfoil underside values
plot(Airfoil_coords(:,1),Airfoil_coords(:,2),'k+');
a=gca();
a.isoview="on"; // isoview mode

plot(Airfoil_coords(1,1),Airfoil_coords(1,2),'ro');
plot(Airfoil_coords(41,1),Airfoil_coords(41,2),'bo');
plot(Airfoil_coords(1:40,1),Airfoil_coords(1:40,2),'r.');
plot(Airfoil_coords(41:$,1),Airfoil_coords(41:$,2),'b.');

upperSurface = Airfoil_coords(1:40,:);
upperSurface = upperSurface($:-1:1,:); // flipud()
//uSx1 = 0.99.*upperSurface(:,1);
//uSy1 = 0.99.*upperSurface(:,2);
//upperSurface1 = [[0; uSx1; 1], interp1(uSx1 , uSy1 ,[0; uSx1; 1] , "spline", "extrap")];
upperSurface1 = [[0; (0.999.*upperSurface(:,1) + 0.0003); 1] [0; (0.999.*upperSurface(:,2)); 0]];

lowerSurface = Airfoil_coords(41:$,:);
lowerSurface1 = [[0; (0.999.*lowerSurface(:,1) + 0.0003)] [0; (0.999.*lowerSurface(:,2))]];
plot(upperSurface1(:,1), upperSurface1(:,2), 'm')
plot(lowerSurface1(:,1), lowerSurface1(:,2), 'm')

lowerSurface1 = lowerSurface1($:-1:1,:); // flipud()
AC2 = [upperSurface1; lowerSurface1];
disp("--PRINT TO SCREEN--");
for i = 1:max(size(AC2))
    mprintf('%1.8f     %1.8f\n',AC2(i,1),AC2(i,2))
end


abort;
//AC2 = AC2 + 1;
plot(0.99.*AC2(:,1),0.99.*AC2(:,2),'g:');
AC3 = 0.99.*AC2(:,1);
ACext = [AC3, interp1(AC3 , 0.99.*AC2(:,2) ,[0; AC3; 1] , "spline", "extrap")]

abort;

shiftLowerSurface = [lowerSurface(1:$-1,1) + 0.001; lowerSurface($,1) - 0.001];
shiftLowerSurface = [shiftLowerSurface interp1(lowerSurface(:,1) , lowerSurface(:,2) ,shiftLowerSurface , "spline", "extrap")];
plot(shiftLowerSurface(:,1),shiftLowerSurface(:,2),'g*');
Airfoil_coords = [upperSurface; shiftLowerSurface];
