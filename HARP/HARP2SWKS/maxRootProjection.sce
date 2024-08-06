//maxRootProjection

//angleI  =  0.8839453;
angleI = atan(hubIz/hubIx)
hubR = 82.5;
maxY = 1.2554217;
minY = -123.13812;  

ZfacingThird = linspace(angleI,%pi/3+%pi/2,50);
hubZTop = [hubR.*sin(ZfacingThird') hubR.*cos(ZfacingThird')];
interG = linspace(angleI+%pi/2,%pi/3+%pi/2,50);
hubITop = [hubR.*sin(interG') hubR.*cos(interG')];


H = maxY-minY;
AA = hubR*(+5*%pi/6-angleI);
alpha = atan(H/AA)
halfSector = 0.5*sqrt(H^2 + AA^2)
chordR = halfSector/cos(%pi/2-alpha);
xR1 = hubR*angleI + AA;
yR2 = maxY-chordR; // check is origin at top or base of hub
// coincides with log spiral / hub intersection
xR2 = hubR*angleI;
yR1 = chordR-H+maxY;

arcR2 = linspace(-2*alpha+%pi/2,%pi/2,50)
arcR1 = linspace(-2*alpha+3*%pi/2,3*%pi/2,50)

R1 = [chordR.*cos(arcR1)'+xR1 chordR.*sin(arcR1)'+yR1];
R2 = [chordR.*cos(arcR2)'+xR2 chordR.*sin(arcR2)'+yR2];

plot([angleI*hubR 5*%pi/6*hubR 5*%pi/6*hubR angleI*hubR angleI*hubR],[maxY maxY maxY-H maxY-H maxY])
//plot([angleI*hubR angleI+2*%pi*hubR angleI+2*%pi*hubR angleI*hubR angleI*hubR],[maxY maxY maxY-H maxY-H maxY],'g')
plot(hubZTop(:,1), hubZTop(:,2));
plot(hubITop(:,1), hubITop(:,2),'r');
plot(R1(:,1),R1(:,2),'g');
plot(R2(:,1),R2(:,2),'g');
plot(hubIz,hubIx,'r+');
plot(hubR.*sin(%pi/3+%pi/2), hubR.*cos(%pi/3+%pi/2),'g+');
plot(hubR.*sin(angleI), hubR.*cos(angleI),'c+');



rootFoil = [hubR.*cos(R1(:,1)./hubR) R1(:,2) hubR.*sin(R1(:,1)./hubR)
            hubR.*cos(R2(:,1)./hubR) R2(:,2) hubR.*sin(R2(:,1)./hubR)];
            
figure();
plot3d(rootFoil(:,1),rootFoil(:,2),rootFoil(:,3))
