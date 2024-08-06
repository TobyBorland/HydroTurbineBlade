function theta = equalCurveRadius(AA,BB,R)
    theta = (-1/BB)*log((AA/R)*sqrt(BB^2 + 1));
endfunction

function R = nRad(R)
    //disp('R1');    disp(R)
    if R >= 2*%pi then
      R = R - 2*%pi;
    elseif R <= -2*%pi then
      R = R + 2*%pi;
    end
    //disp('R2');    disp(R)
endfunction

function at2 = atan2(y,x)
    if x ~= 0 then 
        arctan = atan(y/x);
    end
 
    if x > 0 then 
        at2 = arctan;
    elseif y>=0 & x<0 then
        at2 = %pi+arctan;
    elseif y<0 & x<0 then
        at2 = arctan-%pi;
    elseif y>0 & x==0 then
        at2 = %pi/2;
    elseif y<0 & x==0 then
        at2 = %pi/-2;
    end
endfunction

airfoils2D  =  [249.21737  -0.6860396    126.33333  
                247.1655   -0.6521178    130.04377  
                245.19602  -0.6212149    133.75421  
                243.30488  -0.5929565    137.46465  
                241.48448  -0.5669916    141.17508  
                239.72414  -0.5429905    144.88552  
                238.01052  -0.5206429    148.59596  
                236.32805  -0.4996559    152.3064   
                234.66399  -0.4797999    156.01684  
                233.01605  -0.4609698    159.72727  
                231.38301  -0.4430816    163.43771  
                229.76326  -0.4260557    167.14815  
                228.15538  -0.4098222    170.85859  
                226.5579   -0.3943174    174.56902  
                224.96903  -0.3794800    178.27946  
                223.38722  -0.3652537    181.9899   
                221.81145  -0.3515874    185.70034  
                220.24055  -0.3384322    189.41077  
                218.67337  -0.3257447    193.12121  
                217.10996  -0.3135001    196.83165  
                215.5507   -0.3016797    200.54209  
                213.99588  -0.2902654    204.25253  
                212.44571  -0.2792383    207.96296  
                210.90035  -0.2685788    211.6734   
                209.35987  -0.2582681    215.38384  
                207.82436  -0.2482886    219.09428  
                206.2941   -0.2386256    222.80471  
                204.76939  -0.2292654    226.51515  
                203.25047  -0.2201945    230.22559  
                201.73764  -0.2114007    233.93603  
                200.23123  -0.2028722    237.64646  
                198.73149  -0.1945979    241.3569   
                197.23874  -0.1865672    245.06734  
                195.7533   -0.1787705    248.77778  
                194.27552  -0.1711987    252.48822  
                192.8057   -0.1638427    256.19865  
                191.34414  -0.1566945    259.90909  
                189.89113  -0.1497460    263.61953  
                188.44695  -0.1429895    267.32997  
                187.01194  -0.1364182    271.0404   
                185.58649  -0.1300257    274.75084  
                184.171    -0.1238060    278.46128  
                182.76576  -0.1177530    282.17172  
                181.3707   -0.1118596    285.88215  
                179.98566  -0.1061187    289.59259  
                178.61053  -0.1005238    293.30303  
                177.24626  -0.0950727    297.01347  
                175.89461  -0.0897658    300.72391  
                174.55733  -0.0846036    304.43434  
                173.23576  -0.0795841    308.14478  
                171.93027  -0.0747012    311.85522  
                170.64178  -0.0699478    315.56566  
                169.37077  -0.0653175    319.27609  
                168.11574  -0.0608034    322.98653  
                166.87598  -0.0563992    326.69697  
                165.65075  -0.0520984    330.40741  
                164.43973  -0.0478964    334.11785  
                163.24307  -0.0437900    337.82828  
                162.06102  -0.0397763    341.53872  
                160.89381  -0.0358528    345.24916  
                159.74193  -0.0320173    348.9596   
                158.60588  -0.0282679    352.67003  
                157.48627  -0.0246028    356.38047  
                156.38533  -0.0223585    360.09091  
                155.30668  -0.0211532    363.80135  
                154.25403  -0.0199763    367.51178  
                153.22927  -0.0188271    371.22222  
                152.22955  -0.0177030    374.93266  
                151.2512   -0.0166014    378.6431   
                150.29069  -0.0155198    382.35354  
                149.34597  -0.0144568    386.06397  
                148.41578  -0.0134119    389.77441  
                147.49894  -0.0123846    393.48485  
                146.59753  -0.0113748    397.19529  
                145.71908  -0.0103836    400.90572  
                144.87161  -0.0094118    404.61616  
                144.06159  -0.0084595    408.3266   
                143.28672  -0.0075249    412.03704  
                142.54172  -0.0066057    415.74747  
                141.82135  -0.0056995    419.45791  
                141.12159  -0.0048052    423.16835  
                140.43977  -0.0039222    426.87879  
                139.77324  -0.0030499    430.58923  
                139.12183  -0.0021882    434.29966  
                138.49332  -0.0013371    438.0101   
                137.8971   -0.0004963    441.72054  
                137.34194    0.0           445.43098  
                136.82836    0.0           449.14141  
                136.35145    0.0           452.85185  
                135.90616    0.0           456.56229  
                135.48764    0.0           460.27273  
                135.09142    0.0           463.98316  
                134.71309    0.0           467.6936   
                134.34958    0.0           471.40404  
                134.0071     0.0           475.11448  
                133.69569    0.0           478.82492  
                133.42543    0.0           482.53535  
                133.20639    0.0           486.24579  
                133.04865    0.0           489.95623  
                132.9623     0.0           493.66667];  
airfoils2D(:,1) = airfoils2D(:,1)-airfoils2D(:,2);
airfoils2D(:,2) = airfoils2D(:,2)-airfoils2D(:,2);

//figure;
xy=[linspace(1, 500, 500)' zeros(500,1)];
c = ['r','y','k','m','g','c','b'];
scf();
clf();
a = gca();
a.isoview="on"; // isoview mode

plot(xy(:,1), xy(:,2),'k');

twopi = linspace(0,2*%pi,300);// interp on airfoil set to 100
hubR = 82.5;
hub = [hubR.*cos(twopi') hubR.*sin(twopi') ];
plot(hub(:,1),hub(:,2),'k');

//for A = (2:60)
A =1;
    //for B = (0.014:0.002:0.016)
    B = -3.5
    
        th = log(xy(:,1)./A)./B;
        //th = th - min(th);
        
        lx = xy(:,1).*cos(th);// - xy(:,1).*sin(th); // 
        ly = xy(:,1).*sin(th);// + xy(:,1).*cos(th);
        plot(lx,ly,'r');
        th1 = th - min(th);
        lx = xy(:,1).*cos(th1);// - xy(:,1).*sin(th); // 
        ly = xy(:,1).*sin(th1);// + xy(:,1).*cos(th);
        plot(lx,ly,'g');
        Px = 50;Py = 0;
        plot(Px,Py,'ko');
        Pth = log(250/A)/B - min(th);
        plot(Px*cos(Pth), Px*sin(Pth),'go');
        plot(Px.*cos(twopi'), Px.*sin(twopi'),'k');
        
//        uR = xy(:,2).*B./sqrt(B^2 + 1);
//        the = th - log(B)/B - %pi/2;
//        hx = uR.*cos(the);//.*exp(B.*the);
//        hy = uR.*sin(the);//.*exp(B.*the);
//        plot(hx,hy,'r');
        
        //plot(hub(:,1)+hx,hub(:,2)+hy,'c');


//        th3 = log(xy(:,2)./A*B)./B;
//        th3 = th3.*(%pi/180);
//        th4 = log(xy(:,2)./-A*B)./B;
//        th4 = th4.*(%pi/180);
//        lx = xy(:,2).*cos(th3)-xy(:,1).*sin(th3);
//        ly = xy(:,2).*sin(th3)+xy(:,1).*cos(th3)
//        plot(lx,ly,c(modulo(A+3,7)+1));
        
        th2 = equalCurveRadius(A,B,hubR);//disp(th2);
        ix = A*exp(B*th2)*cos(th2);
        iy = A*exp(B*th2)*sin(th2);
        plot(ix, iy,'r+');
        ixr = A*exp(B*th2)*cos(th2-min(th));
        iyr = A*exp(B*th2)*sin(th2-min(th));
        plot(ixr, iyr,'g+');
        
        th3 = nRad(th2 + atan(1/B) - %pi/2);//disp(th3);
        ccx = hubR*cos(th3) + ix;
        ccy = hubR*sin(th3) + iy;
        plot(ccx, ccy,'ro');
        
        ccR = sqrt(ccx^2 + ccy^2);// != hubR!
        ccA = atan(ccy/ccx)+%pi;disp(atan(ccy/ccx));//_SLOPPISSIMO
        //plot(250*cos(ccA),250*sin(ccA),'ko');
        ccxr = ccR*cos(ccA-min(th));
        ccyr = ccR*sin(ccA-min(th));
        plot(ccxr, ccyr,'go');

        plot(hub(:,1)+ccx,hub(:,2)+ccy,'r');
        plot(hub(:,1)+ccxr,hub(:,2)+ccyr,'g');

        plot([ix ccx],[iy ccy],'r');
        plot([ixr ccxr],[iyr ccyr],'g');
        th4 = nRad(%pi + th3 - (2*%pi/3));//disp(th4);
        jx = hubR*cos(th4) + ccx;
        jy = hubR*sin(th4) + ccy;
        plot(jx, jy,'r+');
        plot([jx ccx],[jy ccy],'r');

        //iR = sqrt(ix^2 + iy^2);// != hubR!
        //iA = atan(iy/ix);//+%pi;disp(atan(ccy/ccx));//_SLOPPISSIMO

        jxr = hubR*cos(th4-min(th)) + ccxr;
        jyr = hubR*sin(th4-min(th)) + ccyr;
        plot(jxr, jyr,'g+');
        plot([jxr ccxr],[jyr ccyr],'g');
        
        plot(lx - ccxr,ly - ccyr,'c:');//_______________________________________________________________
        plot(ixr-ccxr, iyr-ccyr,'c+');
        plot([ixr-ccxr 0],[iyr-ccyr 0],'c');
        plot(jxr-ccxr, jyr-ccyr,'c+');
        plot([jxr-ccxr 0],[jyr-ccyr 0],'c');
        plot(Px*cos(Pth)-ccxr, Px*sin(Pth)-ccyr,'co');

//        Px = 250+ccxr;Py = 0;
//        plot(Px,Py,'k*');
//        Pth = log((250+ccxr)/A)/B - min(th);
//        //Qth = log((250)/A)/B - min(th);
//        //PQy = Px*sin(Qth) - Px*sin(Pth);
//        plot(Px*cos(Pth), Px*sin(Pth),'g*');
//        plot(Px*cos(Pth), Px*sin(Pth),'g*');
//        plot(Px-ccxr, Px*sin(Pth)-ccyr,'c*');
//        plot([250 250], [150 -150],'c');

        Px = 0+ccxr;Py = 0;
        plot(Px,Py,'ko');
        Pth = log(Px/A)/B - min(th);
        Ppx = sqrt((Px*sin(Pth))^2 + Px^2);
        plot(Ppx,Py,'ko');
        Ppth = log(Ppx/A)/B - min(th);
        plot(Px, Ppx*sin(Ppth) ,'r+');
        plot(Px-ccxr, Ppx*sin(Ppth)-ccyr,'r+');
        
        
        function [Tx,Ty,CLSx,CLSy,C3x,C3y] = pTransform(Px,Py,a,b,hR,bR)


//            // bR, blade radius correction for circle translate ccxr 
//            // most of these calculations are repeated below
//            th1 = -log((a/hR)*sqrt(b^2+1))/b;
//            chR = a*exp(b*th1);
//            th2 = th1 + atan(1/b) - %pi/2;
//            ix = chR*cos(th1);
//            iy = chR*sin(th1);
//            ccx = hR*cos(th2) + ix;
//            ccy = hR*sin(th2) + iy;
//            ccR = sqrt(ccx^2 + ccy^2);// != hubR!
//            ccA = atan(ccy/ccx)+%pi;//_SLOPPISSIMO
//            thmin = log(bR/a)/b;
//            ccxr = ccR*cos(ccA-thmin);
//            bR = bR+ccxr;

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
            
            ccxr = ccR*cos(ccA-min(th));
            ccyr = ccR*sin(ccA-min(th));
            
            Tx = Px.*cos(th1-thmin)-ccxr;
            Ty = Px.*sin(th1-thmin)-ccyr;
        endfunction
        
        testX = 350;testY = 0;
        plot(testX,testY,'mo')
        [pTestX, pTestY,clsX,clsY,c3x,c3y] = pTransform(350,0,A,B,hubR,500);//+61.3383);
        plot(pTestX,pTestY,'mo')
        plot(clsX,clsY,'m*')
        plot(c3x,c3y,'m*')
        
    [pX, pY,clsX,clsY,c3x,c3y] = pTransform(airfoils2D(:,3),zeros(size(max(airfoils2D))),A,B,hubR,500);
for i = 1:max(size(airfoils2D)) 
    plot([pX(i) pX(i)],[pY(i) pY(i)-airfoils2D(i,1)],'b');
end
[pX, pY,clsX,clsY,c3x,c3y] = pTransform(linspace(1,200,50),zeros(1,50),A,B,hubR,500);
//[pX, pY] = pTransform2(linspace(1,200,50),zeros(1,50),A,B,hubR,500);
//[pX, pY] = pTransform2(linspace(0.01,hubR,50),zeros(1,50),A,B,hubR,500);
plot(pX, pY,'r:');
//        theta = theta + (2*%pi/3);
//        lx = xy(:,1).*cos(theta)-xy(:,2).*sin(theta);
//        ly = xy(:,1).*sin(theta)+xy(:,2).*cos(theta)
//        plot(lx,ly,c(modulo(A+B,7)+1));
//        theta = theta + (2*%pi/3);
//        lx = xy(:,1).*cos(theta)-xy(:,2).*sin(theta);
//        ly = xy(:,1).*sin(theta)+xy(:,2).*cos(theta)
//        plot(lx,ly,c(modulo(A+B,7)+1));
    //end
//end
                         
//A =20;
//    //for B = (0.014:0.002:0.016)
//    B = 0.02
//        theta = log(xy(:,2)./A)./B;
//        theta = theta.*(%pi/180);
//        lx = xy(:,1).*cos(theta)-xy(:,2).*sin(theta);
//        ly = xy(:,1).*sin(theta)+xy(:,2).*cos(theta)
//        plot(lx,ly,c(modulo(A+B+1,7)+1));
//        theta = theta + (2*%pi/3);
//        lx = xy(:,1).*cos(theta)-xy(:,2).*sin(theta);
//        ly = xy(:,1).*sin(theta)+xy(:,2).*cos(theta)
//        plot(lx,ly,c(modulo(A+B+1,7)+1));
//        theta = theta + (2*%pi/3);
//        lx = xy(:,1).*cos(theta)-xy(:,2).*sin(theta);
//        ly = xy(:,1).*sin(theta)+xy(:,2).*cos(theta)
//        plot(lx,ly,c(modulo(A+B+1,7)+1));
//    //end
////end
//                         

