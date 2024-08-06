function [ShapeError RElm TWIST CHORD PERCENT_THICKNESS DIMENSIONAL_THICKNESS R_CHORD_CP CHORD_CP R_TWIST_CP TWIST_CP THICK_CP] = Define_Blade_Shape(x)

global NumSeg RotorRad HubRad Thickness_values ThickMethod...
       CircleRoot minRootChord maxRootChord RootTranSt RootTranEnd RadialSpacing

%=========================================================================%

%This section defines the blade geometery
ShapeError = 0; %initial value

%Define the Twist and Chord control points
Twist_CP(:,1) = x(1:5)';
Chord_CP(:,1) = x(6:10)';

%Define the radial locations of the control points: using equal or cosine spacing
usecosine = 1;
if usecosine == 0; %Using equal distance spacing for the design variables
   if CircleRoot == 0
   rad_CP = linspace(HubRad,RotorRad,5)';    
   elseif CircleRoot == 1
   rad_CP = linspace(RootTranEnd*RotorRad,RotorRad,5)';       
   end
elseif usecosine == 1; %Using half cosine spacing for the design variables
   if CircleRoot == 0
   rad_CP = hcosspace(HubRad,RotorRad,5,0);    
   elseif CircleRoot == 1
   rad_CP = hcosspace(RootTranEnd*RotorRad,RotorRad,5,0);      
   end
end

%Define the radia locations of the WT_Perf blade elements: radial or cosine spacing
if RadialSpacing == 0; %equal spacing of blade elements  
    %Define the radius of the center of the blade elements
    delElm = (RotorRad - HubRad)/NumSeg;
    RElm = ((HubRad + delElm/2):delElm:((HubRad + delElm/2)+(NumSeg-1)*delElm))'; %Values at the center of each element (as defined in WT_Perf)

elseif RadialSpacing == 1; %cosine spacing of blade elements
    RElm_cs_EndPts = hcosspace(HubRad,RotorRad,NumSeg+1,3); %cosine spacing values at the boundaries of each element
    %calculate the centers of the elements for cosine spacing
    RElm = zeros(NumSeg,1);
    for n = 1:NumSeg
    RElm(n,1) = (RElm_cs_EndPts(n+1) + RElm_cs_EndPts(n))/2;    
    end
end

if CircleRoot == 0; %Non-circular root
    %From the control points, interpolate using a Bezier curve
    B_Chord = Bezier([rad_CP Chord_CP],NumSeg);
    B_Twist = Bezier([rad_CP Twist_CP],NumSeg);
    %Now interpolate to spacing used in RElm
    CHORD = interp1(B_Chord(:,1),B_Chord(:,2),RElm);
    TWIST = interp1(B_Twist(:,1),B_Twist(:,2),RElm);
    R_CHORD_CP = rad_CP;
    R_TWIST_CP = rad_CP;
    CHORD_CP = Chord_CP;
    TWIST_CP = Twist_CP;
    
    %Now define the percent thickness distribution
    if ThickMethod == 1 && length(Thickness_values) > 1;        %Piecewise Constant
        THICK_CP(:,1) = RotorRad*x(11:(10+length(Thickness_values)-1));
        THICK_CP = THICK_CP;
        %create a step function
        Thickness(1:length(RElm),1) = Thickness_values(end);
        a = 1; %intial index value
        for n = 1:length(THICK_CP)
            for k = a:length(RElm)
                if RElm(k) < THICK_CP(n)
                Thickness(k,1) = Thickness_values(n);
                elseif RElm(k) > THICK_CP(n)
                a=k;
                break
                end
            end
        end
        Thickness(end+1:length(RElm)) = Thickness_values(end);
        PERCENT_THICKNESS = Thickness;
        DIMENSIONAL_THICKNESS = CHORD.*PERCENT_THICKNESS./100;
    elseif ThickMethod == 2 && length(Thickness_values) > 1;    %Piecewise Linear
        THICK_CP(:,1) = RotorRad*x(11:(10+length(Thickness_values)));
        THICK_CP = THICK_CP;
        %Make sure that none of the control points are equal, so the interpolation does not crash
        for n = 2:length(THICK_CP);
            if THICK_CP(n) <= THICK_CP(n-1)
                THICK_CP(n) = THICK_CP(n-1) + 0.0001;
            end
        end
        Thickness = interp1(THICK_CP,Thickness_values,RElm);
        Thickness(RElm<THICK_CP(1)) = Thickness_values(1);
        Thickness(RElm>THICK_CP(end)) = Thickness_values(end);
        PERCENT_THICKNESS = Thickness;
        DIMENSIONAL_THICKNESS = CHORD.*PERCENT_THICKNESS./100;
    else
        THICK_CP = Thickness_values(1);
        PERCENT_THICKNESS(1:NumSeg,1) = Thickness_values(1);
        DIMENSIONAL_THICKNESS = CHORD.*PERCENT_THICKNESS./100;
    end
    
elseif CircleRoot == 1; %Circular root
    ShapeError = 1; %set error = 1 initially
    
    %Define new control points for circular root chord
    CP_cr1 = RootTranSt*RotorRad;
    CP_cr2 = RootTranEnd*RotorRad;
    
    if isempty(minRootChord); minChord = min(x(6:10)); else minChord = minRootChord; end
    if isempty(maxRootChord); maxChord = max(x(6:10)); else maxChord = maxRootChord; end
    if minChord > maxChord; minChord = maxChord; end;
    delChord = 0.01; %m

    RtChord = (minChord(1):delChord:maxChord(1))';
    
    for c = 1:length(RtChord)
        RootChord = RtChord(c);
        
        BlendMethod = 2; %I've experimented with different ways of blending the root region, currently this is hardcoded in and not an option for the user to change
        if BlendMethod == 1;
        rad_CP_cr = [CP_cr1;(CP_cr1+CP_cr2)/2;rad_CP];
        Chord_CP_cr = [RootChord;(RootChord+Chord_CP(1))/2;Chord_CP];
        R_CHORD_CP = rad_CP_cr;
        CHORD_CP = Chord_CP_cr;
        elseif BlendMethod == 2;
        b = 0.7;
        d1 = CP_cr2 - b*(CP_cr2 - CP_cr1);
        d2 = CP_cr1 + b*(CP_cr2 - CP_cr1);
        Chord_rad_CP_cr = [CP_cr1;d1;d2;rad_CP];
        Chord_CP_cr = [RootChord;RootChord;Chord_CP(1);Chord_CP];
        R_CHORD_CP = Chord_rad_CP_cr;
        CHORD_CP = Chord_CP_cr;
        end

        B_Chord = Bezier([Chord_rad_CP_cr Chord_CP_cr],NumSeg); %returns equal spaced vector
        %Now interpolate to spacing used in RElm
        Chord = interp1(B_Chord(:,1),B_Chord(:,2),RElm);
        Chord(isnan(Chord)) = RootChord; %only changes values near the root
        CHORD = Chord;
                
        %Now find the location of max chord
        r_ChordMax = RElm(Chord==max(Chord));
        r_ChordMax = r_ChordMax(end); %incase multiple maximums exist, choose the last one
%         hold on; subplot(4,1,1); plot(R_CHORD_CP,CHORD_CP,'sk',RElm,CHORD,'x-r'); legend('Design Variables','Chord');
        if r_ChordMax > rad_CP(1) && c ~= length(RtChord); 
            continue %skip to next iteration with larger chord value to see if this makes a feasible blade shape
        end       
        Twist_rad_CP = [r_ChordMax;rad_CP(2:end)];
        R_TWIST_CP = Twist_rad_CP;
        TWIST_CP = Twist_CP;
        B_Twist = Bezier([R_TWIST_CP TWIST_CP],NumSeg); %returns equal spaced vector
        %Now interpolate to spacing used in RElm
        Twist = interp1(B_Twist(:,1),B_Twist(:,2),RElm);
        Twist(isnan(Twist)) = max(Twist);
        TWIST = Twist;
%         hold on; subplot(4,1,2); plot(R_TWIST_CP,TWIST_CP,'sk',RElm,TWIST,'x-b'); legend('Design Variables','Pre-Twist');
        

        %Now define the percent thickness distribution
     if length(Thickness_values) > 1
        if ThickMethod == 1;        %Piecewise Constant
            THICK_CP(:,1) = RotorRad*x(11:(10+length(Thickness_values)-1));
            %create a step function
            Thickness(1:length(RElm),1) = Thickness_values(end);
            a = 1; %intial index value
            for n = 1:length(THICK_CP)
                for k = a:length(RElm)
                    if RElm(k) < THICK_CP(n)
                    Thickness(k,1) = Thickness_values(n);
                    elseif RElm(k) > THICK_CP(n)
                    a=k;
                    break
                    end
                end
            end
            %Thickness(end+1:length(RElm)) = Thickness_values(end);
            Thickness(RElm<=RotorRad*RootTranSt) = 100;
            dimT = CHORD.*Thickness./100; %Dimensional Thickness
            
        elseif ThickMethod == 2;    %Piecewise Linear
            THICK_CP(:,1) = RotorRad*x(11:(10+length(Thickness_values)));
            %Make sure that none of the control points are equal, so the interpolation does not crash
            for n = 2:length(THICK_CP);
                if THICK_CP(n) <= THICK_CP(n-1)
                    THICK_CP(n) = THICK_CP(n-1) + 0.0001;
                end
            end
            Thickness = interp1(THICK_CP,Thickness_values,RElm);
            Thickness(RElm<THICK_CP(1)) = Thickness_values(1);
            Thickness(RElm>THICK_CP(end)) = Thickness_values(end);
            Thickness(RElm<=RotorRad*RootTranSt) = 100;
            dimT = CHORD.*Thickness./100; %Dimensional Thickness
        end
     else %only one airfoil was used
     THICK_CP = r_ChordMax;
     Thickness(RElm<=RotorRad*RootTranSt,1) = 100;
     Thickness(RElm>RotorRad*RootTranSt,1) = Thickness_values(1);
     dimT = CHORD.*Thickness./100; %Dimensional Thickness
     end
        
            %For the Piecewise Linear case, if the location of Max Chord is inboard of the 1st design
            %variable for thickness, blend to the location of the 1st design variable.  However, for 
            %the Piecewise Constant case, we will blend the location of max chord.
            if ThickMethod == 2; F = false; else F = true; end;
            if (r_ChordMax > THICK_CP(1)) + F
               t_i = RElm<=r_ChordMax & RElm>=RotorRad*RootTranSt;
               x2 = r_ChordMax;
               fx2 = dimT(RElm==r_ChordMax);
               DFDX = FiniteDiff(RElm,dimT);
               dfx2 = DFDX(RElm==r_ChordMax);
            else
               t_i = RElm<=THICK_CP(1) & RElm>=RotorRad*RootTranSt;
               x2 = RElm(FindInd(RElm,THICK_CP(1)));
               fx2 = dimT(FindInd(RElm,THICK_CP(1)));
               DFDX = FiniteDiff(RElm,dimT);
               dfx2 = DFDX(FindInd(RElm,THICK_CP(1)));
            end
            
            t_vals = RElm(t_i);
            x1 = RotorRad.*RootTranSt;
            fx1 = RootChord; 
            dfx1 = 0;
            dimT(t_i) = CubeFit2(x1,x2,fx1,fx2,dfx1,dfx2,t_vals);
            Thickness_Corrected = 100.*dimT./Chord;
            PERCENT_THICKNESS = Thickness_Corrected;
            DIMENSIONAL_THICKNESS = CHORD.*PERCENT_THICKNESS./100;
            
%             hold on; subplot(4,1,3); plot(THICK_CP,Thickness_values,'sk',RElm,PERCENT_THICKNESS,'x-g'); legend('Design Variables','% Thickness (t/c)');
%             hold on; subplot(4,1,4); plot(RElm,DIMENSIONAL_THICKNESS,'x-k'); legend('Dimensional Thickness (t)');
                
        
        %Now calculate the slope of the non-dimensional AND dimensional thickness, and make sure
        %that BOTH are monotonically decreasing
        ddimTdr = FiniteDiff(RElm,DIMENSIONAL_THICKNESS);
        dTdr = FiniteDiff(RElm,PERCENT_THICKNESS);
        tol = 0.0001;
        if any(ddimTdr>tol) == false && any(dTdr>tol) == false;
           ShapeError = 0;
           %break %end the loop, the dimensional & percent thickness are now both monotonically decreasing, hurrah
           return
        end
    
    end
    
%     if ShapeError == 1;     %c == length(RtChord)
%         %We tried all the possible root cord values, but could not make a valid blade shape
%             %ShapeError = 1;
%             RElm = [];
%             TWIST = [];
%             CHORD = [];
%             PERCENT_THICKNESS = [];
%             DIMENSIONAL_THICKNESS = [];
%             %Print a warning:
%             %fprintf('Warning: Could not satisfy constraints for monotonically decreasing\nthickness distributions. Please check for reasonable inputs under\nthe "Circular Root" section, continuing optimization...\n');
%             return
%     end
        
end

% if plot == 1;
% =========================================================================%
% Plots
% figure('name','Initial Population');
% subplot(4,1,1); plot(R_CHORD_CP,CHORD_CP,'sk',RElm,CHORD,'x-r'); legend('Design Variables','Chord');
% subplot(4,1,2); plot(R_TWIST_CP,TWIST_CP,'sk',RElm,TWIST,'x-b'); legend('Design Variables','Pre-Twist');
% if ThickMethod == 2
% subplot(4,1,3); plot(THICK_CP,Thickness_values,'sk',RElm,PERCENT_THICKNESS,'x-g'); legend('Design Variables','% Thickness (t/c)');
% else
%     for k = 1:length(Thickness_values)-1;
%         tcp(k,1) = 0.5*(Thickness_values(k)+Thickness_values(k+1));
%     end
% subplot(4,1,3); plot(THICK_CP,tcp,'sk',RElm,PERCENT_THICKNESS,'x-g'); legend('Design Variables','% Thickness (t/c)');
% end
% subplot(4,1,4); plot(RElm,DIMENSIONAL_THICKNESS,'x-k'); legend('Dimensional Thickness (t)');
% =========================================================================%
% end