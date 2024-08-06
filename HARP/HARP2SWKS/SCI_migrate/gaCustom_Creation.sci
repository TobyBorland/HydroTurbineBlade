function Population = gaCustom_Creation(GenomeLength,unused1,options)
// This function creates the initial population for the genetic algorithm.
// Basically, this file needs to create a random distribution of the genetic
// algorithm variables (control points for the Bezier curve fits), but these 
// variables MUST obey the upper and lower bounds, and linear inequality 
// constraints defined in HARP_Opt.m
// 
// Inputs:    GenomeLength: the number of variables that the genetic algorithm is using. 16 if fixed speed and 15 if variable speed turbine.
//            unused1: apparently an unused argument defined by MATLAB
//            options: this is the options structure created for the genetic algorithm, this options struction is defined in HARP_Opt.m
// 
// Outputs:   Population: a matrix, containing all of the 15 or 16 dimensional vectors of the initial population

global SpdCtrl RotorRad ThickMethod Thickness_values NumSeg...
       TwistLB TwistUB ChordLB ChordUB ThickUB ThickLB OmgMin OmgMax;

// =========================================================================// 
// Temporary variables for debugging
//  LB = [0 0 0 0 0 0 0 0 0 0 0.2 0.4 0.8 20];
//  UB = [10 10 10 10 10 10 10 10 10 10 0.4 0.6 0.9 40];
// =========================================================================// 
totalPopulation = sum(options.PopulationSize);
// initPopProvided = size(options.InitialPopulation,1);
initPopProvided = 0;
// individualsToCreate = totalPopulation - initPopProvided
individualsToCreate = totalPopulation;

//  Initialize Population to be created
Population = zeros(totalPopulation,GenomeLength);
//  Use initial population provided already
if initPopProvided > 0
    Population(1:initPopProvided,:) = options.InitialPopulation;
end

// Preallocate
Twist = zeros(1,5);
Chord = zeros(1,5);
if ThickMethod == 1;
Thick = zeros(1,length(Thickness_values)-1);    
elseif ThickMethod == 2;
Thick = zeros(1,length(Thickness_values));    
end
X = zeros(1,GenomeLength);


fprintf(1,'\nCreating initial population of feasible individuals ... if the\n');
  fprintf(1,'number of feasible individuals found does not increase after a\n');
  fprintf(1,'short period of time, you may have set unfeasible input values\n');
  fprintf(1,'and will need to cancel (CTRL-C) HARP_Opt and input new values:\n\n');

//  //  //  //  //  fprintf(1,'Feasible Individuals   Total Individuals\n');
//  //  //  //  //  fprintf(1,'  Found (need %g)         Created\n',individualsToCreate);
//  //  //  //  //  fprintf(1,'------------------------------------------\n');
//  //  //  //  //  fprintf(1,'                               ');


H = waitbar(0,'','Name','Creating initial population...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(H,'canceling',0)

C=1;
D=1;
while D <= individualsToCreate      
    
    // Twist and Chord
    Twist(1,1) = rand*(TwistUB(1)-TwistLB(1)) + TwistLB(1);
    Chord(1,1) = rand*(ChordUB(1)-ChordLB(1)) + ChordLB(1);
    for k = 2:5
        if Twist(1,k-1) > TwistUB(k)
           ubT = TwistUB(k);
        else
           ubT = Twist(1,k-1); 
        end
        if Chord(1,k-1) > ChordUB(k)
           ubC = ChordUB(k); 
        else
           ubC = Chord(1,k-1); 
        end
    lbT = TwistLB(k);
    lbC = ChordLB(k);
    Twist(1,k) = rand*(ubT-lbT) + lbT;
    Chord(1,k) = rand*(ubC-lbC) + lbC;
    end
    
    // Percent Thickness radial stations
    // if bounds are entered for the Percent Thickness
    if length(Thickness_values) > 1

        if isempty(ThickLB) & isempty(ThickUB)
           ThickLB(1:length(Thickness_values)) = NaN;
           ThickUB = ThickLB;
        end

        if ThickMethod == 1
           NumAFvars = length(Thickness_values)-1; 
        elseif ThickMethod == 2
           NumAFvars = length(Thickness_values);  
        end

        Thick(1,1:NumAFvars) = NaN; // preallocating
        // first take care of the ones with BOTH lower AND upper bounds
        for n = 1:NumAFvars
            if ~isnan(ThickLB(n)) & ~isnan(ThickUB(n))
                Thick(n) = rand*(ThickUB(n) - ThickLB(n)) + ThickLB(n); 
            end
        end

        for n = NumAFvars:-1:1

            //  only upper bound exists
            if isnan(ThickLB(n)) & isfinite(ThickUB(n))

                randUB(n) = ThickUB(n);

                // find a lower bound
                for  N = n:-1:1
                   if isfinite(ThickLB(N)) & isnan(Thick(N))
                       randLB(n) = ThickLB(N);
                       break
                   elseif isfinite(Thick(N))
                       randLB(n) = Thick(N);
                       break
                   elseif isnan(ThickLB(N)) & isfinite(ThickUB(N)) & N ~= n
                       randLB(n) = ThickUB(N);
                       break              
                   elseif N == 1
                       randLB(n) = 0;
                   else
                       continue
                   end
                end
            end

            //  only lower bound exists
            if isfinite(ThickLB(n)) & isnan(ThickUB(n)) 

                randLB(n) = ThickLB(n);

                // find an upper bound
                for  N = n:NumAFvars
                   if isfinite(ThickUB(N)) & isnan(Thick(N))
                       randUB(n) = ThickUB(N);
                       break
                   elseif isfinite(Thick(N))
                       randUB(n) = Thick(N);
                       break
                   elseif N == NumAFvars
                       randUB(n) = 1;
                   else
                       continue
                   end
                end               
            end

            //  no bounds exist
            if isnan(ThickLB(n)) & isnan(ThickUB(n)) 
                // find a lower bound
                for  N = n:-1:1
                   if isfinite(ThickLB(N)) & isnan(Thick(N))
                       randLB(n) = ThickLB(N);
                       break
                   elseif isfinite(Thick(N))
                       randLB(n) = Thick(N);
                       break
                   elseif isnan(ThickLB(N)) & isfinite(ThickUB(N)) & N ~= n
                       randLB(n) = ThickUB(N);
                       break
                   elseif N == 1
                       randLB(n) = 0;
                   else
                       continue
                   end
                end

                // find an upper bound
                for  N = n:NumAFvars
                   if isfinite(ThickUB(N)) & isnan(Thick(N))
                       randUB(n) = ThickUB(N);
                       break
                   elseif isfinite(Thick(N))
                       randUB(n) = Thick(N);
                       break
                   elseif N == NumAFvars
                       randUB(n) = 1;
                   else
                       continue
                   end
                end  
            end

            //  both bounds exist
            if isfinite(ThickLB(n)) & isfinite(ThickUB(n))
            randLB(n) = ThickLB(n);
            randUB(n) = ThickUB(n);
            end

            if isnan(Thick(n))
            Thick(n) = rand*(randUB(n) - randLB(n)) + randLB(n);
            end
        end

       X(1,1:10+NumAFvars) = [Twist(1,:) Chord(1,:) Thick(1,:)];  
    else
       X(1,1:10) = [Twist(1,:) Chord(1,:)];
    end
 
    if SpdCtrl == 0;   
    X(1,GenomeLength) = rand*(OmgMax-OmgMin) + OmgMin;
    end
    
    // evaluate the randomly created individual to see if it produces any kind of error
    F = Main(X);
    if F(1) ~= Inf
        Population(D,:) = X;
        D = D+1;
    end
    
//  //  //  //  //      fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
//  //  //  //  //      fprintf(1,'       %-10.0f\t\t\t%-10.0f',D-1,C);
    
    //  Check for Cancel button press
    if getappdata(H,'canceling')
       delete(H) //  DELETE the waitbar; don't try to CLOSE it.
       fclose('all');
       error('Program Terminated: User aborted program during creation of initial population.');
    end
    //  Report current estimate in the waitbar's message field
    waitbar((D-1)/individualsToCreate,H,sprintf('Feasible Found = %g (Need %g), Total Created = %g',D-1,individualsToCreate,C))

    C = C+1; // update counter for next iteration
end
delete(H) //  DELETE the waitbar; don't try to CLOSE it.
// fprintf(1,' Done. \n');

// =========================================================================// 
// Plot the initial population
r = zeros(NumSeg,size(Population,1));
tw = zeros(NumSeg,size(Population,1));
c = zeros(NumSeg,size(Population,1));
pt = zeros(NumSeg,size(Population,1));
dt = zeros(NumSeg,size(Population,1));
for k = 1:size(Population,1)
    [ShapeError RElm TWIST CHORD PERCENT_THICKNESS DIMENSIONAL_THICKNESS] = Define_Blade_Shape(Population(k,:));
    r(:,k) = RElm;
    tw(:,k) = TWIST;
    c(:,k) = CHORD;
    pt(:,k) = PERCENT_THICKNESS;
    dt(:,k) = DIMENSIONAL_THICKNESS;
end
Fig1 = figure(1);
set(Fig1,'color','white','name','Initial Population','numberTitle','off');

if SpdCtrl == 1;
subplot(4,1,1); plot(r,tw,'-b','LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Pre-twist (deg)'); xlim([0 RotorRad]); box on;
subplot(4,1,2); plot(r,c,'-r','LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Chord (m)'); xlim([0 RotorRad]); box on;
subplot(4,1,3); plot(r,pt,'-','Color',[0 0.5 0],'LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Thickness (%)'); xlim([0 RotorRad]); box on;
subplot(4,1,4); plot(r,dt,'-k','LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Thickness (m)'); xlim([0 RotorRad]); box on;
elseif SpdCtrl == 0; // Fixed Speed
subplot(4,3,[1 2]); plot(r,tw,'-b','LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Pre-twist (deg)'); xlim([0 RotorRad]); box on;
subplot(4,3,[4 5]); plot(r,c,'-r','LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Chord (m)'); xlim([0 RotorRad]); box on;
subplot(4,3,[7 8]); plot(r,pt,'-','Color',[0 0.5 0],'LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Thickness (%)'); xlim([0 RotorRad]); box on;
subplot(4,3,[10 11]); plot(r,dt,'-k','LineWidth',1,'LineSmoothing','off');xlabel('Blade Radius (m)');ylabel('Thickness (m)'); xlim([0 RotorRad]); box on;
subplot(4,3,[3 12]); plot(1:1:size(Population,1),Population(:,GenomeLength),'ok','LineSmoothing','off'); xlabel('Individual');ylabel('Rotor Speed (rpm)'); xlim([0 size(Population,1)]); box on;
end
endfunction
