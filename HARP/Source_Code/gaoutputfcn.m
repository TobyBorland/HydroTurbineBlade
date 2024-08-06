function [state, options,optchanged] = gaoutputfcn(options,state,flag)
%GAOUTPUTFCNTEMPLATE Template to write custom OutputFcn for GA.
%   [STATE, OPTIONS, OPTCHANGED] = GAOUTPUTFCNTEMPLATE(OPTIONS,STATE,FLAG)
%   where OPTIONS is an options structure used by GA. 
%
%   STATE: A structure containing the following information about the state 
%   of the optimization:
%             Population: Population in the current generation
%                  Score: Scores of the current population
%             Generation: Current generation number
%              StartTime: Time when GA started 
%               StopFlag: String containing the reason for stopping
%              Selection: Indices of individuals selected for elite,
%                         crossover and mutation
%            Expectation: Expectation for selection of individuals
%                   Best: Vector containing the best score in each generation
%        LastImprovement: Generation at which the last improvement in
%                         fitness value occurred
%    LastImprovementTime: Time at which last improvement occurred
%
%   FLAG: Current state in which OutputFcn is called. Possible values are:
%         init: initialization state 
%         iter: iteration state
%    interrupt: intermediate state
%         done: final state
% 		
%   STATE: Structure containing information about the state of the
%          optimization.
%
%   OPTCHANGED: Boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH, GA, GAOPTIMSET

%   Copyright 2004-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2007/08/03 21:23:22 $

global SpdCtrl RootDir filename_main NumVars
        
optchanged = false;

switch flag
    case 'init'
        
        fprintf('\n\n');
        disp('Initial population and constraints defined. Beginning iteration...'); %This gets displayed when the genetic algorithm is initiated
        
        %After the ga_CustomCreation function defines the initial population, 
        %this function writes that initial population to a text file, so the
        %Plot_Initial_Pop function can read the data and create a plot of
        %the initial population.  Also, the header for the _Detailed GA
        %Output output file is written in this function.
        
        %need to make a string that tells fprintf to print the right length format
        s = [];
        for n = 1:NumVars
            a = [37    51    46    52   102    92   116]; %ASCII representation of: %3.4f\t
            s = [s a];
        end
        s(end) = 110; %change last character to "n", (next line statment = \n)
        format1 = char(s);
        
        s = [];
        for n = 1:(NumVars - 10 - (SpdCtrl==0))
            a = abs(['Thick' num2str(n) '	']);
            s = [s a];
        end
        format2 = char(s);
            
% % % % %         fid1 = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '_Detailed_GA_Output.dat'],'Wt');
        fid2 = fopen([RootDir '\Output_Files\' filename_main '\' filename_main '_Initial_Population.dat'],'wt');
        
% % % % %         if SpdCtrl == 1
% % % % %         fprintf(fid1,['F:       	Twist1	Twist2	Twist3	Twist4	Twist5	Chord1	Chord2	Chord3	Chord4	Chord5	'	format2 '\n']);
% % % % %         else
% % % % %         fprintf(fid1,['F:       	Twist1	Twist2	Twist3	Twist4	Twist5	Chord1	Chord2	Chord3	Chord4	Chord5	'	format2 'RotorSpd\n']);    
% % % % %         end
% % % % %         fprintf(fid1,'============================================================================================================================================================\n');
        
        if SpdCtrl == 0; %Rotor Speed becomes the last variable
        fprintf(fid2,['Twist1	Twist2	Twist3	Twist4	Twist5	Chord1	Chord2	Chord3	Chord4	Chord5	' format2 'RPM\n']);
        fprintf(fid2,format1,state.Population');     
        else             %Only 15 varaibles
        fprintf(fid2,['Twist1	Twist2	Twist3	Twist4	Twist5	Chord1	Chord2	Chord3	Chord4	Chord5	' format2 '\n']);
        fprintf(fid2,format1,state.Population');   
        end
% % % % %         fclose(fid1);
        fclose(fid2);
        
        %figure('Name','Initial Population','NumberTitle','off'); %Opens the figure for the initial population plot
        %Plot_Initial_Pop; %calls the function to plot the initial population
        
%         maximum=max(state.Population,[],1)';
%         minimum=min(state.Population,[],1)';
%         sigma=std(state.Population,[],1)';
        
    case {'iter','interrupt'}
%This code is not currently being used, but can be used to create
%customized plots while the genetic algorithm is running, instead of using
%the default plots that MATLAB displays.


%         clc;
%         disp('Iterating ...')
%         disp('Generation:');
%          disp(state.Generation);
%          disp('BestVal:');
%          disp(min(state.Score))
%          bestindex=find(state.Score==min(state.Score),1);
%          bestofgen=state.Population(bestindex,:);
%          disp(bestofgen');
%          
                  
% %         subplot(2,2,1), plot(state.Best); title('Best Fitness Value');
% %         subplot(2,2,2), plot(bestofgen(1:5)); title('Twist Distribution');
% %         subplot(2,2,3), plot(bestofgen(6:10)); title('Chord Distribution');
% %         subplot(2,2,4), plot(bestofgen(11:15)); title('%Thickness Distribution');
%         %output generation data to output files
%         if(state.Generation==1)
%             fidmax = fopen('Optimization Output\Maximum.dat', 'wt');
%             fidmin = fopen('Optimization Output\Minimum.dat', 'wt');
%             fidstd = fopen('Optimization Output\Deviation.dat', 'wt');
%             fidbst = fopen('Optimization Output\BestofGen.dat','wt');
%             
%             fprintf(fidmax,'%12s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s\n','Generation','Twist1','Twist2','Twist3','Twist4','Twist5','Chord1','Chord2','Chord3','Chord4','Chord5','Thick1','Thick2','Thick3','Thick4','Thick5');
%             fprintf(fidmin,'%12s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s\n','Generation','Twist1','Twist2','Twist3','Twist4','Twist5','Chord1','Chord2','Chord3','Chord4','Chord5','Thick1','Thick2','Thick3','Thick4','Thick5');
%             fprintf(fidstd,'%12s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s\n','Generation','Twist1','Twist2','Twist3','Twist4','Twist5','Chord1','Chord2','Chord3','Chord4','Chord5','Thick1','Thick2','Thick3','Thick4','Thick5');
%             fprintf(fidbst,'%12s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s\n','Generation','Twist1','Twist2','Twist3','Twist4','Twist5','Chord1','Chord2','Chord3','Chord4','Chord5','Thick1','Thick2','Thick3','Thick4','Thick5');
%         else
%             fidmax = fopen('Optimization Output\Maximum.dat', 'at');
%             fidmin = fopen('Optimization Output\Minimum.dat', 'at');
%             fidstd = fopen('Optimization Output\Deviation.dat', 'at');
%             fidbst = fopen('Optimization Output\BestofGen.dat', 'at');
%         end
%         fprintf(fidmax,'%12.0f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',[state.Generation max(state.Population,[],1)]);
%         fprintf(fidmin,'%12.0f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',[state.Generation min(state.Population,[],1)]);
%         fprintf(fidstd,'%12.0f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',[state.Generation std(state.Population,[],1)]);
%         fprintf(fidbst,'%12.0f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',[state.Generation bestofgen]);
%         fclose('all');
%         %maximum=max(state.Population,[],1)'
%         %minimum=min(state.Population,[],1)'
%         %sigma=std(state.Population,[],1)'

    case 'done'
        disp('Generating output files (please wait)...'); %Displays after the genetic algorithm completes...but other tasks will be performed after this
        fclose('all');
end


