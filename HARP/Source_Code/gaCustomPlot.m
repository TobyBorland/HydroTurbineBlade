function state = gaCustomPlot(options,state,flag)
%options, state, and flag are delivered from the MATLAB genetic algorithm functions

global RotorRad ThickMethod Thickness_values CircleRoot SpdCtrl StructuralOpt OptMethod
persistent best m distance GAFig

%used for plot avg. distance (same as MATLAB gaplotdistance.m)
samples = 20;
choices = ceil(sum(options.PopulationSize) * rand(samples,2));

if OptMethod == 0
    PlotLabel = 'Performance Fitness Value';
    Modifier = 1;
else
    PlotLabel = '(penalized) Annual Energy Production (kWh/yr)';
    Modifier = -1;
end


switch flag
    case 'init'
        GAFig = gcf;
        figure(GAFig);
        set(GAFig,'color','white')
        
        [BestValue,BestIndex] = min(state.Score(:,1));
        BestIndiv = state.Population(BestIndex,:);
        
        %===================== Plot Best Blade Geometry ==================%
        [ShapeError RElm TWIST CHORD PERCENT_THICKNESS DIMENSIONAL_THICKNESS...
         R_CHORD_CP CHORD_CP R_TWIST_CP TWIST_CP THICK_CP] = Define_Blade_Shape(BestIndiv);
   if StructuralOpt == 0
        subplot(4,3,[1 2]);
        hold on;
        plot(R_TWIST_CP,TWIST_CP,'sk',RElm,TWIST,'+-b');
        legend('Control Points','Pre-Twist');
        plot(RElm,TWIST,'-b','LineWidth',1.5,'LineSmoothing','off');
        ylabel('Pre-twist (deg)');
        xlabel('Blade Radius (m)');
        xlim([0 RotorRad]);
        box on
        hold off;
          if SpdCtrl == 0
            title(['Best Individual: Blade Geometry, Rotor Speed=' num2str(BestIndiv(end),'%-3.2f') ' rpm']);
          else
            title('Best Individual: Blade Geometry');
          end
        if CircleRoot == 1
          subplot(4,3,[4 5]);
          hold on;
          plot(R_CHORD_CP(1:3),CHORD_CP(1:3),'ok',R_CHORD_CP(4:end),CHORD_CP(4:end),'sk',RElm,CHORD,'+-r');
          legend('Circular Root Control Points','Control Points','Chord');
          plot(RElm,CHORD,'-r','LineWidth',1.5,'LineSmoothing','off');        
          ylabel('Chord (m)'); 
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]);
          box on
          hold off;
        else
          subplot(4,3,[4 5]);
          hold on;
          plot(R_CHORD_CP,CHORD_CP,'sk',RElm,CHORD,'+-r');
          legend('Control Points','Chord'); 
          plot(RElm,CHORD,'-r','LineWidth',1.5,'LineSmoothing','off');        
          ylabel('Chord (m)'); 
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]);
          box on
          hold off;
        end
        if ThickMethod == 2
          subplot(4,3,[7 8]);
          hold on;
          plot(THICK_CP,Thickness_values,'sk');
          plot(RElm,PERCENT_THICKNESS,'+-','Color',[0 0.5 0]);
          legend('Control Points','% Thickness (t/c)');
          plot(RElm,PERCENT_THICKNESS,'-','Color',[0 0.5 0],'LineWidth',1.5,'LineSmoothing','off');          
          ylabel('Thickness (%)');
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]);
          box on
          hold off;
        else
          subplot(4,3,[7 8]);         
          hold on;
            if length(Thickness_values) > 1
              tcp = zeros(1,length(Thickness_values));
              for k = 1:length(Thickness_values)-1;
                  tcp(k,1) = 0.5*(Thickness_values(k)+Thickness_values(k+1));
              end
            else
                tcp = Thickness_values(1);
            end
          plot(THICK_CP,tcp,'sk');    
          plot(RElm,PERCENT_THICKNESS,'+-','Color',[0 0.5 0]);
          legend('Control Points','% Thickness (t/c)'); 
          plot(RElm,PERCENT_THICKNESS,'-','Color',[0 0.5 0],'LineWidth',1.5,'LineSmoothing','off');          
          ylabel('Thickness (%)');
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]);
          box on
          hold off;
        end
        subplot(4,3,[10 11]);
        hold on;
        plot(RElm,DIMENSIONAL_THICKNESS,'+-k');
        %legend('Dimensional Thickness (t)');
        plot(RElm,DIMENSIONAL_THICKNESS,'-k','LineWidth',1.5,'LineSmoothing','off');           
        ylabel('Thickness (m)');
        xlabel('Blade Radius (m)');
        xlim([0 RotorRad]);
        box on
        hold off;
        %=================================================================%
        
        %=================== Plot Best and Avg. Fitness ==================%
        subplot(4,3,[3 6]);
        best = min(state.Score(:,1));
%         m    = meanf(state.Score);
%         plot(state.Generation,best,'.-k',state.Generation,m,'.b','LineWidth',1.5,'MarkerSize',15,'LineSmoothing','off');
        hold on;
        plot(state.Generation,Modifier.*best,'.k','MarkerSize',15);
        plot(state.Generation,Modifier.*best,'-k','LineWidth',1.5,'LineSmoothing','off');
        hold off;
        set(gca,'xlim',[0 1]);
        set(gca,'xtick',[0 1]);
        xlabel('Generation');
        ylabel(PlotLabel);
        box on
%         legend('Best fitness','Mean fitness');
        title('Best: ');
        %=================================================================%

        %=================== Plot Avg. Distance===========================%
        subplot(4,3,[9 12]);
        population = state.Population;
        distance = 0;
        for i = 1:samples
            d = population(choices(i,1),:) - population(choices(i,2),:);
            distance = distance + sqrt( sum ( d.* d));
        end
        distance = distance/samples;
        hold on;
        plot(state.Generation,distance,'.k','MarkerSize',15);
        plot(state.Generation,distance,'-k','LineWidth',1.5,'LineSmoothing','off');
        hold off;        
        set(gca,'xlimmode','manual','zlimmode','manual','alimmode','manual');
        set(gca,'xlim',[0 1]);
        set(gca,'xtick',[0 1]);
        xlabel('Generation','interp','none');
        ylabel('Average Distance');
        box on
        title('Average Distance Between Individuals')
        %=================================================================%
    elseif StructuralOpt == 1
        xy = mogaCustomPlot(options,state,flag);
        xy(isnan(xy(:,1)),:)=[];  %get rid of any solutions which failed in some way (i.e. NaN)
        POpts = intersect([state.Score(:,1) state.Score(:,2)],xy,'rows');
        INFpts = setxor([state.Score(:,1) state.Score(:,2)],xy,'rows');
        INFpts(isnan(INFpts(:,1)),:)=[];
        hold on
        plot(Modifier.*POpts(:,1),POpts(:,2),'or','MarkerSize',4,'MarkerFaceColor','r','LineSmoothing','off');
        plot(Modifier.*INFpts(:,1),INFpts(:,2),'ok','MarkerSize',3,'LineSmoothing','off');
        hold off
        xlabel(PlotLabel);
        ylabel('Total Blade Mass (kg)');
        box on
        title(['Pareto Front, Generation = ' num2str(state.Generation)]);
        legend('Pareto Optimal','Inferior',2);
    end
    
    case 'iter'
        figure(GAFig);
        
        [BestValue,BestIndex] = min(state.Score(:,1));
        BestIndiv = state.Population(BestIndex,:);
               
        %===================== Plot Best Blade Geometry ==================%
        [ShapeError RElm TWIST CHORD PERCENT_THICKNESS DIMENSIONAL_THICKNESS...
         R_CHORD_CP CHORD_CP R_TWIST_CP TWIST_CP THICK_CP] = Define_Blade_Shape(BestIndiv);
     if StructuralOpt == 0
        subplot(4,3,[1 2]);
        hold on;
        plot(R_TWIST_CP,TWIST_CP,'sk',RElm,TWIST,'+-b');
        legend('Control Points','Pre-Twist');
        plot(RElm,TWIST,'-b','LineWidth',1.5,'LineSmoothing','off');
        ylabel('Pre-twist (deg)');
        xlabel('Blade Radius (m)');
        xlim([0 RotorRad]);
        box on
        hold off;
          if SpdCtrl == 0
            title(['Best Individual: Blade Geometry, Rotor Speed=' num2str(BestIndiv(end),'%-3.2f') ' rpm']);
          else
            title('Best Individual: Blade Geometry');
          end
        if CircleRoot == 1
          subplot(4,3,[4 5]);
          hold on;
          plot(R_CHORD_CP(1:3),CHORD_CP(1:3),'ok',R_CHORD_CP(4:end),CHORD_CP(4:end),'sk',RElm,CHORD,'+-r');
          legend('Circular Root Control Points','Control Points','Chord');
          plot(RElm,CHORD,'-r','LineWidth',1.5,'LineSmoothing','off');        
          ylabel('Chord (m)'); 
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]);
          box on
          hold off;
        else
          subplot(4,3,[4 5]);
          hold on;
          plot(R_CHORD_CP,CHORD_CP,'sk',RElm,CHORD,'+-r');
          legend('Control Points','Chord'); 
          plot(RElm,CHORD,'-r','LineWidth',1.5,'LineSmoothing','off');        
          ylabel('Chord (m)'); 
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]);
          box on
          hold off;
        end
        if ThickMethod == 2
          subplot(4,3,[7 8]);
          hold on;
          plot(THICK_CP,Thickness_values,'sk');
          plot(RElm,PERCENT_THICKNESS,'+-','Color',[0 0.5 0]);
          legend('Control Points','% Thickness (t/c)');
          plot(RElm,PERCENT_THICKNESS,'-','Color',[0 0.5 0],'LineWidth',1.5,'LineSmoothing','off');          
          ylabel('Thickness (%)');
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]); 
          box on
          hold off;
        else
          subplot(4,3,[7 8]);         
          hold on;
            if length(Thickness_values) > 1
              tcp = zeros(1,length(Thickness_values));
              for k = 1:length(Thickness_values)-1;
                  tcp(k,1) = 0.5*(Thickness_values(k)+Thickness_values(k+1));
              end
            else
                tcp = Thickness_values(1);
            end
          plot(THICK_CP,tcp,'sk');    
          plot(RElm,PERCENT_THICKNESS,'+-','Color',[0 0.5 0]);
          legend('Control Points','% Thickness (t/c)'); 
          plot(RElm,PERCENT_THICKNESS,'-','Color',[0 0.5 0],'LineWidth',1.5,'LineSmoothing','off');          
          ylabel('Thickness (%)');
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]); 
          box on
          hold off;
        end
        subplot(4,3,[10 11]);
        hold on;
        plot(RElm,DIMENSIONAL_THICKNESS,'+-k');
        %legend('Dimensional Thickness (t)');
        plot(RElm,DIMENSIONAL_THICKNESS,'-k','LineWidth',1.5,'LineSmoothing','off');           
        ylabel('Thickness (m)');
        xlabel('Blade Radius (m)');
        xlim([0 RotorRad]);
        box on
        hold off;
        %=================================================================%
    
        %=================== Plot Best and Avg. Fitness ==================%
        subplot(4,3,[3 6]);
        best = [best min(state.Score(:,1))];
%         m    = [m meanf(state.Score)];
        X = 0:1:state.Generation;
%         plot(X,best,'.-k',X,m,'.b','LineWidth',1.5,'MarkerSize',15,'LineSmoothing','off');
        hold on;
        plot(X,Modifier.*best,'.k','MarkerSize',15);
        plot(X,Modifier.*best,'-k','LineWidth',1.5,'LineSmoothing','off');
        hold off;
        set(gca,'xlim',[0 state.Generation]);
        xlabel('Generation');
        ylabel(PlotLabel);
        box on
%         legend('Best fitness','Mean fitness');
%         title(['Best: ',num2str(min(state.Score(:,1)),'%3.1f'),' Mean: ',num2str(meanf(state.Score),'%3.1f')])     
        title(['Best: ',num2str(-1.*min(state.Score(:,1)),'%3.1f')]);     
        %=================================================================%
        
        %=================== Plot Avg. Distance===========================%
        subplot(4,3,[9 12]);
        population = state.Population;
        Dist = 0;
        for i = 1:samples
            d = population(choices(i,1),:) - population(choices(i,2),:);
            Dist = Dist + sqrt( sum ( d.* d));
        end
        X = 0:1:state.Generation;
        distance = [distance Dist/samples];
        hold on;
        plot(X,distance,'.k','MarkerSize',15);
        plot(X,distance,'-k','LineWidth',1.5,'LineSmoothing','off');
        hold off;
        set(gca,'xlimmode','manual','zlimmode','manual','alimmode','manual');
        set(gca,'xlim',[0 state.Generation]);
        xlabel('Generation');
        ylabel('Average Distance');
        box on
        title('Average Distance Between Individuals')
        %=================================================================%  
    elseif StructuralOpt == 1
        cla
        xy = mogaCustomPlot(options,state,flag);
        xy(isnan(xy(:,1)),:)=[];  %get rid of any solutions which failed in some way (i.e. NaN)
        POpts = intersect([state.Score(:,1) state.Score(:,2)],xy,'rows');
        INFpts = setxor([state.Score(:,1) state.Score(:,2)],xy,'rows');
        INFpts(isnan(INFpts(:,1)),:)=[];
        hold on
        plot(Modifier.*POpts(:,1),POpts(:,2),'or','MarkerSize',4,'MarkerFaceColor','r','LineSmoothing','off');
        plot(Modifier.*INFpts(:,1),INFpts(:,2),'ok','MarkerSize',3,'LineSmoothing','off');
        hold off
        xlabel(PlotLabel);
        ylabel('Total Blade Mass (kg)');
        box on
        title(['Pareto Front, Generation = ' num2str(state.Generation)]);
        legend('Pareto Optimal','Inferior',2);   
    end
    
    case 'done'
        figure(GAFig);
        
        [BestValue,BestIndex] = min(state.Score(:,1));
        BestIndiv = state.Population(BestIndex,:);
        
        %===================== Plot Best Blade Geometry ==================%
        [ShapeError RElm TWIST CHORD PERCENT_THICKNESS DIMENSIONAL_THICKNESS...
         R_CHORD_CP CHORD_CP R_TWIST_CP TWIST_CP THICK_CP] = Define_Blade_Shape(BestIndiv);
    if StructuralOpt == 0
        subplot(4,3,[1 2]);
        hold on;
        plot(R_TWIST_CP,TWIST_CP,'sk',RElm,TWIST,'+-b');
        legend('Control Points','Pre-Twist');
        plot(RElm,TWIST,'-b','LineWidth',1.5,'LineSmoothing','off');
        ylabel('Pre-twist (deg)');
        xlabel('Blade Radius (m)');
        xlim([0 RotorRad]);
        box on
        hold off;
          if SpdCtrl == 0
            title(['Best Individual: Blade Geometry, Rotor Speed=' num2str(BestIndiv(end),'%-3.2f') ' rpm']);
          else
            title('Best Individual: Blade Geometry');
          end
        if CircleRoot == 1
          subplot(4,3,[4 5]);
          hold on;
          plot(R_CHORD_CP(1:3),CHORD_CP(1:3),'ok',R_CHORD_CP(4:end),CHORD_CP(4:end),'sk',RElm,CHORD,'+-r');
          legend('Circular Root Control Points','Control Points','Chord');
          plot(RElm,CHORD,'-r','LineWidth',1.5,'LineSmoothing','off');        
          ylabel('Chord (m)'); 
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]);
          box on
          hold off;
        else
          subplot(4,3,[4 5]);
          hold on;
          plot(R_CHORD_CP,CHORD_CP,'sk',RElm,CHORD,'+-r');
          legend('Control Points','Chord'); 
          plot(RElm,CHORD,'-r','LineWidth',1.5,'LineSmoothing','off');        
          ylabel('Chord (m)'); 
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]); 
          box on
          hold off;
        end
        if ThickMethod == 2
          subplot(4,3,[7 8]);
          hold on;
          plot(THICK_CP,Thickness_values,'sk');
          plot(RElm,PERCENT_THICKNESS,'+-','Color',[0 0.5 0]);
          legend('Control Points','% Thickness (t/c)');
          plot(RElm,PERCENT_THICKNESS,'-','Color',[0 0.5 0],'LineWidth',1.5,'LineSmoothing','off');          
          ylabel('Thickness (%)');
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]); 
          box on
          hold off;
        else
          subplot(4,3,[7 8]);         
          hold on;
            if length(Thickness_values) > 1
              tcp = zeros(1,length(Thickness_values));
              for k = 1:length(Thickness_values)-1;
                  tcp(k,1) = 0.5*(Thickness_values(k)+Thickness_values(k+1));
              end
            else
                tcp = Thickness_values(1);
            end
          plot(THICK_CP,tcp,'sk');    
          plot(RElm,PERCENT_THICKNESS,'+-','Color',[0 0.5 0]);
          legend('Control Points','% Thickness (t/c)'); 
          plot(RElm,PERCENT_THICKNESS,'-','Color',[0 0.5 0],'LineWidth',1.5,'LineSmoothing','off');          
          ylabel('Thickness (%)');
          xlabel('Blade Radius (m)');
          xlim([0 RotorRad]); 
          box on
          hold off;
        end
        subplot(4,3,[10 11]);
        hold on;
        plot(RElm,DIMENSIONAL_THICKNESS,'+-k');
        %legend('Dimensional Thickness (t)');
        plot(RElm,DIMENSIONAL_THICKNESS,'-k','LineWidth',1.5,'LineSmoothing','off');           
        ylabel('Thickness (m)');
        xlabel('Blade Radius (m)');
        xlim([0 RotorRad]);
        box on
        hold off;
        %=================================================================%
        
        %=================== Plot Best and Avg. Fitness ==================%
        subplot(4,3,[3 6]);
        X = 0:1:state.Generation;
%         plot(X,best,'.-k',X,m,'.b','LineWidth',1.5,'MarkerSize',15,'LineSmoothing','off');
        hold on;
        plot(X,Modifier.*best,'.k','MarkerSize',15);
        plot(X,Modifier.*best,'-k','LineWidth',1.5,'LineSmoothing','off');
        hold off;
        set(gca,'xlim',[0 state.Generation]);
        xlabel('Generation');
        ylabel(PlotLabel);
        box on
%         legend('Best fitness','Mean fitness');
        title(['Best: ',num2str(-1.*min(state.Score(:,1)),'%3.1f')])
        %=================================================================%
        
        %=================== Plot Avg. Distance===========================%
        subplot(4,3,[9 12]);
        X = 0:1:state.Generation;
        hold on;
        plot(X,distance,'.k','MarkerSize',15);
        plot(X,distance,'-k','LineWidth',1.5,'LineSmoothing','off');
        hold off;
        set(gca,'xlimmode','manual','zlimmode','manual','alimmode','manual');
        set(gca,'xlim',[0 state.Generation]);
        xlabel('Generation');
        ylabel('Average Distance');
        box on
        title('Average Distance Between Individuals')
        %=================================================================%
    elseif StructuralOpt == 1
        cla
        xy = mogaCustomPlot(options,state,flag);
        [state.Score(:,1) state.Score(:,2)];
        xy(isnan(xy(:,1)),:)=[]; %get rid of any solutions which failed in some way (i.e. NaN)
        POpts = intersect([state.Score(:,1) state.Score(:,2)],xy,'rows');
        INFpts = setxor([state.Score(:,1) state.Score(:,2)],xy,'rows');
        INFpts(isnan(INFpts(:,1)),:)=[];
        hold on
        plot(Modifier.*POpts(:,1),POpts(:,2),'or','MarkerSize',4,'MarkerFaceColor','r','LineSmoothing','off');
        plot(Modifier.*INFpts(:,1),INFpts(:,2),'ok','MarkerSize',3,'LineSmoothing','off');
        hold off
        xlabel(PlotLabel);
        ylabel('Total Blade Mass (kg)');
        box on
        title(['Pareto Front, Generation = ' num2str(state.Generation)]);
        legend('Pareto Optimal','Inferior',2);
    end
        
end


%-------- Used for Plot best and avg. fitness -----------------------------
function m = meanf(x)
nans = isnan(x);
x(nans) = 0;
n = sum(~nans);
n(n==0) = NaN; % prevent divideByZero warnings
% Sum up non-NaNs, and divide by the number of non-NaNs.
m = sum(x) ./ n;