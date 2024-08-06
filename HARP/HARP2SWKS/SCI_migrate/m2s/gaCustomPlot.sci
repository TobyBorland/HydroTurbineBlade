
function [state] = gaCustomPlot(options,state,flag)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//options, state, and flag are delivered from the MATLAB genetic algorithm functions

global("RotorRad","ThickMethod","Thickness_values","CircleRoot","SpdCtrl","StructuralOpt","OptMethod")
// !! L.5: Matlab function persistent not yet converted, original calling sequence used.
persistent("best","m","distance","GAFig")

//used for plot avg. distance (same as MATLAB gaplotdistance.m)
samples = 20;
choices = ceil(mtlb_sum(mtlb_e(options,"PopulationSize"))*rand(samples,2));

if mtlb_logic(OptMethod,"==",0) then
  PlotLabel = "Performance Fitness Value";
  Modifier = 1;
else
  PlotLabel = "(penalized) Annual Energy Production (kWh/yr)";
  Modifier = -1;
end;


select flag
  case "init" then
    // !! L.21: Matlab function gcf not yet converted.
    GAFig = mtlb(gcf);
    // !! L.22: Matlab function figure not yet converted, original calling sequence used.
    figure(GAFig);
    // !! L.23: Matlab function set not yet converted, original calling sequence used.
    // L.23: (Warning name conflict: function name changed from set to %set).
    %set(GAFig,"color","white")

    %v1_1 = state.Score(:,1);    [BestValue,BestIndex] = mtlb_min(%v1_1,firstnonsingleton(%v1_1));
    BestIndiv = state.Population(BestIndex,:);

    //===================== Plot Best Blade Geometry ==================%

    // !! L.30: Unknown function Define_Blade_Shape not converted, original calling sequence used.
    [ShapeError,RElm,TWIST,CHORD,PERCENT_THICKNESS,DIMENSIONAL_THICKNESS,R_CHORD_CP,CHORD_CP,R_TWIST_CP,TWIST_CP,THICK_CP] = Define_Blade_Shape(BestIndiv);
    if mtlb_logic(StructuralOpt,"==",0) then
      subplot(4,3,[1,2]);
      set(gca(),"auto_clear","off");
      plot(R_TWIST_CP,TWIST_CP,"sk",RElm,TWIST,"+-b");
      // !! L.35: Matlab function legend not yet converted, original calling sequence used.
      legend("Control Points","Pre-Twist");
      plot(RElm,TWIST,"-b","LineWidth",1.5,"LineSmoothing","off");
      ylabel("Pre-twist (deg)");
      xlabel("Blade Radius (m)");
      // !! L.39: Matlab function xlim not yet converted, original calling sequence used.
      xlim([0,RotorRad]);
      %v3_1 = gca();  %v3_1.box = "on";
      set(gca(),"auto_clear","on");
      if mtlb_logic(SpdCtrl,"==",0) then
        title("Best Individual: Blade Geometry, Rotor Speed="+msprintf("%-3.2f",mtlb_e(BestIndiv,$))+" rpm");
      else
        title("Best Individual: Blade Geometry");
      end;
      if mtlb_logic(CircleRoot,"==",1) then
        subplot(4,3,[4,5]);
        set(gca(),"auto_clear","off");
        plot(mtlb_e(R_CHORD_CP,1:3),mtlb_e(CHORD_CP,1:3),"ok",mtlb_e(R_CHORD_CP,4:$),mtlb_e(CHORD_CP,4:$),"sk",RElm,CHORD,"+-r");
        // !! L.51: Matlab function legend not yet converted, original calling sequence used.
        legend("Circular Root Control Points","Control Points","Chord");
        plot(RElm,CHORD,"-r","LineWidth",1.5,"LineSmoothing","off");
        ylabel("Chord (m)");
        xlabel("Blade Radius (m)");
        // !! L.55: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v4_1 = gca();  %v4_1.box = "on";
        set(gca(),"auto_clear","on");
      else
        subplot(4,3,[4,5]);
        set(gca(),"auto_clear","off");
        plot(R_CHORD_CP,CHORD_CP,"sk",RElm,CHORD,"+-r");
        // !! L.62: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","Chord");
        plot(RElm,CHORD,"-r","LineWidth",1.5,"LineSmoothing","off");
        ylabel("Chord (m)");
        xlabel("Blade Radius (m)");
        // !! L.66: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v5_2 = gca();  %v5_2.box = "on";
        set(gca(),"auto_clear","on");
      end;
      if mtlb_logic(ThickMethod,"==",2) then
        subplot(4,3,[7,8]);
        set(gca(),"auto_clear","off");
        plot(THICK_CP,Thickness_values,"sk");
        plot(RElm,PERCENT_THICKNESS,"+-","Color",[0,0.5,0]);
        // !! L.75: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","% Thickness (t/c)");
        plot(RElm,PERCENT_THICKNESS,"-","Color",[0,0.5,0],"LineWidth",1.5,"LineSmoothing","off");
        ylabel("Thickness (%)");
        xlabel("Blade Radius (m)");
        // !! L.79: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v6_1 = gca();  %v6_1.box = "on";
        set(gca(),"auto_clear","on");
      else
        subplot(4,3,[7,8]);
        set(gca(),"auto_clear","off");
        if max(size(Thickness_values))>1 then
          tcp = zeros(1,max(size(Thickness_values)));
          for k = 1:max(size(Thickness_values))-1
            tcp(k,1) = 0.5*mtlb_a(mtlb_e(Thickness_values,k),mtlb_e(Thickness_values,k+1));
          end;
        else
          tcp = mtlb_e(Thickness_values,1);
        end;
        plot(THICK_CP,tcp,"sk");
        plot(RElm,PERCENT_THICKNESS,"+-","Color",[0,0.5,0]);
        // !! L.95: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","% Thickness (t/c)");
        plot(RElm,PERCENT_THICKNESS,"-","Color",[0,0.5,0],"LineWidth",1.5,"LineSmoothing","off");
        ylabel("Thickness (%)");
        xlabel("Blade Radius (m)");
        // !! L.99: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v7_2 = gca();  %v7_2.box = "on";
        set(gca(),"auto_clear","on");
      end;
      subplot(4,3,[10,11]);
      set(gca(),"auto_clear","off");
      plot(RElm,DIMENSIONAL_THICKNESS,"+-k");
      //legend(''Dimensional Thickness (t)'');
      plot(RElm,DIMENSIONAL_THICKNESS,"-k","LineWidth",1.5,"LineSmoothing","off");
      ylabel("Thickness (m)");
      xlabel("Blade Radius (m)");
      // !! L.110: Matlab function xlim not yet converted, original calling sequence used.
      xlim([0,RotorRad]);
      %v8_1 = gca();  %v8_1.box = "on";
      set(gca(),"auto_clear","on");
      //=================================================================%
    
      //=================== Plot Best and Avg. Fitness ==================%
      subplot(4,3,[3,6]);
      %v10_1 = state.Score(:,1);  best = mtlb_min(%v10_1,firstnonsingleton(%v10_1));
      //         m    = meanf(state.Score);
      //         plot(state.Generation,best,''.-k'',state.Generation,m,''.b'',''LineWidth'',1.5,''MarkerSize'',15,''LineSmoothing'',''off'');
      set(gca(),"auto_clear","off");
      plot(mtlb_e(state,"Generation"),Modifier .*best,".k","MarkerSize",15);
      plot(mtlb_e(state,"Generation"),Modifier .*best,"-k","LineWidth",1.5,"LineSmoothing","off");
      set(gca(),"auto_clear","on");
      // !! L.124: Matlab function gca not yet converted, original calling sequence used.
      // !! L.124: Matlab function set not yet converted, original calling sequence used.
      // L.124: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlim",[0,1]);
      // !! L.125: Matlab function gca not yet converted, original calling sequence used.
      // !! L.125: Matlab function set not yet converted, original calling sequence used.
      // L.125: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xtick",[0,1]);
      xlabel("Generation");
      ylabel(PlotLabel);
      %v11_1 = gca();  %v11_1.box = "on";
      //         legend(''Best fitness'',''Mean fitness'');
      title("Best: ");
      //=================================================================%
    
      //=================== Plot Avg. Distance===========================%
      subplot(4,3,[9,12]);
      population = mtlb_e(state,"Population");
      distance = 0;
      for i = 1:samples
        d = population(choices(i,1),:)-population(choices(i,2),:);
        distance = mtlb_a(distance,sqrt(sum(d .*d,firstnonsingleton(d .*d))));
      end;
      distance = distance/samples;
      set(gca(),"auto_clear","off");
      plot(mtlb_e(state,"Generation"),distance,".k","MarkerSize",15);
      plot(mtlb_e(state,"Generation"),distance,"-k","LineWidth",1.5,"LineSmoothing","off");
      set(gca(),"auto_clear","on");
      // !! L.146: Matlab function gca not yet converted, original calling sequence used.
      // !! L.146: Matlab function set not yet converted, original calling sequence used.
      // L.146: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlimmode","manual","zlimmode","manual","alimmode","manual");
      // !! L.147: Matlab function gca not yet converted, original calling sequence used.
      // !! L.147: Matlab function set not yet converted, original calling sequence used.
      // L.147: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlim",[0,1]);
      // !! L.148: Matlab function gca not yet converted, original calling sequence used.
      // !! L.148: Matlab function set not yet converted, original calling sequence used.
      // L.148: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xtick",[0,1]);
      xlabel("Generation","interp","none");
      ylabel("Average Distance");
      %v12_1 = gca();  %v12_1.box = "on";
      title("Average Distance Between Individuals")
      //=================================================================%
    elseif mtlb_logic(StructuralOpt,"==",1) then
      // !! L.155: Unknown function mogaCustomPlot not converted, original calling sequence used.
      xy = mogaCustomPlot(options,state,flag);
      xy(isnan(xy(:,1)),:) = [];  //get rid of any solutions which failed in some way (i.e. NaN)
      // !! L.157: Matlab function intersect not yet converted, original calling sequence used.
      POpts = intersect([state.Score(:,1),state.Score(:,2)],xy,"rows");
      // !! L.158: Matlab function setxor not yet converted, original calling sequence used.
      INFpts = setxor([state.Score(:,1),state.Score(:,2)],xy,"rows");
      INFpts(isnan(INFpts(:,1)),:) = [];
      set(gca(),"auto_clear","off")
      plot(Modifier .*POpts(:,1),POpts(:,2),"or","MarkerSize",4,"MarkerFaceColor","r","LineSmoothing","off");
      plot(Modifier .*INFpts(:,1),INFpts(:,2),"ok","MarkerSize",3,"LineSmoothing","off");
      set(gca(),"auto_clear","on")
      xlabel(PlotLabel);
      ylabel("Total Blade Mass (kg)");
      %v17_2 = gca();  %v17_2.box = "on";
      // !! L.167: string output can be different from Matlab num2str output.
      title("Pareto Front, Generation = "+string(mtlb_e(state,"Generation")));
      // !! L.168: Matlab function legend not yet converted, original calling sequence used.
      legend("Pareto Optimal","Inferior",2);
    end;

  case "iter" then
    // ! L.172: mtlb(GAFig) can be replaced by GAFig() or GAFig whether GAFig is an M-file or not.
    // !! L.172: Matlab function figure not yet converted, original calling sequence used.
    figure(mtlb(GAFig));

    %v19_2 = state.Score(:,1);    [BestValue,BestIndex] = mtlb_min(%v19_2,firstnonsingleton(%v19_2));
    BestIndiv = state.Population(BestIndex,:);

    //===================== Plot Best Blade Geometry ==================%

    // !! L.179: Unknown function Define_Blade_Shape not converted, original calling sequence used.
    [ShapeError,RElm,TWIST,CHORD,PERCENT_THICKNESS,DIMENSIONAL_THICKNESS,R_CHORD_CP,CHORD_CP,R_TWIST_CP,TWIST_CP,THICK_CP] = Define_Blade_Shape(BestIndiv);
    if mtlb_logic(StructuralOpt,"==",0) then
      subplot(4,3,[1,2]);
      set(gca(),"auto_clear","off");
      plot(R_TWIST_CP,TWIST_CP,"sk",RElm,TWIST,"+-b");
      // !! L.184: Matlab function legend not yet converted, original calling sequence used.
      legend("Control Points","Pre-Twist");
      plot(RElm,TWIST,"-b","LineWidth",1.5,"LineSmoothing","off");
      ylabel("Pre-twist (deg)");
      xlabel("Blade Radius (m)");
      // !! L.188: Matlab function xlim not yet converted, original calling sequence used.
      xlim([0,RotorRad]);
      %v21_1 = gca();  %v21_1.box = "on";
      set(gca(),"auto_clear","on");
      if mtlb_logic(SpdCtrl,"==",0) then
        title("Best Individual: Blade Geometry, Rotor Speed="+msprintf("%-3.2f",mtlb_e(BestIndiv,$))+" rpm");
      else
        title("Best Individual: Blade Geometry");
      end;
      if mtlb_logic(CircleRoot,"==",1) then
        subplot(4,3,[4,5]);
        set(gca(),"auto_clear","off");
        plot(mtlb_e(R_CHORD_CP,1:3),mtlb_e(CHORD_CP,1:3),"ok",mtlb_e(R_CHORD_CP,4:$),mtlb_e(CHORD_CP,4:$),"sk",RElm,CHORD,"+-r");
        // !! L.200: Matlab function legend not yet converted, original calling sequence used.
        legend("Circular Root Control Points","Control Points","Chord");
        plot(RElm,CHORD,"-r","LineWidth",1.5,"LineSmoothing","off");
        ylabel("Chord (m)");
        xlabel("Blade Radius (m)");
        // !! L.204: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v22_1 = gca();  %v22_1.box = "on";
        set(gca(),"auto_clear","on");
      else
        subplot(4,3,[4,5]);
        set(gca(),"auto_clear","off");
        plot(R_CHORD_CP,CHORD_CP,"sk",RElm,CHORD,"+-r");
        // !! L.211: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","Chord");
        plot(RElm,CHORD,"-r","LineWidth",1.5,"LineSmoothing","off");
        ylabel("Chord (m)");
        xlabel("Blade Radius (m)");
        // !! L.215: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v23_2 = gca();  %v23_2.box = "on";
        set(gca(),"auto_clear","on");
      end;
      if mtlb_logic(ThickMethod,"==",2) then
        subplot(4,3,[7,8]);
        set(gca(),"auto_clear","off");
        plot(THICK_CP,Thickness_values,"sk");
        plot(RElm,PERCENT_THICKNESS,"+-","Color",[0,0.5,0]);
        // !! L.224: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","% Thickness (t/c)");
        plot(RElm,PERCENT_THICKNESS,"-","Color",[0,0.5,0],"LineWidth",1.5,"LineSmoothing","off");
        ylabel("Thickness (%)");
        xlabel("Blade Radius (m)");
        // !! L.228: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v24_1 = gca();  %v24_1.box = "on";
        set(gca(),"auto_clear","on");
      else
        subplot(4,3,[7,8]);
        set(gca(),"auto_clear","off");
        if max(size(Thickness_values))>1 then
          tcp = zeros(1,max(size(Thickness_values)));
          for k = 1:max(size(Thickness_values))-1
            tcp(k,1) = 0.5*mtlb_a(mtlb_e(Thickness_values,k),mtlb_e(Thickness_values,k+1));
          end;
        else
          tcp = mtlb_e(Thickness_values,1);
        end;
        plot(THICK_CP,tcp,"sk");
        plot(RElm,PERCENT_THICKNESS,"+-","Color",[0,0.5,0]);
        // !! L.244: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","% Thickness (t/c)");
        plot(RElm,PERCENT_THICKNESS,"-","Color",[0,0.5,0],"LineWidth",1.5,"LineSmoothing","off");
        ylabel("Thickness (%)");
        xlabel("Blade Radius (m)");
        // !! L.248: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v25_2 = gca();  %v25_2.box = "on";
        set(gca(),"auto_clear","on");
      end;
      subplot(4,3,[10,11]);
      set(gca(),"auto_clear","off");
      plot(RElm,DIMENSIONAL_THICKNESS,"+-k");
      //legend(''Dimensional Thickness (t)'');
      plot(RElm,DIMENSIONAL_THICKNESS,"-k","LineWidth",1.5,"LineSmoothing","off");
      ylabel("Thickness (m)");
      xlabel("Blade Radius (m)");
      // !! L.259: Matlab function xlim not yet converted, original calling sequence used.
      xlim([0,RotorRad]);
      %v26_1 = gca();  %v26_1.box = "on";
      set(gca(),"auto_clear","on");
      //=================================================================%
    
      //=================== Plot Best and Avg. Fitness ==================%
      subplot(4,3,[3,6]);
      // ! L.266: mtlb(best) can be replaced by best() or best whether best is an M-file or not.
      %v28_1 = state.Score(:,1);  best = [mtlb(best),mtlb_min(%v28_1,firstnonsingleton(%v28_1))];
      //         m    = [m meanf(state.Score)];
      X = mtlb_imp(0,1,mtlb_e(state,"Generation"));
      //         plot(X,best,''.-k'',X,m,''.b'',''LineWidth'',1.5,''MarkerSize'',15,''LineSmoothing'',''off'');
      set(gca(),"auto_clear","off");
      plot(X,Modifier .*best,".k","MarkerSize",15);
      plot(X,Modifier .*best,"-k","LineWidth",1.5,"LineSmoothing","off");
      set(gca(),"auto_clear","on");
      // !! L.274: Matlab function gca not yet converted, original calling sequence used.
      // !! L.274: Matlab function set not yet converted, original calling sequence used.
      // L.274: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlim",[0,mtlb_e(state,"Generation")]);
      xlabel("Generation");
      ylabel(PlotLabel);
      %v29_1 = gca();  %v29_1.box = "on";
      //         legend(''Best fitness'',''Mean fitness'');
      //         title([''Best: '',num2str(min(state.Score(:,1)),''%3.1f''),'' Mean: '',num2str(meanf(state.Score),''%3.1f'')])     
      %v31_1 = state.Score(:,1);  title("Best: "+msprintf("%3.1f",-1 .*mtlb_min(%v31_1,firstnonsingleton(%v31_1))));
      //=================================================================%
    
      //=================== Plot Avg. Distance===========================%
      subplot(4,3,[9,12]);
      population = mtlb_e(state,"Population");
      Dist = 0;
      for i = 1:samples
        d = population(choices(i,1),:)-population(choices(i,2),:);
        Dist = mtlb_a(Dist,sqrt(sum(d .*d,firstnonsingleton(d .*d))));
      end;
      X = mtlb_imp(0,1,mtlb_e(state,"Generation"));
      // ! L.292: mtlb(distance) can be replaced by distance() or distance whether distance is an M-file or not.
      distance = [mtlb(distance),Dist/samples];
      set(gca(),"auto_clear","off");
      plot(X,distance,".k","MarkerSize",15);
      plot(X,distance,"-k","LineWidth",1.5,"LineSmoothing","off");
      set(gca(),"auto_clear","on");
      // !! L.297: Matlab function gca not yet converted, original calling sequence used.
      // !! L.297: Matlab function set not yet converted, original calling sequence used.
      // L.297: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlimmode","manual","zlimmode","manual","alimmode","manual");
      // !! L.298: Matlab function gca not yet converted, original calling sequence used.
      // !! L.298: Matlab function set not yet converted, original calling sequence used.
      // L.298: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlim",[0,mtlb_e(state,"Generation")]);
      xlabel("Generation");
      ylabel("Average Distance");
      %v32_1 = gca();  %v32_1.box = "on";
      title("Average Distance Between Individuals")
      //=================================================================%  
    elseif mtlb_logic(StructuralOpt,"==",1) then
      // !! L.305: All children will be deleted, no HandleVisibility property in Scilab graphics.
      %v33_2 = gca();  delete(%v33_2.children)
      // !! L.306: Unknown function mogaCustomPlot not converted, original calling sequence used.
      xy = mogaCustomPlot(options,state,flag);
      xy(isnan(xy(:,1)),:) = [];  //get rid of any solutions which failed in some way (i.e. NaN)
      // !! L.308: Matlab function intersect not yet converted, original calling sequence used.
      POpts = intersect([state.Score(:,1),state.Score(:,2)],xy,"rows");
      // !! L.309: Matlab function setxor not yet converted, original calling sequence used.
      INFpts = setxor([state.Score(:,1),state.Score(:,2)],xy,"rows");
      INFpts(isnan(INFpts(:,1)),:) = [];
      set(gca(),"auto_clear","off")
      plot(Modifier .*POpts(:,1),POpts(:,2),"or","MarkerSize",4,"MarkerFaceColor","r","LineSmoothing","off");
      plot(Modifier .*INFpts(:,1),INFpts(:,2),"ok","MarkerSize",3,"LineSmoothing","off");
      set(gca(),"auto_clear","on")
      xlabel(PlotLabel);
      ylabel("Total Blade Mass (kg)");
      %v38_2 = gca();  %v38_2.box = "on";
      // !! L.318: string output can be different from Matlab num2str output.
      title("Pareto Front, Generation = "+string(mtlb_e(state,"Generation")));
      // !! L.319: Matlab function legend not yet converted, original calling sequence used.
      legend("Pareto Optimal","Inferior",2);
    end;

  case "done" then
    // ! L.323: mtlb(GAFig) can be replaced by GAFig() or GAFig whether GAFig is an M-file or not.
    // !! L.323: Matlab function figure not yet converted, original calling sequence used.
    figure(mtlb(GAFig));

    %v40_3 = state.Score(:,1);    [BestValue,BestIndex] = mtlb_min(%v40_3,firstnonsingleton(%v40_3));
    BestIndiv = state.Population(BestIndex,:);

    //===================== Plot Best Blade Geometry ==================%

    // !! L.330: Unknown function Define_Blade_Shape not converted, original calling sequence used.
    [ShapeError,RElm,TWIST,CHORD,PERCENT_THICKNESS,DIMENSIONAL_THICKNESS,R_CHORD_CP,CHORD_CP,R_TWIST_CP,TWIST_CP,THICK_CP] = Define_Blade_Shape(BestIndiv);
    if mtlb_logic(StructuralOpt,"==",0) then
      subplot(4,3,[1,2]);
      set(gca(),"auto_clear","off");
      plot(R_TWIST_CP,TWIST_CP,"sk",RElm,TWIST,"+-b");
      // !! L.335: Matlab function legend not yet converted, original calling sequence used.
      legend("Control Points","Pre-Twist");
      plot(RElm,TWIST,"-b","LineWidth",1.5,"LineSmoothing","off");
      ylabel("Pre-twist (deg)");
      xlabel("Blade Radius (m)");
      // !! L.339: Matlab function xlim not yet converted, original calling sequence used.
      xlim([0,RotorRad]);
      %v42_1 = gca();  %v42_1.box = "on";
      set(gca(),"auto_clear","on");
      if mtlb_logic(SpdCtrl,"==",0) then
        title("Best Individual: Blade Geometry, Rotor Speed="+msprintf("%-3.2f",mtlb_e(BestIndiv,$))+" rpm");
      else
        title("Best Individual: Blade Geometry");
      end;
      if mtlb_logic(CircleRoot,"==",1) then
        subplot(4,3,[4,5]);
        set(gca(),"auto_clear","off");
        plot(mtlb_e(R_CHORD_CP,1:3),mtlb_e(CHORD_CP,1:3),"ok",mtlb_e(R_CHORD_CP,4:$),mtlb_e(CHORD_CP,4:$),"sk",RElm,CHORD,"+-r");
        // !! L.351: Matlab function legend not yet converted, original calling sequence used.
        legend("Circular Root Control Points","Control Points","Chord");
        plot(RElm,CHORD,"-r","LineWidth",1.5,"LineSmoothing","off");
        ylabel("Chord (m)");
        xlabel("Blade Radius (m)");
        // !! L.355: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v43_1 = gca();  %v43_1.box = "on";
        set(gca(),"auto_clear","on");
      else
        subplot(4,3,[4,5]);
        set(gca(),"auto_clear","off");
        plot(R_CHORD_CP,CHORD_CP,"sk",RElm,CHORD,"+-r");
        // !! L.362: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","Chord");
        plot(RElm,CHORD,"-r","LineWidth",1.5,"LineSmoothing","off");
        ylabel("Chord (m)");
        xlabel("Blade Radius (m)");
        // !! L.366: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v44_2 = gca();  %v44_2.box = "on";
        set(gca(),"auto_clear","on");
      end;
      if mtlb_logic(ThickMethod,"==",2) then
        subplot(4,3,[7,8]);
        set(gca(),"auto_clear","off");
        plot(THICK_CP,Thickness_values,"sk");
        plot(RElm,PERCENT_THICKNESS,"+-","Color",[0,0.5,0]);
        // !! L.375: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","% Thickness (t/c)");
        plot(RElm,PERCENT_THICKNESS,"-","Color",[0,0.5,0],"LineWidth",1.5,"LineSmoothing","off");
        ylabel("Thickness (%)");
        xlabel("Blade Radius (m)");
        // !! L.379: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v45_1 = gca();  %v45_1.box = "on";
        set(gca(),"auto_clear","on");
      else
        subplot(4,3,[7,8]);
        set(gca(),"auto_clear","off");
        if max(size(Thickness_values))>1 then
          tcp = zeros(1,max(size(Thickness_values)));
          for k = 1:max(size(Thickness_values))-1
            tcp(k,1) = 0.5*mtlb_a(mtlb_e(Thickness_values,k),mtlb_e(Thickness_values,k+1));
          end;
        else
          tcp = mtlb_e(Thickness_values,1);
        end;
        plot(THICK_CP,tcp,"sk");
        plot(RElm,PERCENT_THICKNESS,"+-","Color",[0,0.5,0]);
        // !! L.395: Matlab function legend not yet converted, original calling sequence used.
        legend("Control Points","% Thickness (t/c)");
        plot(RElm,PERCENT_THICKNESS,"-","Color",[0,0.5,0],"LineWidth",1.5,"LineSmoothing","off");
        ylabel("Thickness (%)");
        xlabel("Blade Radius (m)");
        // !! L.399: Matlab function xlim not yet converted, original calling sequence used.
        xlim([0,RotorRad]);
        %v46_2 = gca();  %v46_2.box = "on";
        set(gca(),"auto_clear","on");
      end;
      subplot(4,3,[10,11]);
      set(gca(),"auto_clear","off");
      plot(RElm,DIMENSIONAL_THICKNESS,"+-k");
      //legend(''Dimensional Thickness (t)'');
      plot(RElm,DIMENSIONAL_THICKNESS,"-k","LineWidth",1.5,"LineSmoothing","off");
      ylabel("Thickness (m)");
      xlabel("Blade Radius (m)");
      // !! L.410: Matlab function xlim not yet converted, original calling sequence used.
      xlim([0,RotorRad]);
      %v47_1 = gca();  %v47_1.box = "on";
      set(gca(),"auto_clear","on");
      //=================================================================%
    
      //=================== Plot Best and Avg. Fitness ==================%
      subplot(4,3,[3,6]);
      X = mtlb_imp(0,1,mtlb_e(state,"Generation"));
      //         plot(X,best,''.-k'',X,m,''.b'',''LineWidth'',1.5,''MarkerSize'',15,''LineSmoothing'',''off'');
      set(gca(),"auto_clear","off");
      // ! L.420: mtlb(best) can be replaced by best() or best whether best is an M-file or not.
      plot(X,Modifier .*mtlb(best),".k","MarkerSize",15);
      // ! L.421: mtlb(best) can be replaced by best() or best whether best is an M-file or not.
      plot(X,Modifier .*mtlb(best),"-k","LineWidth",1.5,"LineSmoothing","off");
      set(gca(),"auto_clear","on");
      // !! L.423: Matlab function gca not yet converted, original calling sequence used.
      // !! L.423: Matlab function set not yet converted, original calling sequence used.
      // L.423: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlim",[0,mtlb_e(state,"Generation")]);
      xlabel("Generation");
      ylabel(PlotLabel);
      %v48_1 = gca();  %v48_1.box = "on";
      //         legend(''Best fitness'',''Mean fitness'');
      %v50_1 = state.Score(:,1);  title("Best: "+msprintf("%3.1f",-1 .*mtlb_min(%v50_1,firstnonsingleton(%v50_1))))
      //=================================================================%
    
      //=================== Plot Avg. Distance===========================%
      subplot(4,3,[9,12]);
      X = mtlb_imp(0,1,mtlb_e(state,"Generation"));
      set(gca(),"auto_clear","off");
      // ! L.435: mtlb(distance) can be replaced by distance() or distance whether distance is an M-file or not.
      plot(X,mtlb(distance),".k","MarkerSize",15);
      // ! L.436: mtlb(distance) can be replaced by distance() or distance whether distance is an M-file or not.
      plot(X,mtlb(distance),"-k","LineWidth",1.5,"LineSmoothing","off");
      set(gca(),"auto_clear","on");
      // !! L.438: Matlab function gca not yet converted, original calling sequence used.
      // !! L.438: Matlab function set not yet converted, original calling sequence used.
      // L.438: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlimmode","manual","zlimmode","manual","alimmode","manual");
      // !! L.439: Matlab function gca not yet converted, original calling sequence used.
      // !! L.439: Matlab function set not yet converted, original calling sequence used.
      // L.439: (Warning name conflict: function name changed from set to %set).
      %set(gca(),"xlim",[0,mtlb_e(state,"Generation")]);
      xlabel("Generation");
      ylabel("Average Distance");
      %v51_1 = gca();  %v51_1.box = "on";
      title("Average Distance Between Individuals")
      //=================================================================%
    elseif mtlb_logic(StructuralOpt,"==",1) then
      // !! L.446: All children will be deleted, no HandleVisibility property in Scilab graphics.
      %v52_2 = gca();  delete(%v52_2.children)
      // !! L.447: Unknown function mogaCustomPlot not converted, original calling sequence used.
      xy = mogaCustomPlot(options,state,flag);
      [state.Score(:,1),state.Score(:,2)];
      xy(isnan(xy(:,1)),:) = [];  //get rid of any solutions which failed in some way (i.e. NaN)
      // !! L.450: Matlab function intersect not yet converted, original calling sequence used.
      POpts = intersect([state.Score(:,1),state.Score(:,2)],xy,"rows");
      // !! L.451: Matlab function setxor not yet converted, original calling sequence used.
      INFpts = setxor([state.Score(:,1),state.Score(:,2)],xy,"rows");
      INFpts(isnan(INFpts(:,1)),:) = [];
      set(gca(),"auto_clear","off")
      plot(Modifier .*POpts(:,1),POpts(:,2),"or","MarkerSize",4,"MarkerFaceColor","r","LineSmoothing","off");
      plot(Modifier .*INFpts(:,1),INFpts(:,2),"ok","MarkerSize",3,"LineSmoothing","off");
      set(gca(),"auto_clear","on")
      xlabel(PlotLabel);
      ylabel("Total Blade Mass (kg)");
      %v59_2 = gca();  %v59_2.box = "on";
      // !! L.460: string output can be different from Matlab num2str output.
      title("Pareto Front, Generation = "+string(mtlb_e(state,"Generation")));
      // !! L.461: Matlab function legend not yet converted, original calling sequence used.
      legend("Pareto Optimal","Inferior",2);
    end;

end;


//-------- Used for Plot best and avg. fitness -----------------------------
endfunction

function [m] = meanf(x)

// Output variables initialisation (not found in input variables)
m=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

nans = isnan(x);
x = mtlb_i(x,nans,0);
n = mtlb_sum(~nans);
n = mtlb_i(n,mtlb_logic(n,"==",0),%nan);// prevent divideByZero warnings
// Sum up non-NaNs, and divide by the number of non-NaNs.
m = mtlb_sum(x) ./n;
endfunction
