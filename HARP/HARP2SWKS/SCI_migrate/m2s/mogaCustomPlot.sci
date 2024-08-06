
function [xy] = mogaCustomPlot(options,state,flag)

// Output variables initialisation (not found in input variables)
xy=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//this function is a modified version of GAPLOTPARETO.m created by: 

//   Copyright 2007 The MathWorks, Inc.
//   $Revision: 1.1.6.3 $  $Date: 2008/04/06 19:14:24 $

%v02 = %t;if ~~mtlb_tree(state,"Rank") then %v02 = size(mtlb_e(state,"Score"),2)<2;end;
if %v02 then
  title("Pareto Plot: not available","interp","none");
  xy = [];
  return;
end;

objectivesToPlot = [1,2];

lengths = cumsum(mtlb_e(options,"PopulationSize"),firstnonsingleton(mtlb_e(options,"PopulationSize")));
lengths = lengths(1:$-1);
starts = [1,1+lengths];
ends = mtlb_s(mtlb_a(starts,mtlb_e(options,"PopulationSize")),1);
subPops = [starts;ends];

[unused,c] = size(subPops);
for i = 1:c
  pop = subPops(:,i);
  range = pop(1):pop(2);
  xy = plotFirstFront(state.Score(range,:),state.Rank(range),(objectivesToPlot(:))');
end;

//---------------------------------
endfunction

function [plotData] = plotFirstFront(score,nonDomRank,objectivesToPlot)

// Output variables initialisation (not found in input variables)
plotData=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


objectivesToPlot = [1;2];
score = score(:,objectivesToPlot(1:2));

// Get individual from first front
minRank = 1;
xy = score(mtlb_logic(nonDomRank,"==",minRank),:);
plotData = xy;
endfunction
