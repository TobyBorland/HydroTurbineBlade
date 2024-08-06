function xy = mogaCustomPlot(options,state,flag)
// this function is a modified version of GAPLOTPARETO.m created by: 

//    Copyright 2007 The MathWorks, Inc.
//    $Revision: 1.1.6.3 $  $Date: 2008/04/06 19:14:24 $

if ~isfield(state,'Rank') | size(state.Score,2) < 2
    title('Pareto Plot: not available','interp','none');
    xy = [];
    return;
end

objectivesToPlot = [1 2];

lengths = cumsum(options.PopulationSize);
lengths = lengths(1:($-1));
starts = [1, 1 + lengths];
ends = starts + options.PopulationSize - 1;
subPops = [starts;ends];

[unused,c] = size(subPops);
for i = 1:c
    pop = subPops(:,i);
    range = pop(1):pop(2);
    xy = plotFirstFront(state.Score(range,:),state.Rank(range),objectivesToPlot(:)');
end
endfunction
// ---------------------------------
function plotData = plotFirstFront(score,nonDomRank,objectivesToPlot)

    objectivesToPlot = [1;2];
    score = score(:,objectivesToPlot(1:2));

//  Get individual from first front
minRank = 1;
xy = score((nonDomRank == minRank),:);
plotData = xy;
endfunction
