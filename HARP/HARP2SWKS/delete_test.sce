//series_mod14 = floor((subdiv - hubRoot + 1)/14)modulo((subdiv - hubRoot + 1),14);
for i = 1:floor((subdiv - hubRoot + 1)/14)
    //swb_pathname = fullfile(saveDir, 'AF_'+string(i)+'.swb');
    //createSWKswb(swb_pathname,airfoils,i*14,(i*14)+13);
    disp([((i-1)*14)+1:(i*14)]);
    //pause
end
if modulo((subdiv - hubRoot + 1),14) > 0 then
    //i = i + 1;
    //swb_pathname = fullfile(saveDir, 'AF_'+string(i)+'.swb');
    //createSWKswb(swb_pathname,airfoils,i*14,(subdiv - hubRoot + 1));    
    pause
    disp([i*14:(subdiv - hubRoot + 1)]);
end

