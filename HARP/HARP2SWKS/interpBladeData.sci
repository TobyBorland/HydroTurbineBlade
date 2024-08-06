function [rOverR, r, preTwist, chord, percT, t, pitchAxis] = interpBladeData(fname, numSects, interpType)
    // This program takes a tabluated blade design
    // and interpolates between it to create more airfoil sections
    
    // import data
    // this code assumes that the data is formatted in columns of
    // r/R, r, preTwist, chord, % thickness, thickness, pitch axis location
    
    [LU_xls,HO_SST,SheetN,SheetP] = xls_open(fname);
    // LU_xls: a number, the logical unit on the Excel stream.
    // HO_SST: vector of all character strings which appear in the Excel sheets.
    // SheetN: a vector of strings, the sheet names.
    // SheetP: a vector of numbers, the position of the beginning of sheets in the Excel stream.
    
    [data, iText] = xls_read(LU_xls, SheetP(1));
    // LU_xls: the logical unit on the Excel stream returned by xls_open.
    //SheetP: the position of the beginning of the sheet in the Excel stream.
    //        This position is one of those returned by xls_open.
    // data: numbers matrix, the numerical data found in the sheet. 
    //        Cells without numerical data are represented by NaN values.
    // iText: matrix of indices with the same size as Value. 
    //        The 0 indices indicates that no string exists in the corresponding Excel cell. 
    //        A positive index i points to the string SST(i) where SST is given by xls_open.
    mclose();
    
    // data = importdata(dataFile);
    offsetR = 7;
    offsetC = 22;
    offsetRend = 29;
//    while ~isnan(data(offsetR + offsetRend, 1 + offsetC))
//        offsetRend = offsetRend + 1;
//    end
//    offsetRend = offsetRend - 1;
//    disp(offsetRend);
    harpOpt = struct('r', data(offsetR:offsetRend, 2 + offsetC),...
                     'rOverR', data(offsetR:offsetRend, 1 + offsetC),...
                     'preTwist', data(offsetR:offsetRend, 3 + offsetC),...
                     'chord', data(offsetR:offsetRend, 4 + offsetC),...
                     'percT', data(offsetR:offsetRend, 5 + offsetC),...
                     't', data(offsetR:offsetRend, 6 + offsetC),...
                     'pitchAxis', data(offsetR:offsetRend, 7 + offsetC),...
                     'length', data(offsetR:offsetRend, 7 + offsetC));
                     
    radius = harpOpt.r(1)/harpOpt.rOverR(1); 
    harpOpt.length = radius - harpOpt.r(1);  
    harpOptEnd = harpOpt.r($);     
    r = harpOpt.r(1):harpOpt.length/numSects:harpOpt.r($);      
    
    disp(harpOpt);
    
    // interpolated between the airofil profiles data
    rOverR = r/radius;
    
    preTwist = interp1(harpOpt.r,harpOpt.preTwist,r,interpType);
    chord = interp1(harpOpt.r,harpOpt.chord,r,interpType);
    percT = interp1(harpOpt.r,harpOpt.percT,r,interpType);
    t = interp1(harpOpt.r,harpOpt.t,r,interpType);
    
    
    pitchAxis = interp1(harpOpt.r,harpOpt.pitchAxis,r,interpType);
    lengthB = harpOpt.length;
    disp(size(pitchAxis));
    disp(numSects);
    disp(".. interpolateBladeData() interp1() return");
endfunction
