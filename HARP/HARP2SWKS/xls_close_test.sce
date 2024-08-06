HARP_dir = "C:\Documents and Settings\Foobert\Desktop\Geezer_Hydro\HARP\";
//dataFile =  fullfile(HARP_dir, "\Output_Files\", "\200W_BW3_FP_H2O\", "200W_BW3_FP_H2O_Output_XLS2003.xls");
//dataFile =  fullfile(HARP_dir, "\Output_Files\", "\1W5_BW3_1m_GRP3\", "1W5_BW3_1m_GRP3_Output_XLS3.xls");
dataFile =  fullfile(HARP_dir, "\Output_Files\", "\1W5_BW3_1m\", "1W5_BW3_1m_OutputXLS3.xls");



    [LU_xls,HO_SST,SheetN,SheetP] = xls_open(dataFile);
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
    disp(LU_xls);
    pause
    mclose();
