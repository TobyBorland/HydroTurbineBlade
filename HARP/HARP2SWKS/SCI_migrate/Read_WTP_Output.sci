function [Output, WTP_Error] = Read_WTP_Output(fid,NumCol,CC)
// Reads data from the WT_Perf output file, also checks for WT_Perf output errors
// Input:  fid: the MATLAB file identifier which corresponds to the WT_Perf output file
//         NumCol: the number of output columns in the WT_Perf output file          
//         CC: this is a flag that tells this function if the WT_Perf outfile
//             was from a Parametric Analysis or a Combined Case analysis
// 
// Output: Output:    matrix, the first column is the first column of the
//                    WT_Perf output file (either the flow speeds or TSR's), and the remaining
//                    columns are the columns in the WT_Perf output file
//         WTP_Error: scalar, equal to 0 if no errors are detected in the
//                    output file, and equal to 1 if there are output errors

global NumCases

if CC == 1; // Read the output from a Combined Case Analysis output file
    A = textscan(fid,'% f %f %f %f %f %f %f %f %f %f','HeaderLines',8);
    Output(:,1) = A{1};
    Output(:,2) = A{5};
    
    // Checks for WT_Perf divergence & convergence errors
    if  any(isnan(Output(:,2))) | isempty(Output(:,2)) | any(Output(:,2)<-999) & any(Output(:,2)>-10001) | size(Output,1) < NumCases
        WTP_Error = 1;
                
    else
        WTP_Error = 0;
    end

elseif CC == 0; // Read the output from a Parametric Analysis output file
    
    dbl = zeros(1,3*NumCol+2);
    dbl(1:3:end) = 37; // ASCII character code for "%"
    dbl(2:3:end) = 102; // ASCII character code for "f"
    dbl(3:3:end) = 32; // ASCII character code for " " (empty space)
    chr = char(dbl); // creates a string with a variable number of "%f ", this allows the textscan function to read a variable number of Columns (NumCol)

    // A = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',10,'delimiter','\t','CollectOutput',1);
    A = textscan(fid,chr,'HeaderLines',10,'delimiter','\t','CollectOutput',1);
    Output = cell2mat(A);

    WTP_Error = zeros(1,NumCol);
    // Checks for WT_Perf divergence & convergence errors
    for n = 1:NumCol
    if any(any(isnan(Output(:,n+1)))) | isempty(Output(:,n+1)) |(any(any(Output(:,n+1)<-999)) & any(any(Output(:,n+1)>-10001)))
       WTP_Error(n) = 1;
    end
    
    end

end


// End of Function
endfunction
