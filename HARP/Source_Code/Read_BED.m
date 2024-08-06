function [BED_Error,BED,AoA,AFang,LocVel,Cl,Cd,CpMin,Thrust] = Read_BED(fid,NumCases,NumSeg,NumCol)
%This function reads the WT_Perf BED file and checks for Output errors

global Type

for n = 1:6; %moves cursor through the header of the file
    fgetl(fid);
end

BED = cell(NumCases.*NumCol,1);
BED_Error = zeros(1,NumCol);      %initial dummy values
delim = NumCases:NumCases:NumCases*NumCol;
j = 1;

for n = 1:(NumCases*NumCol)
    last_err = 0;
    
    if Type == 2
    COL = 18;    
    Data = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',4,'delimiter','\t','CollectOutput',0);
    elseif Type == 1
    COL = 17;
    Data = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',4,'delimiter','\t','CollectOutput',0);    
    end
    
    
     if any(cellfun('length',Data)~=NumSeg) %check if any of the cells are empty, this happens when WT_Perf diverges and the tries to output a number which is larger than the write field, resulting in "****'s", and the textscan function wont read the whole line
         BED{n,1} = ones(NumSeg,18).*9999.999;
      else
         BED{n,1} = cell2mat(Data);
      end
    
    
    if n > delim(j);
      j = j+1;
    end
  %Check for errors in the BED file
    if any(size(BED{n,1})~=[NumSeg COL]) || isempty(BED{n,1}) || any(BED{n,1}(:,3)>9999) || any(any(isnan(BED{n,1})))
      BED_Error(1,j) = 1;
      last_err = 1;
      
    end
  
  if last_err == 1;
    LocVel((n-1)*NumSeg+1:NumSeg*(n),1) = zeros(NumSeg,1);
    CpMin = LocVel;
    Cl = LocVel;
    Cd = LocVel;
    Thrust = LocVel;
    AoA = LocVel;
    AFang = LocVel;
  else
    LocVel((n-1)*NumSeg+1:NumSeg*(n),1) = BED{n,1}(:,3);
    CpMin((n-1)*NumSeg+1:NumSeg*(n),1)  = BED{n,1}(:,12);
    Cl((n-1)*NumSeg+1:NumSeg*(n),1)  = BED{n,1}(:,10);
    Cd((n-1)*NumSeg+1:NumSeg*(n),1)  = BED{n,1}(:,11); 
    Thrust((n-1)*NumSeg+1:NumSeg*(n),1)  = BED{n,1}(:,16); 
    AoA((n-1)*NumSeg+1:NumSeg*(n),1)  = BED{n,1}(:,9); 
    AFang((n-1)*NumSeg+1:NumSeg*(n),1)  = BED{n,1}(:,8); 
  end
    
end

%End of Function