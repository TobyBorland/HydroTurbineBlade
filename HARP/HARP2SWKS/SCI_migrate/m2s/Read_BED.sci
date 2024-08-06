function [BED_Error,BED,AoA,AFang,LocVel,Cl,Cd,CpMin,Thrust] = Read_BED(fid,NumCases,NumSeg,NumCol)

// Output variables initialisation (not found in input variables)
BED_Error=[];
BED=[];
AoA=[];
AFang=[];
LocVel=[];
Cl=[];
Cd=[];
CpMin=[];
Thrust=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This function reads the WT_Perf BED file and checks for Output errors

global("Type")

for n = 1:6 //moves cursor through the header of the file
 %v0 = mgetl(fid,1); if meof()~=0 then %v0 = -1;end; %v0;
end;

// ! L.10: real(NumCases .*NumCol) may be replaced by:
// !    --> NumCases .*NumCol if NumCases .*NumCol is Real.
BED = cell(real(NumCases .*NumCol),1);
BED_Error = zeros(1,NumCol);//initial dummy values
delim = mtlb_imp(NumCases,NumCases,NumCases*NumCol);
j = 1;

for n = mtlb_imp(1,NumCases*NumCol)
  last_err = 0;

  if mtlb_logic(Type,"==",2) then
    COL = 18;
    // !! L.20: Matlab function textscan not yet converted, original calling sequence used.
    Data = textscan(fid,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f","HeaderLines",4,"delimiter","\t","CollectOutput",0);
  elseif mtlb_logic(Type,"==",1) then
    COL = 17;
    // !! L.23: Matlab function textscan not yet converted, original calling sequence used.
    Data = textscan(fid,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f","HeaderLines",4,"delimiter","\t","CollectOutput",0);
  end;


  // !! L.27: Matlab function cellfun not yet converted, original calling sequence used.
  // ! L.27: abs(mtlb_logic(cellfun("length",Data),"~=",NumSeg)) may be replaced by:
  // !    --> mtlb_logic(cellfun("length",Data),"~=",NumSeg) if mtlb_logic(cellfun("length",Data),"~=",NumSeg) is Real.

  if mtlb_any(abs(mtlb_logic(cellfun("length",Data),"~=",NumSeg))) then //check if any of the cells are empty, this happens when WT_Perf diverges and the tries to output a number which is larger than the write field, resulting in """"****''s"""", and the textscan function wont read the whole line
   BED(n,1).entries = ones(NumSeg,18) .*9999.999;
  else
    BED(n,1).entries = cell2mat(Data);
  end;


  if n>delim(j) then
    j = j+1;
  end;
  //Check for errors in the BED file
  // ! L.38: abs(size(BED(n,1).entries)~=[NumSeg,COL]) may be replaced by:
  // !    --> size(BED(n,1).entries)~=[NumSeg,COL] if size(BED(n,1).entries)~=[NumSeg,COL] is Real.
  %v02 = %t;  if ~or(abs(size(BED(n,1).entries)~=[NumSeg,COL])) then %v02 = isempty(BED(n,1).entries);end;  // ! L.38: abs(mtlb_logic(BED(n,1).entries(:,3),">",9999)) may be replaced by:
  // !    --> mtlb_logic(BED(n,1).entries(:,3),">",9999) if mtlb_logic(BED(n,1).entries(:,3),">",9999) is Real.
  %v22 = %t;  if ~%v02 then %v22 = mtlb_any(abs(mtlb_logic(BED(n,1).entries(:,3),">",9999)));end;  // ! L.38: abs(isnan(BED(n,1).entries)) may be replaced by:
  // !    --> isnan(BED(n,1).entries) if isnan(BED(n,1).entries) is Real.
  // ! L.38: abs(mtlb_any(abs(isnan(BED(n,1).entries)))) may be replaced by:
  // !    --> mtlb_any(abs(isnan(BED(n,1).entries))) if mtlb_any(abs(isnan(BED(n,1).entries))) is Real.

  if %v22 | mtlb_any(abs(mtlb_any(abs(isnan(BED(n,1).entries))))) then
    BED_Error(1,j) = 1;
    last_err = 1;
  
  end;

  if last_err==1 then
    LocVel(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = zeros(NumSeg,1);
    CpMin = LocVel;
    Cl = LocVel;
    Cd = LocVel;
    Thrust = LocVel;
    AoA = LocVel;
    AFang = LocVel;
  else
    LocVel(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = BED(n,1).entries(:,3);
    CpMin(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = BED(n,1).entries(:,12);
    Cl(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = BED(n,1).entries(:,10);
    Cd(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = BED(n,1).entries(:,11);
    Thrust(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = BED(n,1).entries(:,16);
    AoA(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = BED(n,1).entries(:,9);
    AFang(mtlb_imp(mtlb_a((n-1)*NumSeg,1),NumSeg*n),1) = BED(n,1).entries(:,8);
  end;

end;

//End of Function
endfunction
