function [index] = FindInd(Vector,Value)

// Output variables initialisation (not found in input variables)
index=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//This function finds the index of the closest value to Value within the vector Vector
//Inputs:    Vector: A column or row vector, the values within Vector MUST be unique (i.e. no repeated values)AND monotonically increasing or decresing
//           Value: The scalar value being searched for within Vector
// 
//Outputs:   index: scalar, index of closest value to Value within the vector Vector
// 
//Note:  If the value being searched for lies directly in the middle of its 
//two adjacent values, then it is determined to be closest to the larger
//of the two values if Value is positive, and closest to the smallest value
//if Value is negative.  If the value being searched for lies outside of the
//range of Vector, than index is equal to either 1 or the length of Vector

// Value being searched for is outside of range given

if mtlb_logic(mtlb_s(mtlb_e(Vector,$),mtlb_e(Vector,1)),">",0) then //if Vector is monotonically INCREASING
 if mtlb_logic(Value,"<=",mtlb_e(Vector,1)) then
   index = 1;
 elseif mtlb_logic(Value,">=",mtlb_e(Vector,$)) then
   index = max(size(Vector));
 else //Value being searched for is inside of range given
  distance = abs(mtlb_s(Vector,Value));
  %v0_3 = distance; // ! L.23: abs(mtlb_logic(distance,"==",min(%v0_3,firstnonsingleton(%v0_3)))) may be replaced by:
  // !    --> mtlb_logic(distance,"==",min(%v0_3,firstnonsingleton(%v0_3))) if mtlb_logic(distance,"==",min(%v0_3,firstnonsingleton(%v0_3))) is Real.
  indexes = mtlb_find(abs(mtlb_logic(distance,"==",min(%v0_3,firstnonsingleton(%v0_3)))));
 
  if mtlb_logic(Value,">=",0) then
    if max(size(indexes))==2 then
      index = indexes(2);
    else
      index = indexes(1);
    end;
  
  else
    if max(size(indexes))==2 then
      index = indexes(1);
    else
      index = indexes;
    end;
  end;
 end;
else //if Vector is monotonically DECREASING
 if mtlb_logic(Value,"<=",mtlb_e(Vector,$)) then
   index = max(size(Vector));
 elseif mtlb_logic(Value,">=",mtlb_e(Vector,1)) then
   index = 1;
 else //Value being searched for is inside of range given
  distance = abs(mtlb_s(Vector,Value));
  //indexes = find(distance==min(distance));
  %v1_3 = distance; [minValue,indexes] = min(%v1_3,firstnonsingleton(%v1_3));
 
  if mtlb_logic(Value,">=",0) then
    if max(size(indexes))==2 then
      index = indexes(2);
    else
      index = indexes(1);
    end;
  
  else
    if max(size(indexes))==2 then
      index = indexes(1);
    else
      index = indexes;
    end;
  end;
 end;
end;


//End of Function
endfunction
