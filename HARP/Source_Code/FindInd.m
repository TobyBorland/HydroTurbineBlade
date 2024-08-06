function [index] = FindInd(Vector,Value)
%This function finds the index of the closest value to Value within the vector Vector
%Inputs:    Vector: A column or row vector, the values within Vector MUST be unique (i.e. no repeated values)AND monotonically increasing or decresing
%           Value: The scalar value being searched for within Vector
%
%Outputs:   index: scalar, index of closest value to Value within the vector Vector
%
%Note:  If the value being searched for lies directly in the middle of its 
%two adjacent values, then it is determined to be closest to the larger
%of the two values if Value is positive, and closest to the smallest value
%if Value is negative.  If the value being searched for lies outside of the
%range of Vector, than index is equal to either 1 or the length of Vector

% Value being searched for is outside of range given
    
if Vector(end) - Vector(1) > 0; %if Vector is monotonically INCREASING
    if Value <= Vector(1)
        index = 1;
    elseif Value >= Vector(end)
        index = length(Vector);
    else %Value being searched for is inside of range given
        distance = abs(Vector - Value);
        indexes = find(distance==min(distance));

        if Value >= 0
            if (length(indexes) == 2)
            index = indexes(2);
            else
            index = indexes(1) ;
            end

        else
            if (length(indexes) == 2)
            index = indexes(1);
            else
            index = indexes ;
            end
        end
    end
else %if Vector is monotonically DECREASING
    if Value <= Vector(end)
        index = length(Vector);
    elseif Value >= Vector(1)
        index = 1;
    else %Value being searched for is inside of range given
        distance = abs(Vector - Value);
        %indexes = find(distance==min(distance));
        [minValue indexes] = min(distance);
        
        if Value >= 0
            if (length(indexes) == 2)
            index = indexes(2);
            else
            index = indexes(1) ;
            end

        else
            if (length(indexes) == 2)
            index = indexes(1);
            else
            index = indexes ;
            end
        end
    end    
end
    

%End of Function
