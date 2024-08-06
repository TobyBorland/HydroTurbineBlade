function x = hcosspace(StartPoint, EndPoint, NumPoints, Flag)
%if Flag = 0, half cosing spacing, points are packed more dense at the beginning
%if Flag = 1, half cosing spacing, points are packed more dense at the end
%if Flag = 3, full cosing spacing, points are packed more dense at the beginning AND end

L = EndPoint - StartPoint;

if Flag == 0;
   Ang = linspace(pi,pi/2,NumPoints)';
   x = StartPoint + L*(1+cos(Ang));
   x(1) = StartPoint; %fixes rounding errors from floating point
   x(end) = EndPoint; %fixes rounding errors from floating point
end

if Flag == 1;
   Ang = linspace(pi/2,0,NumPoints)';
   x = StartPoint + L*cos(Ang);
   x(1) = StartPoint; %fixes rounding errors from floating point
   x(end) = EndPoint; %fixes rounding errors from floating point
end

if Flag == 3;
   Ang = linspace(pi,0,NumPoints)';
   x = StartPoint + L*(1+cos(Ang))/2;
   x(1) = StartPoint; %fixes rounding errors from floating point
   x(end) = EndPoint; %fixes rounding errors from floating point
end



