function [x] = hcosspace(StartPoint,EndPoint,NumPoints,Flag)

// Output variables initialisation (not found in input variables)
x=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//if Flag = 0, half cosing spacing, points are packed more dense at the beginning
//if Flag = 1, half cosing spacing, points are packed more dense at the end
//if Flag = 3, full cosing spacing, points are packed more dense at the beginning AND end

L = mtlb_s(EndPoint,StartPoint);

if mtlb_logic(Flag,"==",0) then
  Ang = mtlb_linspace(%pi,%pi/2,NumPoints)';
  x = mtlb_a(StartPoint,L*(1+cos(Ang)));
  x = mtlb_i(x,1,StartPoint);  //fixes rounding errors from floating point
  x = mtlb_i(x,$,EndPoint);  //fixes rounding errors from floating point
end;

if mtlb_logic(Flag,"==",1) then
  Ang = mtlb_linspace(%pi/2,0,NumPoints)';
  x = mtlb_a(StartPoint,L*cos(Ang));
  x = mtlb_i(x,1,StartPoint);  //fixes rounding errors from floating point
  x = mtlb_i(x,$,EndPoint);  //fixes rounding errors from floating point
end;

if mtlb_logic(Flag,"==",3) then
  Ang = mtlb_linspace(%pi,0,NumPoints)';
  x = mtlb_a(StartPoint,(L*(1+cos(Ang)))/2);
  x = mtlb_i(x,1,StartPoint);  //fixes rounding errors from floating point
  x = mtlb_i(x,$,EndPoint);  //fixes rounding errors from floating point
end;



endfunction
