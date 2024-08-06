function [Alpha, CL_3D, CD_3D, C3_2D, C4_2D] = Apply_3D_Corrections(Aero_Coefficients,r_over_R,Chord,Vavg,Omg)
// This function applies the stall-delay equations to the 2D aerodynamic
// coefficients. This function only applies the stall-delay corrections for
// a single radial station at a time.
    
// Inputs:  Aero_Coefficients: a cell containing the entire 2D airfoil table, this
//                             cell can contain data for multiple Reynolds numbers as well.
//          r_over_R: scalar value, the non-dimensional radius
//          Chord:    scalar value, the local chord
//          Vavg:     scalar value, the average of the SpdSt and SpdEnd speeds
//          Omg:      scaler value, the rotor speed
// Outputs: Alpha: vector of the angles of attack, can include multiple Reynolds numbers
//          CL_3D: vector of the 3D corrected lift coefficients, can include multiple Reynolds numbers
//          CD_3D: vector of the 3D corrected drag coefficients, can include multiple Reynolds numbers
//          CP_2D: vector of the 2D minimum pressure coefficients, can include multiple Reynolds numbers

global Correct_3D RotorRad

// Separate Aerodynamic Coefficients Matrix into intuitive variables
    Alpha = Aero_Coefficients(:,1);
    CL_2D = Aero_Coefficients(:,2);
    CD_2D = Aero_Coefficients(:,3);
    
    // Checks to see if a third variable is included in the table (if there is, a new variable is created)
    if size(Aero_Coefficients,2) == 5;
      C3_2D = Aero_Coefficients(:,4);
      C4_2D = Aero_Coefficients(:,5);
    elseif size(Aero_Coefficients,2) == 4;
      C3_2D = Aero_Coefficients(:,4);
    end
//     if any(any(isnan(Aero_Coefficients(:,4)))); 
//         // don't do anything
//     else
//         CP_2D = Aero_Coefficients(:,4);
//     end

   
    
// // Apply 3D corrections to the lift coefficients (Selig-Du method)
// for all r/R <= 0.35, just use the 3D corrected values for r/R = 0.35.  This is
// done because the 3D Stall delay models act very strange for small r/R values.
if round(100*r_over_R) <= 35;
    r_over_R = 0.35;
end

d2r = pi/180; // converts degrees to radians
R_over_r = 1./r_over_R;
c_over_r = R_over_r.*Chord./RotorRad;

// modified tip speed ratio (from Selig-Du stall delay method)
mod_TSR = (Omg*RotorRad*pi/30)/sqrt(Vavg^2 + (Omg*RotorRad*pi/30)^2);

// calculate the lift slope from AoA = -2 to 4 degrees
x1 = -2;
x2 = 4;
y1 = interp1(Alpha,CL_2D,x1,'linear');
y2 = interp1(Alpha,CL_2D,x2,'linear');
Lift_Slope = (y2 - y1)/(x2 - x1); // units [1/deg]
CL_Yintercept = interp1(Alpha,CL_2D,0,'linear');
Alpha_zero = -CL_Yintercept/Lift_Slope; // units [deg]
Drag_zero = interp1(Alpha,CD_2D,0,'linear');

// Selig-Du scaling parameters (have been made a function of radius for better correlation)
// a = -11.847*(r_over_R)^3 + 15.184*(r_over_R)^2 - 4.6298*(r_over_R) + 1.0522;
a = 1;
// b = 0.8;
b = 1;
d = 1;
Alpha_End = 15; // units [deg] (may need to experiment with AlphaEnd value)


CL_3D = CL_2D; // initial dummy values
CD_3D = CD_2D; // initial dummy values

ii = FindInd(Alpha,90);
adj = zeros(ii,1);
DelCL = zeros(ii,1);
DelCD = zeros(ii,1);
// Make corrections only in the range of Alpha = -10 to 90 degrees
for n = FindInd(Alpha,-10):FindInd(Alpha,90);
    
  // Apply 3D corrections to the lift coefficients (Selig-Du method)
  // =======================================================================// 
    // Original definition of Cl_p defined by Selig-Du
    // CL_p = 2*pi*(Alpha(n)*pi/180 - Alpha_zero*pi/180);
    // Alternate definition of Cl_p used by AirfoilPrep
    CL_p = Lift_Slope*(Alpha(n) - Alpha_zero);
    
    Expon = d/(mod_TSR*r_over_R);
    // FL = ((1.6*c_over_r/0.1267)*(a - c_over_r^Expon)/(b + c_over_r^Expon) - 1)/(2*pi);
    // Alternate definition from AirfoilPrep divides by the lift slope
    FL = ((1.6*c_over_r/0.1267)*(a - c_over_r^Expon)/(b + c_over_r^Expon) - 1)/(Lift_Slope/d2r);
    
    if Alpha(n) < Alpha_End;
        CL_3D(n) = CL_2D(n) + FL*(CL_p - CL_2D(n));
    else
        // This adjustment factor was defined by C. Hansen (Windward Engineering)
        // for the AirfoilPrep worksheet
        adj(n) = ((90 - Alpha(n))/(90 - Alpha_End))^2;
        CL_3D(n) = CL_2D(n) + FL*(CL_p - CL_2D(n))*adj(n);
    end
  // =======================================================================// 
    
    if Correct_3D == 1
  // // Apply 3D corrections to the drag coefficients (Selig-Du method)
  // =======================================================================// 
        dExpon = Expon/2;   
        FD = ((1.6*c_over_r/0.1267)*(a - c_over_r^dExpon)/(b + c_over_r^dExpon) - 1)/(Lift_Slope/d2r);
        CD_3D(n) = CD_2D(n) - FD*(CD_2D(n) - Drag_zero);
        
        // Sometimes the stall model goes crazy at a small r/R values and
        // creates negative drags, this forces the drag to be zero if
        // negative
        if CD_3D(n) < 0
            CD_3D(n) = 0;
        end
  // =======================================================================// 
    else 
  // // Apply 3D corrections to the drag coefficients (Eggars method)
  // =======================================================================//     
        DelCL(n) = CL_3D(n) - CL_2D(n);
        DelCD(n) = DelCL(n) * (sin(Alpha(n)*d2r) - 0.12 * cos(Alpha(n)*d2r)) / (cos(Alpha(n)*d2r) + 0.12 * sin(Alpha(n)*d2r));
        CD_3D(n) = CD_2D(n) + DelCD(n);
        if CD_3D(n) < 0 // Make sure the drag coefficient never becomes negative, which sometimes happens when using these stall-delay models
            CD_3D(n) = 0;
        end
  // =======================================================================// 
    end
    
end

CL_3D(isnan(CL_3D)) = 0;
CD_3D(isnan(CD_3D)) = 0;

endfunction
// End of Function
