function [Strain, Iner, ShellArea] = getStrains(ST,STmax,CHORD,Mtan,Mnorm,Scaled_AF_Coordinates)

global NumAFPoints STdel SFstruct Einput

// For the inside shell, need to make the coordinates in CCW direction for
// the polygeom.m function to work properly
xx = flipud(Scaled_AF_Coordinates(:,1));
yy = flipud(Scaled_AF_Coordinates(:,2));

// separate into upper and lower coordinates
Xlow = xx(1:NumAFPoints/2);
Xup = xx(NumAFPoints-1:-1:NumAFPoints/2);
Ylow = yy(1:NumAFPoints/2);
Yup = yy(NumAFPoints-1:-1:NumAFPoints/2);

AngUp = atand(-1./FiniteDiff(Xup,Yup));
AngLow = atand(-1./FiniteDiff(Xlow,Ylow));
// Make sure surface normal vectors are pointing towards inside of airfoil
for n = 1:NumAFPoints/2
    if AngUp(n) > 0
       AngUp(n) = AngUp(n) - 180;
    end
    if AngLow(n) < 0
       AngLow(n) = AngLow(n) + 180;
    end
end

    if ST <= STmax & ST > (STmax-STdel)
        Xshell = Scaled_AF_Coordinates(:,1);    
        Yshell = Scaled_AF_Coordinates(:,2);
        // fprintf(1,'reached STmax\n');
    else
        // Calculate the inner x,y coordinates based on a shell thickness measured normal to the surface    
        InnerXup = Xup + ST.*cosd(AngUp);
        InnerYup = Yup + ST.*sind(AngUp);
        InnerXlow = Xlow + ST.*cosd(AngLow);
        InnerYlow = Ylow + ST.*sind(AngLow);

        // Find any intersections of upper and lower inside surfaces
        [Xo Yo I J] = intersections(InnerXup,InnerYup,InnerXlow,InnerYlow,1);

            if isempty(Xo)
                // Blade is completely solid material
                innerX = [];
                innerY = [];
            else
                innerX = [Xo(1);InnerXlow( ceil(J(1)):floor(J(end)) );Xo(end);...
                          InnerXup( floor(I(end)):-1:ceil(I(1)) );Xo(1)];
                innerY = [Yo(1);InnerYlow( ceil(J(1)):floor(J(end)) );Yo(end);...
                          InnerYup( floor(I(end)):-1:ceil(I(1)) );Yo(1)];        
            end

        Xshell = [Scaled_AF_Coordinates(:,1);innerX];
        Yshell = [Scaled_AF_Coordinates(:,2);innerY];
    end 
        
    // calculate the new values for Area, Inertia, and Centroids for the shell    
    [shellGEOM shellINER] = polygeom(Xshell,Yshell);
    Area = shellGEOM(1);
    Xcen_sh = shellGEOM(2);
    Ycen_sh = shellGEOM(3);
    Iuu_sh = shellINER(4);
    Ivv_sh = shellINER(5);
   
    // NOTE: for stress and strain calculations, the upper and lower points
    // only consider the normal moment component and the LE and TE only
    // consider the tangential moment component. In reality, both moment
    // components (normal AND tangential) will contribute to the
    // stress/strain at these points.
    
    // Calculate strains  
    cLE = Xcen_sh;
    cTE = CHORD - Xcen_sh;
    cUpper = abs(max(Scaled_AF_Coordinates(1:NumAFPoints/2,2)) - Ycen_sh);
    cLower = abs(min(Scaled_AF_Coordinates(NumAFPoints/2:NumAFPoints-1,2)) - Ycen_sh);
    StrainLE = (SFstruct.*Mtan.*1000.*cLE)./(1e9.*Einput.*Ivv_sh); // (m/m)
    StrainTE = (-SFstruct.*Mtan.*1000.*cTE)./(1e9.*Einput.*Ivv_sh); // (m/m)
    StrainUpper = (SFstruct.*Mnorm.*1000.*cUpper)./(1e9.*Einput.*Iuu_sh); // (m/m)
    StrainLower = (-SFstruct.*Mnorm.*1000.*cLower)./(1e9.*Einput.*Iuu_sh); // (m/m)
    
    // function outputs
    Strain.LE = StrainLE;
    Strain.TE = StrainTE;
    Strain.Upper = StrainUpper;
    Strain.Lower = StrainLower;
    Iner.uu = Iuu_sh;
    Iner.vv = Ivv_sh;
    ShellArea = Area;
    
endfunction
    
