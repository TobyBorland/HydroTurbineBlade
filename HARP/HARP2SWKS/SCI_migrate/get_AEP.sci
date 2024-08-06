function [AEP] = get_AEP(PvsV,pUvars)
// This function calculate the annual energy production
// Inputs:    PvsV: a two column vector, first column is the flow speeds, second column is the power values
//            pUvars: vector which contains the appropriate variables needed for the selected probability distribution, pUvars  is defined in Main.m
// 
// Outputs:   AEP: scalar, the AEP value

global user_pU_interp ProbDist

U = PvsV(:,1);
P = PvsV(:,2);

if ProbDist == 1;
    // Rayleigh Distribution
    U_mean = pUvars;
    pU = (pi./2).*(U./U_mean.^2).*exp(-(pi./4).*(U./U_mean).^2);
    
elseif ProbDist == 2;
    // Weibull Distribution
    k = pUvars(1);
    c = pUvars(2);    
    pU = (k./c).*((U./c).^(k-1)).*exp(-(U./c).^k);
    
elseif ProbDist == 3;
    // User Defined Distribution, this was already read in and interpolated in the HARP_Opt.m file
    // Make sure the Power curve and flow distribution vectors are the same length
    pU = user_pU_interp(1:length(U));
elseif ProbDist == 0;
    AEP = 0;
    return
end

// Average Power
Pavg = SimpInt(U,P.*pU,1,length(U));

// Annual Energy Production
Hours = 8760;
AEP = Pavg*Hours;

// End of Function
endfunction
