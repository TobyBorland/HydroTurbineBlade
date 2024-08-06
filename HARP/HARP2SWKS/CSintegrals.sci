function CS = CSintegrals(theta0,theta_,theta1,z)
// Calculate the complex integral in Eq 36, with substitution from Eq's 27-28,
// for scalars theta0, theta_, theta1 (theta values at t = 0, 1/2 and 1, resp., Eq's
// 24-26).  z can be a matrix.
dtheta0 = -3*theta0+4*theta_-theta1; // Eq 27
ddtheta0 = 4*theta0-8*theta_+4*theta1; // Eq 28
if abs(ddtheta0)^3>4*eps*dtheta0^2
    // Eq's 37 and 44
    sqrt_ = sqrt(abs(ddtheta0)/pi);
    CS = (FresnelCS(sqrt_*(z+dtheta0/ddtheta0))-...
        FresnelCS(sqrt_*dtheta0/ddtheta0))/sqrt_;
    if ddtheta0<0
        CS = conj(CS);
    end
    CS = exp(-.5*%i*dtheta0^2/ddtheta0+1*%i*theta0)*CS;
else
    // Eq's 38, 45
    exp_itheta0 = exp(1*%i*theta0);
    if dtheta0^2>2*eps
        CS = ((exp(1*%i*dtheta0*z)-1)/(1*%i*dtheta0))*exp_itheta0; // Eq 46
    else
        CS = z*exp_itheta0; // Eq 47
    end
    if dtheta0^4>24*eps
        // Eq's 45, 48
        CS = CS+(.5*%i*ddtheta0)*(2/(1*%i*dtheta0)^3)*...
            ((1-1*%i*dtheta0*z-.5*(dtheta0*z).^2).*exp(1*%i*dtheta0*z)-1)*exp_itheta0;
    else
        // Eq's 45, 49
        CS = CS+(.5*%i*ddtheta0)*(z.^3/3)*exp_itheta0;
    end
end
return
