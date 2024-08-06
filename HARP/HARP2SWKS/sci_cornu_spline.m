function [x,y,spline_spec] = cornu_spline(spline_spec,t)
% Cornu Spline: two-point x, y data fit with slope constraints
%
% Requires erf_.m: erf_(x) returns the error function for a complex matrix x:
%   http://www.mathworks.com/matlabcentral/fileexchange/33577-complex-erf-error-function-fresnel-integrals
% cornu_spline also uses MATLAB's fsolve function.
%
% See demo example at the end of this file.
%
% syntax:
%   [x,y] = cornu_spline(spline_spec,t);
%   [~,~,spline_spec] = cornu_spline(spline_spec);
%   [x,y,spline_spec] = cornu_spline(spline_spec,t);
%
% inputs:
%
%   spline_spec: struct scalar, spline specification
%   spline_spec must contain the following data fields:
%
%     x0, y0: real scalar, spline's first end point
%
%     x1, y1: real scalar, spline's second end point
%
%     theta0, theta1: real scalar, tangency angles (radians) at first and second end
%     points, resp. (theta0 and theta1 are zero in the x-axis direction, positive toward
%     the y axis.)  theta0 and theta1 must be in the range
%       abs(theta0-Theta)+abs(theta1-Theta)<=1.5*pi
%     where
%       Theta = mod(atan2(y1-y0,x1-x0)-.5*(theta0+theta1)+pi,2*pi)+.5*(theta0+theta1)-pi
%
%     theta_ (optional): real scalar or [], curve angle at arc-length midpoint (t = 1/2)
%
%     L (optional): real scalar or [], curve's arc length between end points
%
%   The data fields theta_ and L may be initialized from a previous call to cornu_spline
%   (see outputs).  Re-use of the prior initialization can significantly reduce
%   computation time.
%
% optional input:
%
%   t: real matrix (optional, defaults to []): parameter values identifying points at
%   which spline points (x, y) are to be calculated.  t is proportional to arc length
%   from (x0,y0) and takes on the values 0 and 1 at points (x0,y0) and (x1,y1), resp.
%   The t values need not be limited to the range 0<=t<=1, although accuracy may be
%   compromized if they are very far outside of this range.
%
% outputs:
%
%   x, y: real matrices, size-matched to t
%   spline points corresponding to t
%
% optional output:
%
%   spline_spec: struct scalar, same as input, but with the theta_ and L data fields in
%   spline_spec initialized.  theta_ and L will be calculated if nargout>=3, or if
%   either value is initially unspecified or empty.
%
%
% Author: Kenneth C. Johnson, kjinnovation@earthlink.net, kjinnovation.com
% Version: November 2, 2011
%
% BSD Copyright notice:
%
% Copyright 2011 by Kenneth C. Johnson (kjinnovation.com)
%
% Redistribution and use in source and binary forms, with or without modification, are
% permitted provided that this entire copyright notice is duplicated in all such copies.
%
% This software is provided "as is" and without any express or implied warranties,
% including, without limitation, the implied warranties of merchantibility and fitness
% for any particular purpose.
%

% Documentation reference for "Eq" citations:
%   K. Johnson, "Cornu Spline," ver. Nov. 2, 2011
%
% subfunctions:  CSintegrals, FresnelCS

if ~isstruct(spline_spec) || ~isscalar(spline_spec)
    error('cornu_spiral:validation','spline_spec must be a length-1 struct.')
end
if ~all(isfield(spline_spec,{'x0','y0','x1','y1','theta0','theta1'}))
    error('cornu_spiral:validation',...
        'Validation error in spline_spec (x0, x0, x1, y1, theta0 or theta1).')
end
x0 = spline_spec.x0;
x1 = spline_spec.x1;
y0 = spline_spec.y0;
y1 = spline_spec.y1;
theta0 = spline_spec.theta0;
theta1 = spline_spec.theta1;
if ...
        ~isscalar(x0) || ~isnumeric(x0) || ~isreal(x0) || ...
        ~isscalar(y0) || ~isnumeric(y0) || ~isreal(y0) || ...
        ~isscalar(x1) || ~isnumeric(x1) || ~isreal(x1) || ...
        ~isscalar(y1) || ~isnumeric(y1) || ~isreal(y1) || ...
        ~isscalar(theta0) || ~isnumeric(theta0) || ~isreal(theta0) || ...
        ~isscalar(theta1) || ~isnumeric(theta1) || ~isreal(theta1)
    error('cornu_spiral:validation',...
        'Validation error in spline_spec (x0, x0, x1, y1, theta0 or theta1).')
end
theta_ = [];
L = [];
if nargout<3
    if isfield(spline_spec,'theta_')
        theta_ = spline_spec.theta_;
        if ~isempty(theta_) && ...
                (~isscalar(theta_) || ~isnumeric(theta_) || ~isreal(theta_))
            error('cornu_spiral:validation',...
                'Validation error in spline_spec (theta_).')
        end
    end
    if isfield(spline_spec,'L')
        L = spline_spec.L;
        if ~isempty(L) && ...
                (~isscalar(L) || ~isnumeric(L) || ~isreal(L))
            error('cornu_spiral:validation',...
                'Validation error in spline_spec (L).')
        end
    end
end
if nargin<2
    t = [];
end

if isempty(theta_) || isempty(L)
    Z = complex(x1-x0,y1-y0);
    absZ = abs(Z);
    Z = Z/absZ;
    theta_ = .5*(theta0+theta1);
    Theta = mod(angle(Z)-theta_+pi,2*pi)+theta_-pi; % Eq 34
    if abs(theta0-Theta)+abs(theta1-Theta)>1.5*pi % Eq 33
        error('cornu_spiral:validation','Invalid theta0, theta1.')
    end
    
    re_im = @(z) [real(z),imag(z)];
    %~ function check = re_im(z)
      %~ check = return([real(z) imag(z)]);
    %~ return
    %endfunction
    
    % thetaL = [theta_,L/absZ]:
    %~ function check = exp_fsolve(thetaL)
      %~ check = re_im(thetaL(2)*CSintegrals(theta0,thetaL(1),theta1,1)-Z) 
    %~ return

    
    [thetaL,~,exit_flag] =  fsolve(@(thetaL) ...
        re_im(thetaL(2)*CSintegrals(theta0,thetaL(1),theta1,1)-Z), ...
        [theta_,2],optimset('Display','off')); % Eq 31
	
	
    %~ [thetaL,check,exit_flag] =  fsolve(exp_fsolve(thetaL), ...
         %~ [theta_,2],optimset('Display','off')); % Eq 31

	
    theta_ = thetaL(1);
    L = thetaL(2)*absZ;
    if exit_flag==1 && (L<0 || ...
            abs(thetaL(2)*CSintegrals(theta0,theta_,theta1,1)-Z)>1e-4)
        exit_flag = -2;
    end
    if exit_flag~=1
        error('cornu_spiral:algorithm',...
            'Convergence failure. (There may be no solution.)')
    end
    spline_spec.theta_ = theta_;
    spline_spec.L = L;
end
if isempty(t)
    x = [];
    y = [];
    return
end
% Eq 30:
CS = CSintegrals(theta0,theta_,theta1,t);
x = L*real(CS)+x0;
y = L*imag(CS)+y0;
return

%~ function CS = CSintegrals(theta0,theta_,theta1,z)
%~ % Calculate the complex integral in Eq 36, with substitution from Eq's 27-28,
%~ % for scalars theta0, theta_, theta1 (theta values at t = 0, 1/2 and 1, resp., Eq's
%~ % 24-26).  z can be a matrix.
%~ dtheta0 = -3*theta0+4*theta_-theta1; % Eq 27
%~ ddtheta0 = 4*theta0-8*theta_+4*theta1; % Eq 28
%~ if abs(ddtheta0)^3>4*eps*dtheta0^2
    %~ % Eq's 37 and 44
    %~ sqrt_ = sqrt(abs(ddtheta0)/pi);
    %~ CS = (FresnelCS(sqrt_*(z+dtheta0/ddtheta0))-...
        %~ FresnelCS(sqrt_*dtheta0/ddtheta0))/sqrt_;
    %~ if ddtheta0<0
        %~ CS = conj(CS);
    %~ end
    %~ CS = exp(-.5*%i*dtheta0^2/ddtheta0+1*%i*theta0)*CS;
%~ else
    %~ % Eq's 38, 45
    %~ exp_itheta0 = exp(1*%i*theta0);
    %~ if dtheta0^2>2*eps
        %~ CS = ((exp(1*%i*dtheta0*z)-1)/(1*%i*dtheta0))*exp_itheta0; % Eq 46
    %~ else
        %~ CS = z*exp_itheta0; % Eq 47
    %~ end
    %~ if dtheta0^4>24*eps
        %~ % Eq's 45, 48
        %~ CS = CS+(.5*%i*ddtheta0)*(2/(1*%i*dtheta0)^3)*...
            %~ ((1-1*%i*dtheta0*z-.5*(dtheta0*z).^2).*exp(1*%i*dtheta0*z)-1)*exp_itheta0;
    %~ else
        %~ % Eq's 45, 49
        %~ CS = CS+(.5*%i*ddtheta0)*(z.^3/3)*exp_itheta0;
    %~ end
%~ end
%~ return

%~ function CS = FresnelCS(x)
%~ % Fresnel integrals (CS = C+1i*S)  (See erf_.m comment header.)
%~ CS = (.5+.5*%i)*erf_(sqrt(pi)*(.5-.5*%i)*x);
%~ return

%--------------------------------------------------------------------------------------
% Demo example:
figure, hold on
t = 0:.001:1;
clear spline_spec
spline_spec.x0 = -1;
spline_spec.y0 = 0;
spline_spec.x1 = 1;
spline_spec.y1 = 0;
Theta = atan2(spline_spec.y1-spline_spec.y0,spline_spec.x1-spline_spec.x0);
a=-1.5*pi; % a==theta0+theta1
for b = (-1.5:.3:1.5)*pi; % b==theta1-theta0
    spline_spec.theta0 = (a-b)/2+Theta;
    spline_spec.theta1 = (a+b)/2+Theta;
    [x,y] = cornu_spline(spline_spec,t);
    plot(x,y,'LineWidth',2), axis equal
end
set(gca,'XTick',[])
set(gca,'YTick',[])
axis([-2,2,-1.5,1.5])


