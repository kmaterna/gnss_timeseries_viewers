function [e_phiphi,e_thetaphi,e_thetatheta,omega_r,U_theta,U_phi,...
          s_omega_r,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,s_U_theta,...
          s_U_phi,chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p,r_PHITHETA,u_phi_p,u_theta_p]=...
          strain_sphere(phi,theta,u_phi,u_theta,s_phi,s_theta,weight,paramsel)
% [e_phiphi,e_thetaphi,e_thetatheta,omega_r,U_theta,U_phi,...
%          s_omega_r,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,s_U_theta,...
%          s_U_phi,chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p,r_PHITHETA,u_phi_p,u_theta_p]=...
%          strain_sphere(phi,theta,u_phi,u_theta,s_phi,s_theta,weight,paramsel)
%
% INPUTS:
%   phi is longitude, entered in DEGREES
%   theta is COlatitude, entered in DEGREES (note: same as latitude minus 90 degrees).
%   u_phi is velocity in phi direction, entered in meters per yr.
%   u_theta is velocity in theta direction, entered in meters per yr.
%   s_phi is uncertainty in velocity in phi direction, entered in meters per yr.
%   s_theta is uncertainty in velocity in theta direction, entered in meters per yr.
%   weight = 1 if you want to use uncertainties as weights in the inversion.
%   paramsel = 1 if you want to invert for solid body rotations only (not strains).
%   paramsel = 2 if you want to invert for strains only.
%   paramsel = anything else, if you want to invert for both
%
% OUTPUTS:
%   e_phiphi = extensional strain rate along phi direction
%   e_thetaphi = shear strain rate 
%   e_thetatheta = extensional strain rate along theta direction
%   omega_r = average rotation rate around vertical axis
%   U_phi = average translation rate phi-ward (meters!)            
%   U_theta = average translation rate theta-ward (meters!)          
%   s_e_phiphi = uncertainty in extensional strain rate along phi direction
%   s_e_thetaphi = uncertainty in shear strain rate 
%   s_e_thetatheta = uncertainty in extensional strain rate along theta direction
%   s_omega_r = uncertainty in average rotation rate along vertical axis
%   chi2 is the misfit between the strain rate estimate and the data
%   OMEGA  = total rotation rate associated with Euler pole
%   THETA_p = colatitude of Euler Pole
%   PHI_p = longitude of Euler Pole
%   s_OMEGA  = unc. in total rotation rate associated with Euler pole
%   s_THETA_p = unc. in colatitude of Euler Pole
%   s_PHI_p = unc. in longitude of Euler Pole
%   r_PHITHETA = is correllation coefficient between PHI & THETA for
%               plotting uncertainty ellipses in gmt.
%   u_phi_p, u_theta_p are the predicted velocities at each location
%
% note that theta is COlatitude AND that thetaward is Southing
%                    ^^
%
% see Savage et al., JGR October 2001, p.22,005 for derivation
%
% Version Oct. 12, 2004

if nargin == 7
    paramsel = 0;
end

% convert to radians
theta = theta*pi/180;
phi = phi*pi/180;

r0 = 6.378e6;  % mean equitorial Earth radius

theta_0 = mean(theta);
phi_0 = mean(phi);

[~,m]=size(u_phi); 
[n,m2]=size(u_theta); 
if (m~=1 || m2~=1) 
   error('u_phi and u_theta must be column vectors'); 
end;

% check for colinearity
colin=0;
dc = 90-theta;
Gc = [ones(length(phi),1) phi];
mc = lscov(Gc,dc);
resc = Gc*mc - dc;
if (max(abs(resc))<1e-5)
    colin=1;
end
if (all(abs(phi-mean(phi))<1e-5) || all(abs(theta-mean(theta))<1e-5))
    colin=2;
end
if (colin==1 || colin==2)
    disp(['Warning:  Points are collinear. Type=' num2str(colin)]);
end

%% Make the matrix needed for least squares inversion
d = reshape([u_phi';u_theta'],2*n,1);

G = [];
if (paramsel~=2)
    for i = 1:n
       del_phi = phi(i)-phi_0; 
       del_theta =  theta(i)-theta_0;

       G(2*i-1,1) = -r0;
       G(2*i-1,2) = -r0*cos(theta_0)*del_phi;
       G(2*i-1,3) =  r0*del_theta;
       if (paramsel ~=1)
           G(2*i-1,4) =  r0*sin(theta_0)*del_phi;  
           G(2*i-1,5) =  r0*del_theta;  
           G(2*i-1,6) =  0; 
       end;

       G(2*i,1) = -r0*cos(theta_0)*del_phi;
       G(2*i,2) =  r0;
       G(2*i,3) = -r0*sin(theta_0)*del_phi;
       if (paramsel ~=1) 
           G(2*i,4) =  0;
           G(2*i,5) =  r0*sin(theta_0)*del_phi;  
           G(2*i,6) =  r0*del_theta;  
       end;
    end;
else
    for i = 1:n
       del_phi = phi(i)-phi_0; 
       del_theta =  theta(i)-theta_0;

       G(2*i-1,1) =  r0*sin(theta_0)*del_phi;  
       G(2*i-1,2) =  r0*del_theta;  
       G(2*i-1,3) =  0; 

       G(2*i,1) =  0;
       G(2*i,2) =  r0*sin(theta_0)*del_phi;  
       G(2*i,3) =  r0*del_theta;  
    end;
end



% perform the inversion as in an overdetermined system

covd = diag(reshape([(s_phi.^2)';(s_theta.^2)'],2*n,1));

if (weight == 1)
   W = diag(1./diag(covd));
   M = (G'*W*G)\G'*W;
else
   M = (G'*G)\G';
end;
m = M*d;
dpred = G*m;

% the predicted u_phi_p, u_theta_p;
u_phi_p=dpred(1:2:end);
u_theta_p=dpred(2:2:end);

if paramsel==2
    omega_theta=NaN;
    omega_phi=NaN;
    omega_r=NaN;
    e_phiphi = m(1);
	e_thetaphi = m(2);
	e_thetatheta = m(3);
elseif (paramsel==1)
    omega_theta = m(1);
    omega_phi = m(2);
    omega_r = m(3);
    e_phiphi = NaN;
	e_thetaphi = NaN;
	e_thetatheta = NaN;
else
 	omega_theta = m(1);
    omega_phi = m(2);
    omega_r = m(3);
    e_phiphi = m(4);
	e_thetaphi = m(5);
	e_thetatheta = m(6);
end;

U_theta = r0*omega_phi;
U_phi = -r0*omega_theta;

covm = M*covd*M';

% obtain the uncertainties
if (paramsel==2)
	s_e_phiphi = sqrt(covm(1,1));
	s_e_thetaphi = sqrt(covm(2,2));
	s_e_thetatheta = sqrt(covm(3,3));
    s_U_phi = 0;
    s_U_theta = 0;
    s_omega_r = NaN;
    s_omega_phi = NaN;
    s_omega_theta = NaN;
elseif (paramsel==1)
    s_e_phiphi = 0;
	s_e_thetaphi = 0;
	s_e_thetatheta = 0;
    s_U_phi = r0*sqrt(covm(1,1));
    s_U_theta = r0*sqrt(covm(2,2));
    s_omega_r = sqrt(covm(3,3));
    s_omega_phi = sqrt(covm(2,2));
    s_omega_theta = sqrt(covm(1,1));
else
    s_U_phi = r0*sqrt(covm(1,1));
    s_U_theta = r0*sqrt(covm(2,2));
    s_omega_r = sqrt(covm(3,3));
    s_omega_phi = sqrt(covm(2,2));
    s_omega_theta = sqrt(covm(1,1));
    s_e_phiphi = sqrt(covm(4,4));
    s_e_thetaphi = sqrt(covm(5,5));
    s_e_thetatheta = sqrt(covm(6,6));
end

N=length(d);
if (paramsel == 1 || paramsel==2)
   
   chi2 = ( (((d-G*m)')/covd)*(d-G*m))/(N-3);

else
   
   chi2 = ( (((d-G*m)')/covd)*(d-G*m))/(N-6);

end


% Compute the Euler Vectors for the solid body rotation
% from (A7) of Savage et al., October 2001, JGR appendix
if paramsel~=2
    OMEGA = sqrt(omega_r.^2 + omega_phi.^2 + omega_theta.^2);
    y = omega_r*cos(theta_0) - omega_theta*sin(theta_0);
    THETA_p = acos(y./OMEGA);
    A = omega_r.*sin(theta_0).*sin(phi_0) + omega_theta.*cos(theta_0).*sin(phi_0) + omega_phi.*cos(phi_0);
    B = omega_r.*sin(theta_0).*cos(phi_0) + omega_theta.*cos(theta_0).*cos(phi_0) - omega_phi.*sin(phi_0);
    PHI_p=atan2(A,B);
    z = 1./sqrt(1 - (y./OMEGA)^2);

    %cast the problem of finding the uncertainties as a linear
    %inverse problem (linearized) of 
    J(1,1) = omega_r/OMEGA;
    J(1,2) = omega_phi/OMEGA;
    J(1,3) = omega_theta/OMEGA;

    J(2,1) = (sin(theta_0)/(A^2 + B^2))*(B*sin(phi_0) - A*cos(phi_0));
    J(2,2) = (1/(A^2 + B^2))*(B*cos(phi_0) + A*sin(phi_0));
    J(2,3) = (cos(theta_0)/(A^2 + B^2))*(B*sin(phi_0) - A*cos(phi_0));

    J(3,1) = -(z/OMEGA)*(cos(theta_0) - y*omega_r/OMEGA^2);
    J(3,2) = z*y*omega_phi/OMEGA^3;
    J(3,3) = (z/OMEGA)*(sin(theta_0) + y*omega_theta/OMEGA^2);

    covp = [covm(3,3) covm(3,2) covm(3,1);
            covm(3,2) covm(2,2) covm(1,2);
            covm(3,1) covm(1,2) covm(1,1)];
    covE = J*covp*J';

    s_OMEGA     = sqrt(covE(1,1));
    s_PHI_p   = sqrt(covE(2,2));
    s_THETA_p     = sqrt(covE(3,3));

%     s_PHITHETA = covE(2,3);
    r_PHITHETA = covE(2,3)/sqrt(covE(2,2)*covE(3,3));

else
    OMEGA=NaN;
    THETA_p=NaN;
    PHI_p=NaN;
    r_PHITHETA=NaN;
    s_OMEGA   = NaN;
    s_PHI_p   = NaN;
    s_THETA_p = NaN;
end;


