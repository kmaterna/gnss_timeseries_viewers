function [alphaaz,e1,e2,salpha,se1,se2]=PrincipalStrains(e_phiphi,e_thetatheta,e_thetaphi,s_e_phiphi,s_e_thetatheta,s_e_thetaphi)
%  function [alphaaz,e1,e2,salpha,se1,se2]=PrincipalStrains(e_phiphi,e_thetatheta,e_thetaphi,s_e_phiphi,s_e_thetatheta,s_e_thetaphi)
%
% takes the outputs from strain_sphere to get principal strain (rates)
% alphaaz is the azimuth of the maximum extension direction clockwise from north.
% se1, se2 are the uncertainties in the principal strains
% salpha is uncertainty in alphaaz
%

A = [e_phiphi -e_thetaphi;
     -e_thetaphi e_thetatheta];
[V,D]=eig(A);
[e1,imax] = max(diag(D));
[e2,imin] = min(diag(D));
vmax = V(:,imax);

alpha = atan2(vmax(2),vmax(1))*180/pi;
alphaaz = 90-alpha; % make it into a geographic azimuth.

A = (e_phiphi - e_thetatheta).^2 + 4*e_thetaphi.^2;
dF_dethetaphi = -2*(e_phiphi- e_thetatheta)./A;
dF_dephiphi =  2*e_thetatheta./A;
dF_dethetatheta = -1*dF_dephiphi;
salpha = (180/(2*pi))*sqrt((s_e_thetaphi.*dF_dethetaphi).^2 + ...
                            (s_e_phiphi.*dF_dephiphi).^2 + ...
                            (s_e_thetatheta.*dF_dethetatheta).^2);

% compute uncertainties in the principal strain rates.
B = .5*(e_phiphi - e_thetatheta);
A = .5*(e_thetaphi.^2 + B.^2).^(-.5);
se1= sqrt( ... 
     (s_e_phiphi.*(.5 + A.*B)).^2 + ...
     (s_e_thetatheta.*(.5 - A.*B)).^2 + ...
     (s_e_thetaphi.*2.*A.*e_thetaphi).^2);
se2= sqrt( ... 
     (s_e_phiphi.*(.5 - A.*B)).^2 + ...
     (s_e_thetatheta.*(.5 + A.*B)).^2 + ...
     (s_e_thetaphi.*2.*A.*e_thetaphi).^2);
 
gamma = e1 - e2;
maxshear = gamma/2;
dilat = e1 + e2;
sdilat = sqrt(se1^2 + se2^2);