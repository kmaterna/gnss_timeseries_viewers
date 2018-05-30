function DispStrainRate(e_phiphi,e_thetaphi,e_thetatheta,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,...
    chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p)
%DispStrainRate(e_phiphi,e_thetaphi,e_thetatheta,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,...
%    chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p)
%
% prints the result of the strain rate calculation to the screen

lonP=PHI_p*180/pi;
latP=90-(THETA_p*180/pi);    
if (90-(THETA_p*180/pi))<0
   lonP=lonP-180;
   latP=-latP;
   OMEGA = -OMEGA;
end

slonP=s_PHI_p*180/pi;
slatP=s_THETA_p*180/pi;


if ~((e_phiphi==0 && e_thetaphi==0 && e_thetatheta==0) ...
        || (isnan(e_phiphi) && isnan(e_thetaphi) && isnan(e_thetatheta)))

    [alphaaz,e1,e2,salphaaz,se1,se2]=...
        PrincipalStrains(e_phiphi,e_thetatheta,e_thetaphi,s_e_phiphi,s_e_thetatheta,s_e_thetaphi);

    disp(['e1= ' num2str(e1)  ' +/- ' num2str(se1)]);
    disp(['e2= ' num2str(e2)  ' +/- ' num2str(se2)]);
    disp(['alphaaz= ' num2str(alphaaz)  ' +/- ' num2str(salphaaz)]);
    edelta= e1+e2;
    sedelta= sqrt(se1^2 + se2^2);
    emaxshear = (e1-e2)/2;
    semaxshear = sedelta/2;
    disp(['edelta= ' num2str(edelta)  ' +/- ' num2str(sedelta)]);
    disp(['emaxshear= ' num2str(emaxshear)  ' +/- ' num2str(semaxshear)]);
end;

disp(['lonP= ' num2str(lonP) '+/- ' num2str(slonP)]);
disp(['latP= ' num2str(latP) '+/- ' num2str(slatP)]);
disp(['omega= ' num2str(OMEGA) '+/- ' num2str(s_OMEGA) ]);
disp(['omega= ' num2str(OMEGA*1e6*180/pi) '+/-' num2str(s_OMEGA*1e6*180/pi) ' Degrees per million years']);
disp(['chi2perdof= ' num2str(chi2)]);


