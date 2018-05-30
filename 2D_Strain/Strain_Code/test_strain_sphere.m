% Testing strain_sphere from a simple example. 
% The first triangle of my rectangular test. 

phi=[-124, -123, -123]';
theta=[-50.0, -50.0, -49.0]';
u_phi=[0, 0, 0]';
u_theta=[-10, -0, -0]';
s_phi=[0.5, 0.5, 0.5]';
s_theta=[0.5, 0.5, 0.5]';
weight=1;
paramsel=0;


[e_phiphi,e_thetaphi,e_thetatheta,omega_r,U_theta,U_phi,...
          s_omega_r,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,s_U_theta,...
          s_U_phi,chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p,r_PHITHETA,u_phi_p,u_theta_p]=...
          strain_sphere(phi,theta,u_phi,u_theta,s_phi,s_theta,weight,paramsel);
e_phiphi
e_thetaphi
e_thetatheta
omega_r
U_theta
U_phi
s_omega_r
s_e_phiphi
s_e_thetaphi
s_e_thetatheta
s_U_theta
s_U_phi
chi2
OMEGA
THETA_p
PHI_p
s_OMEGA
s_THETA_p
s_PHI_p
r_PHITHETA
u_phi_p
u_theta_p
      
DispStrainRate(e_phiphi,e_thetaphi,e_thetatheta,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,...
    chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p);