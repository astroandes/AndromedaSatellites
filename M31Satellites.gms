$Title Simulación de satélites de M31

$Ontext
================================================================================
                         SEMINARIO DE ASTRONOMÍA
                                 2015-20

                         Simulación de satélites de M31

                         David Esteban Bernal Neira
================================================================================
$Offtext

Sets
         time      /0*3200/
         Sat       /III,V/
;

Scalar
         deltaT
         T "units: 9.8*10e8 yr/h (approx 1Gyr)"  /8/
         G "in internal units (Gadget units)"    /43007.1/
         DM31_obs "+19 -18 Conn et al 2012"      /779/

*From internal kinematics (http://arxiv.org/pdf/1205.6864v1.pdf)
         VM31_x_obs      /162.8/
         VM31_y_obs      /-117.2/

         VM31_z_obs      /-281.1/

         RascenM31
         DeclM31

         pi              /3.14159/

         m_MW            /100/
         oneoverm_MW

         zxz_alpha_rz "rotation angle in the trasformation to M31galactic coordinates (see Conn et al. 2013 for a description)"
         zxz_beta_rx "rotation angle in the trasformation to M31galactic coordinates"
         zxz_gamma_rz "rotation angle in the trasformation to M31galactic coordinates"

         tiny "Minimal value for comparison with the coordinates" /1e-6/


*Potential parameters M31 parameters (from Geehan et al 2006 - via Magda)
         e       /2.71828183/
         rbulge  /0.61/
         rdisk   /5.4/
         rhalo   /13.5/
         Mbulge  /2.86144/
         disksurface     /0.046/
         halodensity     /120794.0/
         criticaldensity /277.72e-10/

         Mhalo_MW        /91.36/
         rhalo_MW        /24.54/
         Mdisk_MW        /10.0/
         rdisk_MW        /6.65/
         bdisk_MW        /0.26/
         Mbulge_MW       /3.4/
         rbulge_MW       /0.7/
;

deltaT=T/card(time);

RascenM31=(00.+42/60+44.3/3600)*pi/12;
DeclM31=(41+16/60+09/3600)*pi/180;
oneoverm_MW=1/m_MW;

zxz_alpha_rz=(90-39.8)*pi/180;
zxz_beta_rx=(-77.5)*pi/180;
zxz_gamma_rz=(90)*pi/180;

Parameters
         Dist(Sat)
         /
         III     723
         V       742
         /

         errDistplus(Sat)
         /
         III     +18
         V       +21
         /

         errDistminus(Sat)
         /
         III     -24
         V       -22
         /

         rascenHH(Sat)
         /
         III     00
         V       01
         /

         rascenMM(Sat)
         /
         III     35
         V       10
         /

         rascenSS(Sat)
         /
         III     33.8
         V       17.1
         /

         declG(Sat)
         /
         III     36
         V       47
         /

         declMM(Sat)
         /
         III     29
         V       34
         /

         declSS(Sat)
         /
         III     52
         V       41
         /

         Vls(Sat)
         /
         III     -344.3
         V       -391.5
         /

         errVlsplus(Sat)
         /
         III     +1.7
         V       +2.7
         /

         errVlsminus(Sat)
         /
         III     -1.7
         V       -2.7
         /

         m0(Sat)
         /
         III     0.96
         V       2.6
         /

         PlaneYN(Sat)
         /
         III     1
         V       0
         /

         Vtany_param(Sat)
         /
         III     -330
         V       0
         /

         rascen(Sat)
         decl(Sat)

*getM31centeredCoordinates
         sdecz
         sdec(Sat)
         cdecz
         cdec(Sat)
         redif(Sat)
         sradif(Sat)
         cradif(Sat)
         denom(Sat)
** compute tangent plane coords
         xi(Sat)
         eta(Sat)
** costheta is the cos of the angular separation between the satelite and M31
         costheta(Sat)
         theta(Sat)
** equations from Conn et al 2013 to get x, y and z and the velocities in a cartesian M31-centered coordinate system.
         old_Xsat_M31(Sat)
         old_Ysat_M31(Sat)
         old_Zsat_M31(Sat)

         old_VXsat_M31(Sat)
         old_VYsat_M31(Sat)
         old_VZsat_M31(Sat)



;

rascen(Sat)=(rascenHH(Sat)+rascenMM(Sat)/60+rascenSS(Sat)/3600)*pi/12;
decl(Sat)=(declG(Sat)+declMM(Sat)/60+declSS(Sat)/3600)*pi/180

sdecz=sin

display decl;

*Función objetivo
Variable
         fobj
;

Equation
         eq_fobj
;

eq_fobj..        fobj=e=1;

Model M31Satellites /all/;
Options NLP=IPOPT;
Solve M31Satellites using NLP min fobj;
