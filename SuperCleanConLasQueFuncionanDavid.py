import numpy as np 			# paquete numerico de python
import scipy as sp
import pdb 				# python debugger
import matplotlib.pyplot as plt		# paquete para graficar
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D


DeltaT        = 0.0025
T	      = 8. #5.                     # units: 9.8*10e8 yr/h (approx 1Gyr , according to Geraint and Magda 0.98=1).
G 	      = 43007.1 		# in internal units (Gadget units) 
DM31_obs      = 779. 			# +19 -18 Conn et al 2012

# From internal kinematics (http://arxiv.org/pdf/1205.6864v1.pdf)
VM31_x_obs    =  162.8  
VM31_y_obs    =  -117.2

## From satellite LOS kinematics (http://arxiv.org/pdf/1205.6864v1.pdf)
#VM31_x_obs    =  -97  
#VM31_y_obs    =  -45.1 

# Weighted average (http://arxiv.org/pdf/1205.6864v1.pdf)
#VM31_x_obs    =  125.2 
#VM31_y_obs    =  -73.8

# Ignoring M31 tangential velocity
#VM31_x_obs    =  0.
#VM31_y_obs    =  0.

VM31_z_obs    = - 281.1

RascenM31     = (00. + 42./60 + 44.3/3600.)* np.pi / 12. 
DeclM31       = (41.+ 16./60. + 09./3600.)*np.pi/180 
m_MW          =	100.	
oneoverm_MW   = 1./m_MW
#SatSatInter   = 'SSI'                   # satellite-satellite interaction 
SatSatInter   = 'NSI'                   # NO satellite-satellite interaction   
MW_POT        = 'MW'                    # Includes MW potential   
#MW_POT        = 'NoMW'                  # does NOT include MW potential
NumIntegr     = 'Fwd'                   # Forward in time
#NumIntegr     = 'Bwd'                   # Backwards in time

TestSatName   = 'All'


zxz_alpha_rz = (90.-39.8)*np.pi/180     # rotation angle in the trasformation to M31galactic coordinates (see Conn et al. 2013 for a description)
zxz_beta_rx  = -77.5*np.pi/180          # rotation angle in the trasformation to M31galactic coordinates
zxz_gamma_rz = 90.0*np.pi/180           # rotation angle in the trasformation to M31galactic coordinates



class body(object):
	def __init__( self, Dist, errDistplus , errDistminus , rascenHH , rascenMM , rascenSS , declG , declMM , declSS , Vls , errVlsplus , errVlsminus , m0, name, PlaneYN, Vtany_param ):


		# mass
		# input m0 in the initialization is in 10^7 solar masses
 
		m0 =  m0 * 1e-3 #convert imput mass to 10^10 solar masses 
		#self.m	= 10. * m0  #assuming total mass is 10 times literature stellar mass
		self.m	= 0.5  #test mass
		self.oneoverm = 1/self.m
		# Satellite name	
		self.name = name
		# satelite belongs to the plane Yes/No/NA
		self.PlaneYN = PlaneYN
		# Line of sight velocity	
		self.Vls = Vls 
	
		# convert declination and right ascension into radians
		rascen    = (rascenHH + rascenMM/60 + rascenSS/3600)* np.pi / 12.
		decl  =  (declG + declMM/60. + declSS/3600.)*np.pi/180  
		print decl


		# calculate x, y, z, Vx, Vy, Vz in a M31-centered cartesian coordinates from the right ascension, declination and distance of the satellite. 
		old_Xsat_M31, old_Ysat_M31, old_Zsat_M31, old_VXsat_M31, old_VYsat_M31, old_VZsat_M31, eta, xi = getM31centeredCoorndinates(decl, rascen, Dist, Vls)


		# store the line of sight velocity before rotation and WHITOUT substracting M31 velocity
		self.Vls_xbr = old_VXsat_M31 
		self.Vls_ybr = old_VYsat_M31
		self.Vls_zbr = old_VZsat_M31 

		# rotation of the cartesian M31-centered coordinate system to get an M31galactic coordinate system (M31 disk in the xy plane). 
		# (Conn et al 2013)

		# 1) position
		#Xsat_M31, Ysat_M31, Zsat_M31 = rotateZXZ( old_Xsat_M31, old_Ysat_M31, old_Zsat_M31, zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
		Xsat_M31, Ysat_M31, Zsat_M31 = rotateZXZ( old_Xsat_M31, old_Ysat_M31, old_Zsat_M31, zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
		self.x = Xsat_M31
		self.y = Ysat_M31
		self.z = Zsat_M31
		R= np.sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
		self.R = R

		# 2) Line of sight velocity
		VXsat_M31, VYsat_M31, VZsat_M31 = rotateZXZ( old_VXsat_M31, old_VYsat_M31, old_VZsat_M31, zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
		#Satellite velocity: only line of sight velocity
		self.Vls_x = VXsat_M31 
		self.Vls_y = VYsat_M31
		self.Vls_z = VZsat_M31
		Vls= np.sqrt(self.Vls_x*self.Vls_x+self.Vls_y*self.Vls_y+self.Vls_z*self.Vls_z)
		self.Vls = Vls

		# 2.a) Construct vector perpendicular to the plane that contains the satellite, M31 and the MW.
		old_x_MW = 0.
		old_y_MW = 0.	
		old_z_MW = - DM31_obs
		
		x_MW0, y_MW0, z_MW0 = rotateZXZ( old_x_MW, old_y_MW, old_z_MW, zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
		#vector sat-M31
		x_satM31 = -self.x  #self.Vls_x
		y_satM31 = -self.y  #self.Vls_y   
		z_satM31 = -self.z  #self.Vls_z   
		#vector sat-MW
		x_satMW = x_MW0-self.x	    #unit_V_pl_x
		y_satMW = y_MW0-self.y      #unit_V_pl_y
		z_satMW = z_MW0-self.z      #unit_V_pl_z

		Vec_satM31 = np.array([x_satM31, y_satM31, z_satM31])
		Vec_satMW = np.array([x_satMW, y_satMW, z_satMW])

		Cross_VSM31VSMW = np.cross(Vec_satM31,Vec_satMW)
		#array([-1,  0,  0])
		#print "cross product", Cross_VSM31VSMW
		
		MagCross_VSM31VSMW = np.sqrt(np.dot(Cross_VSM31VSMW,Cross_VSM31VSMW))
		unit_plane_x = Cross_VSM31VSMW[0]/MagCross_VSM31VSMW
		unit_plane_y = Cross_VSM31VSMW[1]/MagCross_VSM31VSMW
		unit_plane_z = Cross_VSM31VSMW[2]/MagCross_VSM31VSMW
		Mag_unit_plane = unit_plane_x*unit_plane_x + unit_plane_y*unit_plane_y + unit_plane_z*unit_plane_z
				  
		VM31x, VM31y, VM31z = rotateZXZ( VM31_x_obs , VM31_y_obs , VM31_z_obs,  zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
				  
		# 3) Construct tangential velocity
		# conditions: 
		# a) Vtan perpendicular to the line of sight velocity
		# b) Vtot (Vls + Vtan) contained in the plane of satellites
				  
		# Vector perpendicular to the plane that contains M31, MW and the satellite:
		unit_V_pl_x = unit_plane_x 
		unit_V_pl_y = unit_plane_y 
		unit_V_pl_z = unit_plane_z 

		# step 1: choose one of the components of Vtan
		#print "self.PlaneYN", self.PlaneYN, self.name 
		if self.PlaneYN == "Yes": 
			Vtan_y= Vtany_param 
			if self.name == TestSatName: 
				Vtan_y= Vtany_param #- 98   
			#print "Vtan_y, Vtany_param", Vtan_y, Vtany_param 
		else:	
			#auxVtan_y= np.random.uniform(-1.732*self.Vls, 1.732*self.Vls, size=1) 
			auxVtan_y= np.random.uniform(-0.1*self.Vls, 0.1*self.Vls, size=1) 
			Vtan_y= auxVtan_y[0] 
		# step 2: construct the other 2 components so that Vtan satisfies conditions a) and b)
		
		auxVtan= self.Vls_x/(self.Vls_z*unit_V_pl_x-self.Vls_x*unit_V_pl_z)
		auxVtan2= self.Vls_x*unit_V_pl_x + self.Vls_y*unit_V_pl_y + self.Vls_z*unit_V_pl_z - VM31x*unit_V_pl_x - VM31y*unit_V_pl_y - VM31z*unit_V_pl_z -  self.Vls_y*Vtan_y*unit_V_pl_x/self.Vls_x + Vtan_y*unit_V_pl_y
		
		Vtan_z = auxVtan * auxVtan2		
		
		Vtan_x= -((self.Vls_y*Vtan_y + self.Vls_z*Vtan_z)/self.Vls_x)
		
		self.Vtan_x = Vtan_x
		self.Vtan_y = Vtan_y
		self.Vtan_z = Vtan_z
		Vtan = np.sqrt(Vtan_x*Vtan_x+Vtan_y*Vtan_y+Vtan_z*Vtan_z)
		self.Vtan = Vtan

		#M31 velocity after rotation
		VM31x, VM31y, VM31z = rotateZXZ( VM31_x_obs , VM31_y_obs , VM31_z_obs,  zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
		
		#Satellite velocity: line of sight velocity + tangential component - M31 velocity after rotation
		self.Vx = self.Vls_x+ self.Vtan_x - VM31x
		self.Vy = self.Vls_y+ self.Vtan_y - VM31y
		self.Vz = self.Vls_z+ self.Vtan_z - VM31z
		Vel= np.sqrt(self.Vx*self.Vx+self.Vy*self.Vy+self.Vz*self.Vz)
		self.Vel = Vel

		#test compute the component of the velocity vector in the unit_V_ls direction
		oneoverV_ls = 1/ Vls
		unit_V_ls_x = self.Vls_x*oneoverV_ls 
		unit_V_ls_y = self.Vls_y*oneoverV_ls 
		unit_V_ls_z = self.Vls_z*oneoverV_ls 

		self.StoreSpatialCoordsInit ()

	def StoreSpatialCoordsInit (self):
		#stores spatial coordinates			
		self.data_x = []
		self.data_y = []
		self.data_z = []
		self.data_Vx = []
		self.data_Vy = []
		self.data_Vz = []
		self.data_R = []
		self.data_Vel = []
		self.data_distPlane = []
		
		self.data_xbr = []
		self.data_ybr = []
		self.data_zbr = []
		self.data_Vxbr = []
		self.data_Vybr = []
		self.data_Vzbr = []
		
		self.data_x.append( self.x )
		self.data_y.append( self.y )
		self.data_z.append( self.z )
		self.data_R.append( self.R )
		self.data_Vx.append( self.Vx )
		self.data_Vy.append( self.Vy )
		self.data_Vz.append( self.Vz )
		self.data_Vel.append( self.Vel )
		#distance to the plane of satellites Conn et al. 2013
		self.distPlane = 0.158*self.x+0.769*self.y+0.620*self.z
		#print "self.distPlane", self.distPlane
		self.data_distPlane.append( self.distPlane )
		
		self.xbr, self.ybr, self.zbr = rotateZXZ( self.x,self.y,self.z, -zxz_gamma_rz, - zxz_beta_rx, -zxz_alpha_rz )
		self.Vxbr, self.Vybr, self.Vzbr = rotateZXZ(self.Vx,self.Vy,self.Vz, -zxz_gamma_rz, - zxz_beta_rx, -zxz_alpha_rz )
		
		self.data_xbr.append( self.xbr )
		self.data_ybr.append( self.ybr )
		self.data_zbr.append( self.zbr )
		self.data_Vxbr.append( self.Vxbr )
		self.data_Vybr.append( self.Vybr )
		self.data_Vzbr.append( self.Vzbr )

	def StoreSpatialCoords (self):
		self.data_x.append( self.x )
		self.data_y.append( self.y )
		self.data_z.append( self.z )
		self.data_R.append( self.R )
		self.data_Vx.append( self.Vx )
		self.data_Vy.append( self.Vy )
		self.data_Vz.append( self.Vz )
		self.data_Vel.append( self.Vel )
		#distance to the plane of satellites Conn et al. 2013
		self.distPlane = 0.158*self.x+0.769*self.y+0.620*self.z
		
		self.data_distPlane.append( self.distPlane )
		self.xbr, self.ybr, self.zbr = rotateZXZ( self.x,self.y,self.z, -zxz_gamma_rz, - zxz_beta_rx, -zxz_alpha_rz )
		self.Vxbr, self.Vybr, self.Vzbr = rotateZXZ(self.Vx,self.Vy,self.Vz, -zxz_gamma_rz, - zxz_beta_rx, -zxz_alpha_rz )
		
		self.data_xbr.append( self.xbr )
		self.data_ybr.append( self.ybr )
		self.data_zbr.append( self.zbr )
		self.data_Vxbr.append( self.Vxbr )
		self.data_Vybr.append( self.Vybr )
		self.data_Vzbr.append( self.Vzbr )

	def compute_potential( self , x_MW, y_MW, z_MW):

		#M31 parameters (from Geehan et al 2006 - via Magda)
		e=2.71828183
		rbulge=0.61
		rdisk=5.4
		rhalo=13.5
		Mbulge=2.86144
		disksurface=0.046
		halodensity=120794.0
		criticaldensity=277.72e-10
		sq_dist = (self.x*self.x+self.y*self.y+self.z*self.z)
		rsqt = np.sqrt( sq_dist )
		devidewith=(rbulge+rsqt)*(rbulge+rsqt)*rsqt

		BulgeConst=-Mbulge/devidewith
		#DiskConst=2.0*np.pi*disksurface*rdisk*rdisk*(np.power(e,-rsqt/rdisk)*(1./rdisk+1./rsqt)-1./rsqt)/(rsqt*rsqt)
		DiskConst=2.0*np.pi*disksurface*rdisk*rdisk*(np.power(e,-rsqt/rdisk)*(1.- np.power(e,-rsqt/rdisk)+rsqt/rdisk))/(rsqt*rsqt*rsqt)
		HaloConst=4.0*np.pi*halodensity*criticaldensity*(rhalo*rhalo*rhalo)*(1./(rsqt*rsqt))*(1./(rsqt+rhalo)-(np.log(rsqt/rhalo+1.)/rsqt))

		if (MW_POT =='MW'):
		
			Mhalo_MW  =91.36    #    =119.257454
			rhalo_MW  =24.54    #    =28.759113
			Mdisk_MW  =10.0     #    =5.5
			rdisk_MW  =6.65     #    =3.5
			bdisk_MW  =0.26     #    =rdisk_MW/5.0
			Mbulge_MW =3.4      #    =1.0
			rbulge_MW =0.7      #    =0.7
			#df_MW     =0.428
			#/*--------HALO_MW--------------------*/
			sq_dist_MW = ((self.x-x_MW)*(self.x-x_MW)+(self.y-y_MW)*(self.y-y_MW)+(self.z-z_MW)*(self.z-z_MW))
			rsqt_MW = np.sqrt( sq_dist_MW )
			HaloConst_MW = Mhalo_MW*(1./(sq_dist_MW*(rhalo_MW+rsqt_MW)))-(Mhalo_MW/(rsqt_MW*rsqt_MW*rsqt_MW))*(np.log((rsqt_MW/rhalo_MW)+1.) )
			#*---------DISK_MW--------------------*/
			Dxy_MW = rdisk_MW+np.sqrt((self.z-z_MW)*(self.z-z_MW)+bdisk_MW*bdisk_MW)
			div1_MW = (self.x-x_MW)*(self.x-x_MW)+(self.y-y_MW)*(self.y-y_MW)+Dxy_MW*Dxy_MW
			div1_MW = np.power(div1_MW,1.5)
			DiskConstXY_MW = -Mdisk_MW/div1_MW
			Dz_MW = np.sqrt((self.z-z_MW)*(self.z-z_MW)+bdisk_MW*bdisk_MW)
			DiskConstZ_MW = DiskConstXY_MW*Dxy_MW/Dz_MW
			#/*--------BULGE_MW-------------------*/
			divBulge_MW = rsqt_MW*(rsqt_MW+rbulge_MW)*(rsqt_MW+rbulge_MW)
			BulgeConst_MW = -Mbulge_MW/divBulge_MW

		if (MW_POT == 'NoMW'):
			HaloConst_MW   = 0.
			DiskConstXY_MW = 0.
			DiskConstZ_MW  = 0.
			BulgeConst_MW  = 0.
			x_MW           = 0.
			y_MW           = 0.  
			z_MW           = 0.  
		
		FixPotx = G*self.m*(HaloConst+DiskConst+BulgeConst)*self.x + G*self.m*(HaloConst_MW+DiskConstXY_MW+BulgeConst_MW)*(self.x-x_MW)
		FixPoty = G*self.m*(HaloConst+DiskConst+BulgeConst)*self.y + G*self.m*(HaloConst_MW+DiskConstXY_MW+BulgeConst_MW)*(self.y-y_MW)
		FixPotz = G*self.m*(HaloConst+DiskConst+BulgeConst)*self.z + G*self.m*(HaloConst_MW+DiskConstZ_MW+BulgeConst_MW)*(self.z-z_MW)
		Fx = FixPotx
		Fy = FixPoty
		Fz = FixPotz
		#print self.name,": Fx Fy Fz pot", Fx/self.m, Fy/self.m, Fz/self.m 
		return Fx, Fy, Fz

	def compute_force( self, satellite_list , x_MW, y_MW, z_MW):
		self.Fx, self.Fy, self.Fz = self.compute_potential(x_MW, y_MW, z_MW)
		if (SatSatInter == 'SSI'):
			for sat in satellite_list:
				if sat != self:
					dx = self.x- sat.x
					dy = self.y- sat.y
					dz = self.z- sat.z
					dist_sat_sat = dx*dx+dy*dy+dz*dz
					aux = - G*self.m*sat.m*np.power(dist_sat_sat,-1.5)
					self.Fx += aux*dx
					self.Fy += aux*dy
					self.Fz += aux*dz

	def LeapFrog_start( self, dt=DeltaT ):	# dt optional parameter.
		self.Vx-=(self.Fx*self.oneoverm)*0.5*dt 
		self.Vy-=(self.Fy*self.oneoverm)*0.5*dt 
		self.Vz-=(self.Fz*self.oneoverm)*0.5*dt 

	def LeapFrog_move( self, dt=DeltaT ):
		self.Vx+=self.Fx*self.oneoverm*dt
		self.Vy+=self.Fy*self.oneoverm*dt
		self.Vz+=self.Fz*self.oneoverm*dt
		self.x+=self.Vx*dt
		self.y+=self.Vy*dt
		self.z+=self.Vz*dt
		R= np.sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
		self.R = R
		Vel= np.sqrt(self.Vx*self.Vx+self.Vy*self.Vy+self.Vz*self.Vz)
		self.Vel = Vel
		self.StoreSpatialCoords()

	def LeapFrog_startBackwards( self, dt=DeltaT ):	# dt optional parameter.
		self.Vx+=(self.Fx*self.oneoverm)*0.5*dt 
		self.Vy+=(self.Fy*self.oneoverm)*0.5*dt 
		self.Vz+=(self.Fz*self.oneoverm)*0.5*dt 

	def LeapFrog_moveBackwards( self, dt=DeltaT ):
		self.Vx-=self.Fx*self.oneoverm*dt
		self.Vy-=self.Fy*self.oneoverm*dt
		self.Vz-=self.Fz*self.oneoverm*dt
		self.x-=self.Vx*dt
		self.y-=self.Vy*dt
		self.z-=self.Vz*dt
		R= np.sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
		self.R = R
		Vel= np.sqrt(self.Vx*self.Vx+self.Vy*self.Vy+self.Vz*self.Vz)
		self.Vel = Vel
		self.StoreSpatialCoords()

	def print_data( self, data_realTime):
		#print satellite data to a file 	
		fname = 'OrbitsSatellites_'+NumIntegr+'_'+SatSatInter+MW_POT+'_'+self.name+'.dat'
		f = open( fname, 'w' )
		for i in xrange(len(self.data_x)):
			print >> f, '{0:.4f}'.format(data_realTime[i]).rjust(10),'{0:.4f}'.format(self.data_x[i]).rjust(10), '{0:.4f}'.format(self.data_y[i]).rjust(10), '{0:.4f}'.format(self.data_z[i]).rjust(10) ,'{0:.4f}'.format(self.data_Vx[i]).rjust(10), '{0:.4f}'.format(self.data_Vy[i]).rjust(10), '{0:.4f}'.format(self.data_Vz[i]).rjust(10)
		f.close()

def getM31centeredCoorndinates(decl, rascen, Dist, Vls):	
	#calculate x, y and z M31-centered cartesian coordinates from right ascension, declination and distance of the satellite. 
	#Nick Bate's Code (adapted) to get eta and xi 
	tiny   =1.e-6
	#/* Trig functions */
	sdecz  = np.sin(DeclM31)
	sdec   = np.sin(decl)
	cdecz  = np.cos(DeclM31)
	cdec   = np.cos(decl)
	radif  = rascen-RascenM31
	sradif = np.sin(radif)
	cradif = np.cos(radif)
	
	#/* reciprocal of star vector length to tangent plane */
	denom = sdec*sdecz + cdec*cdecz*cradif
	
	#/* handle vectors too far from axis */
	if(denom > tiny):
		j = 0
	elif(denom >= 0):
		j = 1
		denom = tiny
	elif(denom > -1.*tiny): 
		j = 2
		denom = -1. * tiny
	else:
		j = 3
	if j > 0:
		print "j = ",j, " (0 if okay)"
	
	#/* compute tangent plane coords (even in dubious cases) */
	xi =cdec*sradif/denom
	eta = (sdec*cdecz - cdec*sdecz*cradif)/denom
	#print "xi = ", xi*180./np.pi, "eta = ", eta*180./np.pi
	
	#costheta is the cos of the angular separation between the satelite and M31		
	costheta = np.sin(decl)*np.sin(DeclM31) + np.cos(decl)*np.cos(DeclM31)*np.cos(rascen-RascenM31)	
	theta = np.arccos(costheta)
		
	#equations from Conn et al 2013 to get x, y and z and the velocities in a cartesian M31-centered coordinate system.
	#Positions
	old_Xsat_M31 = Dist*costheta*np.tan(xi)
	old_Ysat_M31 = Dist*np.sin(eta)
	old_Zsat_M31 = Dist*costheta - DM31_obs
		
	#Velocities
	#Case 1 (to get the same velocities as Ibata el al 2013): substract M31 velocity from Vls, before x, y and z components are calculated for the M31-centered coordinate system 
	#Vls = Vls - VM31_z_obs
	#signVz= 1*np.sign(Vls)	
	#old_VZsat_M31 = -signVz*(- np.abs(Vls) * costheta )     #- VM31_z_obs          
	
	#Case 2 (What I think is the correct way of doing it): substract M31 velocity only to the z-component of the line of sight velocity in the M31-centered coordinate system 
	#old_VZsat_M31 = - np.abs(Vls) * costheta - VM31_z_obs          
	
	old_VZsat_M31 = Vls*costheta
	old_VYsat_M31 = Vls*np.sin(eta)
	old_VXsat_M31 = Vls*costheta*np.tan(xi)

	return old_Xsat_M31, old_Ysat_M31, old_Zsat_M31, old_VXsat_M31, old_VYsat_M31 , old_VZsat_M31 , xi , eta 


def rotateZXZ( x0, y0, z0, alpha_rz, beta_rx, gamma_rz ): 
	cos_alpha=np.cos(alpha_rz)
	cos_beta =np.cos(beta_rx) 
	cos_gamma=np.cos(gamma_rz)
	sin_alpha=np.sin(alpha_rz)
	sin_beta =np.sin(beta_rx) 
	sin_gamma=np.sin(gamma_rz)
	rot = np.zeros([3,3])
	rot[0,0] = cos_alpha*cos_gamma - sin_alpha*cos_beta*sin_gamma
	rot[0,1] = -cos_alpha*sin_gamma - sin_alpha*cos_beta*cos_gamma
	rot[0,2] = sin_alpha*sin_beta
	rot[1,0] = sin_alpha*cos_gamma + cos_alpha*cos_beta*sin_gamma
	rot[1,1] = -sin_alpha*sin_gamma + cos_alpha*cos_beta*cos_gamma
	rot[1,2] = -cos_alpha*sin_beta
	rot[2,0] = sin_beta*sin_gamma
	rot[2,1] = sin_beta*cos_gamma
	rot[2,2] = cos_beta
	coor = np.array([x0,y0,z0])
	coor_new = np.dot(coor,rot)
	alpha_rot = np.arctan2(rot[0,2] , -rot[1,2])
	beta_rot  = np.arccos(rot[2,2])
	gamma_rot = np.arctan2(rot[2,0] , rot[2,1])
	return coor_new[0], coor_new[1], coor_new[2]


def ComputeDispersion (SatList, Dist_aver_All, Dist_aver_PlaneSat, Dist_aver_NotPlaneSat):
	DispPlane_All=0.
	SumDistPlane_All=0.
	Dist_aver_plusDisp_All=0.
	Dist_aver_minusDisp_All=0.
	DispPlane_PlaneSat=0.
	SumDistPlane_PlaneSat=0.
	Dist_aver_plusDisp_PlaneSat=0.
	Dist_aver_minusDisp_PlaneSat=0.
	DispPlane_NotPlaneSat=0.
	SumDistPlane_NotPlaneSat=0.
	Dist_aver_plusDisp_NotPlaneSat=0.
	Dist_aver_minusDisp_NotPlaneSat=0.

	for sat in SatList:
		SumDistPlane_All=SumDistPlane_All + (np.abs(sat.distPlane)-Dist_aver_All)*(np.abs(sat.distPlane)-Dist_aver_All)
		if sat.PlaneYN == 'Yes':
			SumDistPlane_PlaneSat=SumDistPlane_PlaneSat + (np.abs(sat.distPlane)-Dist_aver_PlaneSat)*(np.abs(sat.distPlane)-Dist_aver_PlaneSat)
		if sat.PlaneYN == 'No':
			SumDistPlane_NotPlaneSat=SumDistPlane_NotPlaneSat + (np.abs(sat.distPlane)-Dist_aver_NotPlaneSat)*(np.abs(sat.distPlane)-Dist_aver_NotPlaneSat)
	DispPlane_All=np.sqrt(SumDistPlane_All)/N_Sat_All
	DispPlane_PlaneSat=np.sqrt(SumDistPlane_PlaneSat)/N_Sat_PlaneSat
	DispPlane_NotPlaneSat=np.sqrt(SumDistPlane_NotPlaneSat)/N_Sat_NotPlaneSat

	Dist_aver_plusDisp_All =Dist_aver_All + DispPlane_All
	Dist_aver_minusDisp_All=Dist_aver_All - DispPlane_All
	data_DispPlane_All.append(DispPlane_All)
	data_DistAver_All.append(Dist_aver_All)
	data_DistAver_plusDisp_All.append(Dist_aver_plusDisp_All)
	data_DistAver_minusDisp_All.append(Dist_aver_minusDisp_All)

	Dist_aver_plusDisp_PlaneSat =Dist_aver_PlaneSat + DispPlane_PlaneSat
	Dist_aver_minusDisp_PlaneSat=Dist_aver_PlaneSat - DispPlane_PlaneSat
	data_DispPlane_PlaneSat.append(DispPlane_PlaneSat)
	data_DistAver_PlaneSat.append(Dist_aver_PlaneSat)
	data_DistAver_plusDisp_PlaneSat.append(Dist_aver_plusDisp_PlaneSat)
	data_DistAver_minusDisp_PlaneSat.append(Dist_aver_minusDisp_PlaneSat)

	Dist_aver_plusDisp_NotPlaneSat =Dist_aver_NotPlaneSat + DispPlane_NotPlaneSat
	Dist_aver_minusDisp_NotPlaneSat=Dist_aver_NotPlaneSat - DispPlane_NotPlaneSat
	data_DispPlane_NotPlaneSat.append(DispPlane_NotPlaneSat)
	data_DistAver_NotPlaneSat.append(Dist_aver_NotPlaneSat)
	data_DistAver_plusDisp_NotPlaneSat.append(Dist_aver_plusDisp_NotPlaneSat)
	data_DistAver_minusDisp_NotPlaneSat.append(Dist_aver_minusDisp_NotPlaneSat)

	Dist_aver_All=0.
	Dist_aver_PlaneSat=0.
	Dist_aver_NotPlaneSat=0.


#############################################################################################################################################
#routines to calculate the position of the Milky Way in the M31 coordinate system and the changes in this positons due to the potential of M31. 

def compute_potentialM31(x_MW, y_MW, z_MW):

	#Compute M31 potential in the location of the Milky Way
	
	#M31 parameters (from Geehan et al 2006 - via Magda)
	e=2.71828183
	rbulge=0.61
	rdisk=5.4
	rhalo=13.5
	Mbulge=2.86144
	disksurface=0.046
	halodensity=120794.0
	criticaldensity=277.72e-10
	sq_dist = (x_MW*x_MW+y_MW*y_MW+z_MW*z_MW)
	rsqt = np.sqrt( sq_dist )
	devidewith=(rbulge+rsqt)*(rbulge+rsqt)*rsqt

	BulgeConst=-Mbulge/devidewith
	#DiskConst=2.0*np.pi*disksurface*rdisk*rdisk*(np.power(e,-rsqt/rdisk)*(1./rdisk+1./rsqt)-1./rsqt)/(rsqt*rsqt)
	DiskConst=2.0*np.pi*disksurface*rdisk*rdisk*(np.power(e,-rsqt/rdisk)*(1.- np.power(e,-rsqt/rdisk)+rsqt))/(rsqt*rsqt*rsqt)
	HaloConst=4.0*np.pi*halodensity*criticaldensity*(rhalo*rhalo*rhalo)*(1./(rsqt*rsqt))*(1./(rsqt+rhalo)-(np.log(rsqt/rhalo+1.)/rsqt))
	
	FixPotx = G*m_MW*(HaloConst+DiskConst+BulgeConst)*x_MW
	FixPoty = G*m_MW*(HaloConst+DiskConst+BulgeConst)*y_MW 
	FixPotz = G*m_MW*(HaloConst+DiskConst+BulgeConst)*z_MW 
	Fx = FixPotx
	Fy = FixPoty
	Fz = FixPotz
	#print self.name,": Fx Fy Fz pot", Fx/self.m, Fy/self.m, Fz/self.m 
	return Fx, Fy, Fz

def LeapFrog_startMW(x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW, dt=DeltaT ):	# dt optional parameter.
	MWFx, MWFy, MWFz = compute_potentialM31(x_MW, y_MW, z_MW) 
	Vx_MW-=(MWFx*oneoverm_MW)*0.5*dt 
	Vy_MW-=(MWFy*oneoverm_MW)*0.5*dt 
	Vz_MW-=(MWFz*oneoverm_MW)*0.5*dt 
	return Vx_MW, Vy_MW, Vz_MW

def LeapFrog_moveMW(x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW, dt=DeltaT ):
	MWFx, MWFy, MWFz = compute_potentialM31(x_MW, y_MW, z_MW) 
	Vx_MW+=(MWFx*oneoverm_MW)*0.5*dt
	Vy_MW+=(MWFy*oneoverm_MW)*0.5*dt
	Vz_MW+=(MWFz*oneoverm_MW)*0.5*dt
	x_MW+=Vx_MW*dt
	y_MW+=Vy_MW*dt
	z_MW+=Vz_MW*dt
	return x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW

def LeapFrog_startBackwardsMW(x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW, dt=DeltaT ):	# dt optional parameter.
	MWFx, MWFy, MWFz = compute_potentialM31(x_MW, y_MW, z_MW) 
	Vx_MW+=(MWFx*oneoverm_MW)*0.5*dt 
	Vy_MW+=(MWFy*oneoverm_MW)*0.5*dt 
	Vz_MW+=(MWFz*oneoverm_MW)*0.5*dt 
	return Vx_MW, Vy_MW, Vz_MW

def LeapFrog_moveBackwardsMW(x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW, dt=DeltaT ):
	MWFx, MWFy, MWFz = compute_potentialM31(x_MW, y_MW, z_MW) 
	Vx_MW-=(MWFx*oneoverm_MW)*0.5*dt
	Vy_MW-=(MWFy*oneoverm_MW)*0.5*dt
	Vz_MW-=(MWFz*oneoverm_MW)*0.5*dt
	x_MW-=Vx_MW*dt
	y_MW-=Vy_MW*dt
	z_MW-=Vz_MW*dt
	return x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW

#end of routines to calculate the position of the Milky Way in the M3 coordinate system and the changes in this positons due to the potential of M31. 
#####################################################################################################################################################

##############################PLOTING ROUTINES#######################################################################################################################
def plotSnapshotsAll_beforeRotation_VelMin300 (SatList, tSnap, dM31obs, vM31x_obs, vM31y_obs, vM31z_obs):

	#calculate initial MW position in the M31-centric coordinate system
	old_x_MW = 0.
	old_y_MW = 0.	
	old_z_MW = - dM31obs

	plt.figure()
	plt.xlim( -500., 500. )
	plt.ylim( -500., 500. )
	ax = plt.gca()
	circle1 = Circle( (old_x_MW,old_y_MW), 3.5, facecolor='none', edgecolor='k', linewidth=2 )
	circle2 = Circle( (0.,0.), 3.5, facecolor='none', edgecolor='k', linewidth=2 )
	ax.add_artist(circle1)
	ax.add_artist(circle2)
	for sat in Satellite:
		Vls_xbr, Vls_ybr, Vls_zbr = rotateZXZ(sat.Vls_x, sat.Vls_y, sat.Vls_z, -zxz_gamma_rz, -zxz_beta_rx, -zxz_alpha_rz)
		Vtan_xbr, Vtan_ybr, Vtan_zbr = rotateZXZ(sat.Vtan_x, sat.Vtan_y, sat.Vtan_z, -zxz_gamma_rz, -zxz_beta_rx, -zxz_alpha_rz)
		if sat.PlaneYN == "Yes":
	        	ax = plt.gca()
			plt.plot(sat.data_xbr, sat.data_ybr, color="grey", linewidth=0.5)
			if sat.name == TestSatName: 
				plt.plot(sat.data_xbr, sat.data_ybr, color="black", linewidth=1.5)
			plt.plot(sat.data_xbr[0], sat.data_ybr[0], '.', markersize=10, color= "red")
			ax.quiver(sat.data_xbr[0], sat.data_ybr[0], Vls_xbr, Vls_ybr ,angles='xy',scale_units='xy',scale=1, color= 'red' , linewidth=0.15)
		if sat.PlaneYN == "YesNO":
	        	ax = plt.gca()
			#plt.plot(sat.data_xbr, sat.data_ybr, color="grey", linewidth=0.5)
			plt.plot(sat.data_xbr[0], sat.data_ybr[0], '*', markersize=10, color= "red")
			ax.quiver(sat.data_xbr[0], sat.data_ybr[0], Vls_xbr, Vls_ybr ,angles='xy',scale_units='xy',scale=1, color= 'red' , linewidth=0.15)
		if sat.PlaneYN == "No":
			plt.plot(sat.data_xbr[0], sat.data_ybr[0], '.', markersize=10, color= "blue")
	plt.xlabel('X (kpc)', fontsize=15) 
	plt.ylabel('Y (kpc)', fontsize=15)
	plt.savefig('SnapshotSatellites_All_BeforeRotationM300_xy_'+TestSatName+'_'+NumIntegr+'_'+SatSatInter+MW_POT+'_'+tSnap+'.pdf')

	plt.figure()
	plt.xlim( -500., 500. )
	plt.ylim( -500., 500. )
	ax = plt.gca()
	circle1 = Circle( (old_x_MW,old_z_MW), 3.5, facecolor='none', edgecolor='k', linewidth=2 )
	circle2 = Circle( (0.,0.), 3.5, facecolor='none', edgecolor='k', linewidth=2 )
	ax.add_artist(circle1)
	ax.add_artist(circle2)
	for sat in Satellite:
		Vls_xbr, Vls_ybr, Vls_zbr = rotateZXZ(sat.Vls_x, sat.Vls_y, sat.Vls_z, -zxz_gamma_rz, -zxz_beta_rx, -zxz_alpha_rz)
		Vtan_xbr, Vtan_ybr, Vtan_zbr = rotateZXZ(sat.Vtan_x, sat.Vtan_y, sat.Vtan_z, -zxz_gamma_rz, -zxz_beta_rx, -zxz_alpha_rz)
		if sat.PlaneYN == "Yes":
	        	ax = plt.gca()
			plt.plot(sat.data_xbr, sat.data_zbr, color="grey", linewidth=0.5)
			if sat.name == TestSatName: 
				plt.plot(sat.data_xbr, sat.data_zbr, color="black", linewidth=1.5)
			plt.plot(sat.data_xbr[0], sat.data_zbr[0], '.', markersize=10, color= "red")
			ax.quiver(sat.data_xbr[0], sat.data_zbr[0], Vls_xbr, Vls_zbr-vM31z_obs ,angles='xy',scale_units='xy',scale=1, color= 'red' , linewidth=0.15)
		if sat.PlaneYN == "YesNO":
	        	ax = plt.gca()
			#plt.plot(sat.data_xbr, sat.data_zbr, color="grey", linewidth=0.5)
			plt.plot(sat.data_xbr[0], sat.data_zbr[0], '*', markersize=10, color= "red")
			ax.quiver(sat.data_xbr[0], sat.data_zbr[0], Vls_xbr, Vls_zbr-vM31z_obs ,angles='xy',scale_units='xy',scale=1, color= 'red' , linewidth=0.15)
		if sat.PlaneYN == "No":
			plt.plot(sat.data_xbr[0], sat.data_zbr[0], '.', markersize=10, color= "blue")
	plt.xlabel('X (kpc)', fontsize=15) 
	plt.ylabel('Z (kpc)', fontsize=15)
	plt.savefig('SnapshotSatellites_All_BeforeRotationM300_xz_'+TestSatName+'_'+NumIntegr+'_'+SatSatInter+MW_POT+'_'+tSnap+'.pdf')

	plt.figure()
	plt.xlim( -500., 500. )
	plt.ylim( -500., 500. )
	ax = plt.gca()
	circle1 = Circle( (old_y_MW,old_z_MW), 3.5, facecolor='none', edgecolor='k', linewidth=2 )
	circle2 = Circle( (0.,0.), 3.5, facecolor='none', edgecolor='k', linewidth=2 )
	ax.add_artist(circle1)
	ax.add_artist(circle2)
	for sat in Satellite:
		Vls_xbr, Vls_ybr, Vls_zbr = rotateZXZ(sat.Vls_x, sat.Vls_y, sat.Vls_z, -zxz_gamma_rz, -zxz_beta_rx, -zxz_alpha_rz)
		Vtan_xbr, Vtan_ybr, Vtan_zbr = rotateZXZ(sat.Vtan_x, sat.Vtan_y, sat.Vtan_z, -zxz_gamma_rz, -zxz_beta_rx, -zxz_alpha_rz)
		if sat.PlaneYN == "Yes":
	        	ax = plt.gca()
			plt.plot(sat.data_ybr, sat.data_zbr, color="grey", linewidth=0.5)
			if sat.name == TestSatName: 
				plt.plot(sat.data_ybr, sat.data_zbr, color="black", linewidth=1.5)
			plt.plot(sat.data_ybr[0], sat.data_zbr[0], '.', markersize=10, color= "red")
			ax.quiver(sat.data_ybr[0], sat.data_zbr[0], Vls_ybr, Vls_zbr-vM31z_obs ,angles='xy',scale_units='xy',scale=1, color= 'red' , linewidth=0.15)
		if sat.PlaneYN == "YesNO":
	        	ax = plt.gca()
			#plt.plot(sat.data_ybr, sat.data_zbr, color="grey", linewidth=0.5)
			plt.plot(sat.data_ybr[0], sat.data_zbr[0], '*', markersize=10, color= "red")
			ax.quiver(sat.data_ybr[0], sat.data_zbr[0], Vls_ybr, Vls_zbr-vM31z_obs ,angles='xy',scale_units='xy',scale=1, color= 'red' , linewidth=0.15)
		if sat.PlaneYN == "No":
			plt.plot(sat.data_ybr[0], sat.data_zbr[0], '.', markersize=10, color= "blue")
	plt.xlabel('Y (kpc)', fontsize=15) 
	plt.ylabel('Z (kpc)', fontsize=15)
	plt.savefig('SnapshotSatellites_All_BeforeRotationM300_yz_'+TestSatName+'_'+NumIntegr+'_'+SatSatInter+MW_POT+'_'+tSnap+'.pdf')




#****************************************************************************
#************************  Main Program  ************************************
#****************************************************************************


#Initialization of the Satellites: Andromeda (Distances Paper Conn et al. 2012, velocities Collins et al 2013 and Tollerud et al 2012) 
Satellite = []
#And I
#Satellite.append( body( Dist=727, errDistplus= +18 , errDistminus= -17 ,rascenHH=00. , rascenMM=45. , rascenSS=39.8 , declG=38. , declMM=02. , declSS=28.0 , Vls=-376.3, errVlsplus= +2.2 , errVlsminus= -2.2, m0=7.0  , name = 'AndI'          , PlaneYN = 'Yes', Vtany_param = -265.0) )
###And II ********************************************************************
#Satellite.append( body( Dist=630, errDistplus= +15 , errDistminus= -15 ,rascenHH=  . , rascenMM=  . , rascenSS=     , declG=  . , declMM=  . , declSS=     , Vls=      , errVlsplus=      , errVlsminus=     , m0=6.1  , name = 'AndII'         , PlaneYN = 'No', Vtany_param = ) ) 
#And III 
Satellite.append( body( Dist=723, errDistplus= +18 , errDistminus= -24 ,rascenHH=00. , rascenMM=35. , rascenSS=33.8 , declG=36. , declMM=29. , declSS=52.0 , Vls=-344.3, errVlsplus= +1.7 , errVlsminus= -1.7, m0=0.96 , name = 'AndIII'        , PlaneYN = 'Yes', Vtany_param = -330.0) ) 
#And V    
Satellite.append( body( Dist=742, errDistplus= +21 , errDistminus= -22 ,rascenHH=01. , rascenMM=10. , rascenSS=17.1 , declG=47. , declMM=34. , declSS=41.0 , Vls=-391.5, errVlsplus= +2.7 , errVlsminus= -2.7, m0=2.6  , name = 'AndV'          , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And VI    
#Satellite.append( body( Dist=783, errDistplus= +28 , errDistminus= -28 ,rascenHH=23. , rascenMM=51. , rascenSS=39.0 , declG=24. , declMM=35. , declSS=42.0 , Vls=-339.8, errVlsplus= +1.8 , errVlsminus= -1.8, m0=4.7  , name = 'AndVI'         , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And VII    
#Satellite.append( body( Dist=762, errDistplus= +25 , errDistminus= -24 ,rascenHH=23. , rascenMM=26. , rascenSS=31.7 , declG=50. , declMM=40. , declSS=33.0 , Vls=-307.2, errVlsplus= +1.3 , errVlsminus= -1.3, m0=6.9  , name = 'AndVII'        , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And IX # este se va lejos solo  ********************************************* cuidado masa inventada   
#Satellite.append( body( Dist=600, errDistplus= +91 , errDistminus= -23 ,rascenHH=00. , rascenMM=52. , rascenSS=53.0 , declG=43. , declMM=11. , declSS=45.0 , Vls=-209.4, errVlsplus= +2.5 , errVlsminus= -2.5, m0= 1.  , name = 'AndIX'         , PlaneYN = 'YesNO', Vtany_param = -130.) ) 
#And X                                                           
#Satellite.append( body( Dist=670, errDistplus= +24 , errDistminus= -39 ,rascenHH=01. , rascenMM=06. , rascenSS=33.7 , declG=44. , declMM=48. , declSS=16.0 , Vls=-164.1, errVlsplus= +1.7 , errVlsminus= -1.7, m0=0.47 , name = 'AndX'          , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XI                                                                                       
#Satellite.append( body( Dist=763, errDistplus= +29 , errDistminus= -106,rascenHH=00. , rascenMM=46. , rascenSS=20.0 , declG=33. , declMM=48. , declSS=05.0 , Vls=-427.5, errVlsplus= +3.5 , errVlsminus= -3.4, m0=0.53 , name = 'AndXI'         , PlaneYN = 'Yes', Vtany_param = -190.0) ) 
#And XII *********************************************** cuidado masa cambiada (de 0.0 a 0.01) 
#Satellite.append( body( Dist=928, errDistplus= +40 , errDistminus= -136,rascenHH=00. , rascenMM=47. , rascenSS=27.0 , declG=34. , declMM=22. , declSS=29.0 , Vls=-557.1, errVlsplus= +1.7 , errVlsminus= -1.7, m0=0.01 , name = 'AndXII'        , PlaneYN = 'YesNO', Vtany_param = 0.) ) 
#And XIII  *********************************************** cuidado masa cambiada (de 0.0 a 0.01)
#Satellite.append( body( Dist=760, errDistplus= +126, errDistminus= -154,rascenHH=00. , rascenMM=51. , rascenSS=51.0 , declG=33. , declMM=00. , declSS=16.0 , Vls=-204.8, errVlsplus= +4.9 , errVlsminus= -4.9, m0=0.01 , name = 'AndXIII'       , PlaneYN = 'YesNO', Vtany_param = 50.0) ) 
#And XIV               
#Satellite.append( body( Dist=793, errDistplus= +23 , errDistminus= -179,rascenHH=00. , rascenMM=51. , rascenSS=35.0 , declG=29. , declMM=41. , declSS=49.0 , Vls=-480.6, errVlsplus= +1.2 , errVlsminus= -1.2, m0=0.92 , name = 'AndXIV'        , PlaneYN = 'Yes', Vtany_param = -90) ) 
#And XV *********************************************** cuidado masa inventada  
#Satellite.append( body( Dist=626, errDistplus= +79 , errDistminus= -35 ,rascenHH=01. , rascenMM=14. , rascenSS=18.7 , declG=38. , declMM=07. , declSS=03.0 , Vls=-323.0, errVlsplus= +1.4 , errVlsminus= -1.4, m0= 1.  , name = 'AndXV'         , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XVI   *********************************************** cuidado masa inventada   
#Satellite.append( body( Dist=476, errDistplus= +44 , errDistminus= -29 ,rascenHH=00. , rascenMM=59. , rascenSS=29.8 , declG=32. , declMM=22. , declSS=36.0 , Vls=-367.3, errVlsplus= +2.8 , errVlsminus= -2.8, m0= 1.  , name = 'AndXVI'        , PlaneYN = 'YesNO', Vtany_param = 5.) ) 
#And XVII #este se va a la miercoles solo #                
#Satellite.append( body( Dist=727, errDistplus= +39 , errDistminus= -25 ,rascenHH=00. , rascenMM=37. , rascenSS=07.0 , declG=44. , declMM=19. , declSS=20.0 , Vls=-251.6, errVlsplus= +1.8 , errVlsminus= -2.0, m0=0.13 , name = 'AndXVII'       , PlaneYN = 'Yes', Vtany_param = -240.0) ) 
#And XVIII  *********************************************** cuidado masa cambiada (de 0.0 a 0.01) #puede no estar asociada con M31 
#Satellite.append( body( Dist=1214, errDistplus= +40 , errDistminus= -43 ,rascenHH=00. , rascenMM=02. , rascenSS=14.5 , declG=45. , declMM=05. , declSS=20.0 , Vls=-346.8, errVlsplus= +2.0 , errVlsminus= -2.0, m0=0.01 , name = 'AndXVIII'      , PlaneYN = 'No', Vtany_param = 0.) ) 
#And XIX                              
#Satellite.append( body( Dist=821, errDistplus= +32 , errDistminus= -148,rascenHH=00. , rascenMM=19. , rascenSS=32.1 , declG=35. , declMM=02. , declSS=37.1 , Vls=-111.6, errVlsplus= +1.6 , errVlsminus= -1.4, m0=1.9  , name = 'AndXIX'        , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XX                                
#Satellite.append( body( Dist=741, errDistplus= +42 , errDistminus= -52 ,rascenHH=00. , rascenMM=07. , rascenSS=30.7 , declG=35. , declMM=07. , declSS=56.4 , Vls=-456.2, errVlsplus= +3.1 , errVlsminus= -3.6, m0=0.33 , name = 'AndXX'         , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XXI                                
#Satellite.append( body( Dist=827, errDistplus= +23 , errDistminus= -25 ,rascenHH=23. , rascenMM=54. , rascenSS=47.7 , declG=42. , declMM=28. , declSS=15.0 , Vls=-362.5, errVlsplus= +0.9 , errVlsminus= -0.9, m0=0.99 , name = 'AndXXI'        , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XXII (Tri)                          
#Satellite.append( body( Dist=920, errDistplus= +32 , errDistminus= -139,rascenHH=01. , rascenMM=27. , rascenSS=40.0 , declG=28. , declMM=05. , declSS=25.0 , Vls=-129.8, errVlsplus= +2.0 , errVlsminus= -2.0, m0=0.11 , name = 'AndXXII-Tri'   , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XXIII                                
#Satellite.append( body( Dist=748, errDistplus= +31 , errDistminus= -21 ,rascenHH=01. , rascenMM=29. , rascenSS=21.8 , declG=38. , declMM=43. , declSS=08.0 , Vls=-237.7, errVlsplus= +1.2 , errVlsminus= -1.2, m0=2.9  , name = 'AndXXIII'      , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XXIV                                  
#Satellite.append( body( Dist=898, errDistplus= +28 , errDistminus= -42 ,rascenHH=01. , rascenMM=18. , rascenSS=30.0 , declG=46. , declMM=21. , declSS=58.0 , Vls=-128.2, errVlsplus= +5.2 , errVlsminus= -5.2, m0=0.4  , name = 'AndXXIV'       , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XXV                                    
#Satellite.append( body( Dist=736, errDistplus= +23 , errDistminus= -69 ,rascenHH=00. , rascenMM=30. , rascenSS=08.9 , declG=46. , declMM=51. , declSS=07.0 , Vls=-107.8, errVlsplus= +1.0 , errVlsminus= -1.0, m0=0.34 , name = 'AndXXV'        , PlaneYN = 'Yes', Vtany_param = -195.0) ) 
#And XXVI                                    
#Satellite.append( body( Dist=754, errDistplus= +218, errDistminus= -164,rascenHH=00. , rascenMM=23. , rascenSS=45.6 , declG=47. , declMM=54. , declSS=58.0 , Vls=-261.6, errVlsplus= +3.0 , errVlsminus= -2.8, m0=0.96 , name = 'AndXXVI'       , PlaneYN = 'YesNO', Vtany_param = -50.) ) 
##And XXVII                                 
#Satellite.append( body(Dist=1255, errDistplus= +42 , errDistminus= -474,rascenHH=00. , rascenMM=37. , rascenSS=27.2 , declG=45. , declMM=23. , declSS=13.0 , Vls=-539.6, errVlsplus= +4.7 , errVlsminus= -4.5, m0=8.3  , name = 'AndXXVII'      , PlaneYN = 'YesNO', Vtany_param = -110.) ) 
#And XXVIII                                    
#Satellite.append( body( Dist=650, errDistplus= +150, errDistminus= -80 ,rascenHH=22. , rascenMM=32. , rascenSS=41.2 , declG=31. , declMM=12. , declSS=51.2 , Vls=-326.2, errVlsplus= +2.7 , errVlsminus= -2.7, m0=0.53 , name = 'AndXXVIII'     , PlaneYN = 'No' , Vtany_param = 0.) ) 
#And XXX (CassII)                               
#Satellite.append( body( Dist=681, errDistplus= +32 , errDistminus= -78 ,rascenHH=00. , rascenMM=36. , rascenSS=34.9 , declG=49. , declMM=38. , declSS=48.0 , Vls=-139.8, errVlsplus= +6.0 , errVlsminus= -6.6, m0=2.2  , name = 'AndXXX-CassII'                , PlaneYN = 'YesNO' , Vtany_param = -70) ) 
# NGC 147                                         
#Satellite.append( body( Dist=712, errDistplus= +21 , errDistminus= -19 , rascenHH=00. , rascenMM=33. , rascenSS=12.1 , declG=48. , declMM=30. , declSS=32, Vls= -193.1 ,errVlsplus=0.8 , errVlsminus=0.8 , m0=56. ,  name = 'NGC147'                           , PlaneYN = 'Yes' , Vtany_param = -230) ) 
# NGC 185                                        
#Satellite.append( body( Dist=620, errDistplus= +19 , errDistminus= -18 ,rascenHH=00. , rascenMM= 38. , rascenSS=58.0 , declG=48. , declMM=20. , declSS=15     , Vls= -203.8      ,errVlsplus= 1.1     , errVlsminus= 1.1 , m0=72.   ,  name = 'NGC185'         , PlaneYN = 'YesNO' , Vtany_param = -98.) ) 


#initialization of some variables:
cont = 0
Dist_aver_All=0.
data_DispPlane_All = []
data_DistAver_All = []
data_DistAver_plusDisp_All = []
data_DistAver_minusDisp_All = []

Dist_aver_PlaneSat=0.
data_DispPlane_PlaneSat = []
data_DistAver_PlaneSat = []
data_DistAver_plusDisp_PlaneSat = []
data_DistAver_minusDisp_PlaneSat = []

Dist_aver_NotPlaneSat=0.
data_DispPlane_NotPlaneSat = []
data_DistAver_NotPlaneSat = []
data_DistAver_plusDisp_NotPlaneSat = []
data_DistAver_minusDisp_NotPlaneSat = []

N_Sat_All=0
N_Sat_PlaneSat=0
N_Sat_NotPlaneSat= 0

data_MWx=[]
data_MWy=[]
data_MWz=[]

data_realTime=[]

#calculate initial MW position and velocity in the new M31-centric coordinate system
old_x_MW = 0.
old_y_MW = 0.	
old_z_MW = - DM31_obs
old_Vx_MW = 0.#-VM31_x_obs #0.
old_Vy_MW = 0.#-VM31_y_obs #0.
old_Vz_MW = 109 #- VM31_z_obs

x_MW0, y_MW0, z_MW0 = rotateZXZ( old_x_MW, old_y_MW, old_z_MW, zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
Vx_MW0, Vy_MW0, Vz_MW0 = rotateZXZ( old_Vx_MW, old_Vy_MW, old_Vz_MW, zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )
x_MW, y_MW, z_MW = x_MW0, y_MW0, z_MW0

#initialization of the numerical integration
for sat in Satellite:
	sat.compute_force( Satellite, x_MW=x_MW0, y_MW=y_MW0, z_MW=z_MW0 )
	#sat.CalculateDistanceSatSat( Satellite) 
	if NumIntegr == 'Fwd' :
		sat.LeapFrog_start()
	if NumIntegr == 'Bwd' :
		sat.LeapFrog_startBackwards()
	
	#count number of satellites (all, plane sats, not-plane sats) 	
	N_Sat_All=N_Sat_All+1	
	if sat.PlaneYN == 'Yes' :
		N_Sat_PlaneSat=N_Sat_PlaneSat+1	
	if sat.PlaneYN == 'No' :
		N_Sat_NotPlaneSat=N_Sat_NotPlaneSat+1	
	
#initialization of the numerical integration to get the Milky Way's location in The M31 coordinate system
if NumIntegr == 'Fwd' :
	Vx_MW, Vy_MW, Vz_MW = LeapFrog_startMW(x_MW=x_MW0, y_MW=y_MW0, z_MW=z_MW0, Vx_MW=Vx_MW0, Vy_MW=Vy_MW0, Vz_MW=Vz_MW0, dt=DeltaT)
if NumIntegr == 'Bwd' :
	Vx_MW, Vy_MW, Vz_MW = LeapFrog_startBackwardsMW(x_MW=x_MW0, y_MW=y_MW0, z_MW=z_MW0, Vx_MW=Vx_MW0, Vy_MW=Vy_MW0, Vz_MW=Vz_MW0, dt=DeltaT)


#calculate distance average
for sat in Satellite:
	Dist_aver_All=Dist_aver_All + np.sqrt(sat.distPlane*sat.distPlane)
	if sat.PlaneYN == 'Yes':
		Dist_aver_PlaneSat=Dist_aver_PlaneSat + np.sqrt(sat.distPlane*sat.distPlane)
	if sat.PlaneYN == 'No':
		Dist_aver_NotPlaneSat=Dist_aver_NotPlaneSat + np.sqrt(sat.distPlane*sat.distPlane)
Dist_aver_All=Dist_aver_All/N_Sat_All
Dist_aver_PlaneSat=Dist_aver_PlaneSat/N_Sat_PlaneSat
Dist_aver_NotPlaneSat=Dist_aver_NotPlaneSat/N_Sat_NotPlaneSat
#call routine to compute dispersion around distance average
ComputeDispersion (Satellite, Dist_aver_All, Dist_aver_PlaneSat, Dist_aver_NotPlaneSat)

realTime = 0.
data_realTime.append(0.0000001*realTime)
tSnap = "0000"


#numerical integration
for t in xrange(int(T/DeltaT)):

	# Case 1: ignoring M31's potential
	#old_x_MW = 0.
	#old_y_MW = 0.
	#if NumIntegr == 'Fwd' :
	#	old_z_MW = old_z_MW -VM31*DeltaT
	#if NumIntegr == 'Bwd' :
	#	old_z_MW = old_z_MW +VM31*DeltaT
	#data_oldMWx.append(old_x_MW)
	#data_oldMWy.append(old_y_MW)
	#data_oldMWz.append(old_z_MW)
	#x_MW, y_MW, z_MW = rotateZXZ( old_x_MW, old_y_MW, old_z_MW, zxz_alpha_rz, zxz_beta_rx, zxz_gamma_rz )

	# Case2: considering M31's potential
	if NumIntegr == 'Fwd' :
		x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW = LeapFrog_moveMW(x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW, dt=DeltaT)
	if NumIntegr == 'Bwd' :
		x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW = LeapFrog_moveBackwardsMW(x_MW, y_MW, z_MW, Vx_MW, Vy_MW, Vz_MW, dt=DeltaT)

	data_MWx.append(x_MW)
	data_MWy.append(y_MW)
	data_MWz.append(z_MW)


	#orbit integration
	for sat in Satellite:
		sat.compute_force( Satellite, x_MW, y_MW, z_MW )
	for sat in Satellite:
		if NumIntegr == 'Fwd' :
			sat.LeapFrog_move()
		if NumIntegr == 'Bwd' :
			sat.LeapFrog_moveBackwards()

	#calculate distance average
		Dist_aver_All=Dist_aver_All + np.sqrt(sat.distPlane*sat.distPlane)
		if sat.PlaneYN == 'Yes':
			Dist_aver_PlaneSat=Dist_aver_PlaneSat + np.sqrt(sat.distPlane*sat.distPlane)
		if sat.PlaneYN == 'No':
			Dist_aver_NotPlaneSat=Dist_aver_NotPlaneSat + np.sqrt(sat.distPlane*sat.distPlane)
	Dist_aver_All=Dist_aver_All/N_Sat_All
	Dist_aver_PlaneSat=Dist_aver_PlaneSat/N_Sat_PlaneSat
	Dist_aver_NotPlaneSat=Dist_aver_NotPlaneSat/N_Sat_NotPlaneSat
	
	#call routine to compute dispersion around distance average
	ComputeDispersion (Satellite, Dist_aver_All, Dist_aver_PlaneSat, Dist_aver_NotPlaneSat)

	if NumIntegr == 'Fwd' :
		realTime = realTime + DeltaT
		data_realTime.append(realTime)
	if NumIntegr == 'Bwd' :
		realTime = realTime - DeltaT
		data_realTime.append(realTime)


plotSnapshotsAll_beforeRotation_VelMin300 (Satellite, tSnap, DM31_obs, VM31_x_obs, VM31_y_obs, VM31_z_obs)

















