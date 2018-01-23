# Functions extracting emulator and any other data from LES output NetCDF files,
# and collection of functions for generating LES inputs.
#
#	Tomi Raatikanen 8.1.2018
#
# Functions
# =========
# Use Python import to make these functions available, e.g. from LES2emu import GetEmuVars, get_netcdf_variable
# a) Functions for extracting data from the LES outputs
#	GetEmuVars(fname,tstart,tend)
#	get_netcdf_variable(fname,var_name,target_time,[end_time])
#	extract_write_data(fname_template,specs,[name_out,nmax])
#
# b) Functions for generating LES inputs
#	calc_psat_w(T)
#	calc_sat_mixr(p,T)
#	calc_rh(rw,T,press)
#	calc_cloud_base(p_surf,theta,rw)
#	calc_lwc_altitude(p_surf,theta,rw,zz)
#	solve_rw(p_surf,theta,lwc,zz)
#
# Notes
# =====
# 1) Input file name should contain complete path in addition to the file name
# 	e.g. '/ibrix/arch/ClimRes/aholaj/case_emulator_DESIGN_v1.4.0_LES_cray.dev20170324_LVL4/emul01/emul01.ts.nc'
# 2) Voima requires python 2.7.10, so execute "module load Python/2.7.10"
#


def GetEmuVars(fname,tstart,tend,ttol=3600.,start_offset=0,end_offset=0):
	# Function calculates LES output variables for the emulator as defined in the ECLAIR proof-of-concept document
	#	https://docs.google.com/document/d/1L-YyJLhtmLYg4rJYo5biOW96eeRC7z_trZsow_8TbeE/edit
	# Inputs:
	# 	fname			Complete path and name of a time statistics file (*.ts.nc)
	#	start, tend	Time window (s)
	# Optional inputs
	#	ttol			Time tolelance (s) for finding averaging window
	#	start_offset	Point offset to start index (time in the NetCDF files is the save time, so typically should ignore the first points)
	#	end_offset	-||- end ind index
	#
	# Example:
	#	file='/ibrix/arch/ClimRes/aholaj/case_emulator_DESIGN_v1.4.0_LES_cray.dev20170324_LVL4/emul01/emul01.ts.nc'
	#	tstart=2.5*3600
	#	tend=3.5*3600
	#	cfrac, CDNC, prcp, dn, we, cfrac_std, CDNC_std, prcp_std, dn_std, we_std = GetEmuVars(file,tstart,ttol=10.,start_offset=1)
	#
	import os
	import netCDF4 as netcdf
	import numpy
	#
	# Outputs
	cfrac=-999.	# Cloud fraction
	CDNC=-999.	# Cloud droplet number concentration in cloudy columns (#/kg)
	prcp=-999.		# Precipitation tendency = domain mean surface precipitation (kg/m^2/s)
	dn=-999.		# In-cloud aerosol number loss tendency = change in interstitial aerosol+cloud droplet number concentration in cloudy columns (#/kg/s)
	we=-999.		# Mixing between FT and BL = domeain mean entrainment velocity (m/s)
	#
	# ... and their standard deviations
	cfrac_std=-999.; CDNC_std=-999.; prcp_std=-999.; dn_std=-999.; we_std=-999.
	#
	# File must exist
	if not os.path.lexists(fname):
		print fname+' not found!'
		return cfrac, CDNC, prcp, dn, we, cfrac_std, CDNC_std, prcp_std, dn_std, we_std,
	#
	# Open the target NetCDF file
	ncid = netcdf.Dataset(fname,'r')
	#
	# Time
	# ====
	if 'time' not in ncid.variables:
		print 'Time not found from '+fname+'!'
		return cfrac, CDNC, prcp, dn, we, cfrac_std, CDNC_std, prcp_std, dn_std, we_std,
	times = ncid.variables['time']
	#
	# Find the closest matching time points
	ind_tstart=0
	ind_tend=0
	i=0
	for t in times:
		if abs(t-tstart)<abs(times[ind_tstart]-tstart): ind_tstart=i
		if abs(t-tend)<abs(times[ind_tend]-tend): ind_tend=i
		i+=1
	#
	if abs(times[ind_tstart]-tstart)>ttol or abs(times[ind_tend]-tend)>ttol:
		print 'Matching start or end time not found from '+fname+'!'
		return cfrac, CDNC, prcp, dn, we, cfrac_std, CDNC_std, prcp_std, dn_std, we_std,
	#
	# Apply offset (typically the first point is ignored
	ind_tstart+=start_offset
	ind_tend+=end_offset
	if ind_tstart<0 or ind_tstart>ind_tend or ind_tend>=len(times):
		print 'Invalid data range for '+fname+': ',ind_tstart,ind_tend,len(times)
		return cfrac, CDNC, prcp, dn, we, cfrac_std, CDNC_std, prcp_std, dn_std, we_std,
	#
	#
	# Outputs
	# ========
	# Cloud fraction
	if 'cfrac' not in ncid.variables:
		print 'Cloud fraction not found from '+fname+'!'
		return cfrac, CDNC, prcp, dn, we, cfrac_std, CDNC_std, prcp_std, dn_std, we_std,
	#
	# Need cloud fractions for normalizing domain mean interstitial and cloud droplet number concentrations
	cfrac_ts=ncid.variables['cfrac'][ind_tstart:ind_tend]
	ncfrac = sum( cfrac_ts>0.0 ) # The number of non-zero cloud fractions
	#
	cfrac = numpy.mean( cfrac_ts )
	cfrac_std = numpy.std( cfrac_ts )
	#
	if 'Nc_ic' in ncid.variables:	# Level 4 = SALSA microphysics
		# Cloud droplet number concentration averaged over cloudy columns (#/kg)
		CDNC,CDNC_std=average_scaled(ncid.variables['Nc_ic'][ind_tstart:ind_tend],cfrac_ts)
		#
		# Surface precipitation (kg/m^2/s)
		if ind_tstart < ind_tend:
			prcp = numpy.mean( ncid.variables['rmH2Opr'][ind_tstart:ind_tend] )
			prcp_std = numpy.std( ncid.variables['rmH2Opr'][ind_tstart:ind_tend] )
		else:
			prcp = ncid.variables['rmH2Opr'][ind_tstart]
			prcp_std = -999.
		#
		# Change in in-cloud aerosol+cloud droplet number concentration
		if ncfrac>=2: # Linear fit needs at least two data points
			tt = ncid.variables['time'][ind_tstart:ind_tend] 	# Time (s) vector
			nc = ncid.variables['Nc_ic'][ind_tstart:ind_tend]	# Cloud droplets (domain mean)
			nc += ncid.variables['Na_int'][ind_tstart:ind_tend]	# + interstitial aerosol (domain mean)
			# Normalize by cloud fraction => concentrations for cloudy columns
			i=0
			for cf in cfrac_ts:
				if cf>0:
					nc[i]/=cf
				else:
					# Divide-by-zero => NaN
					nc[i]=float('nan')
				i+=1
			#
			a,dn,a_std,dn_std=ls_fit(tt,nc)	# Least squares fit (nc=a+b*tt)
		else:
			dn=-999.
			dn_std=-999.
	else:		# Level 3 = saturation adjustment method (given CDNC)
		# Cloud droplet number concentration (#/kg): fixed
		if ind_tstart < ind_tend:
			CDNC = numpy.mean( ncid.variables['CCN'][ind_tstart:ind_tend] )
			CDNC_std = numpy.std( ncid.variables['CCN'][ind_tstart:ind_tend] )
		else:
			CDNC = ncid.variables['CCN'][ind_tstart]
			CDNC_std = -999.
		#
		# Surface precipitation (kg/m^2/s): variable prcp is in W/m^2=J/s/m^2, which can be
		# converted to kg using latent heat of vaporization (2.5e+06 J/kg)
		if ind_tstart < ind_tend:
			prcp = numpy.mean( ncid.variables['prcp'][ind_tstart:ind_tend] )/2.5e6
			prcp_std = numpy.std( ncid.variables['prcp'][ind_tstart:ind_tend] )/2.5e6
		else:
			prcp = ncid.variables['prcp'][ind_tstart]/2.5e6
			prcp_std = -999.
		#
		# Change in in-cloud aerosol+cloud droplet number concentration: N/A
	#
	# Entrainment velocity (m/s)
	#	we=dz/dt+D*z, where z is PBL height and D is large scale divergence (1.5e-6 1/s) (see e.g. Kazil et al., ACP, 2016).
	if ind_tstart < ind_tend:
		# Must have at least two points for the slope, but should haev more than that
		zz = ncid.variables['zi1_bar'][ind_tstart:ind_tend]	# PBL height (m) vector
		tt = ncid.variables['time'][ind_tstart:ind_tend]	# Time (s) vector
		a,dzdt,a_std,dzdt_std=ls_fit(tt,zz)	# Least squares fit (zz=a+b*tt)
		z=numpy.mean(zz)	# Mean PBL height
		we=dzdt+1.5e-6*z
		we_std=dzdt_std
	else:
		we = -999.
		we_std = -999.
	#
	# Close file
	ncid.close()
	#
	# All done
	return cfrac, CDNC, prcp, dn, we, cfrac_std, CDNC_std, prcp_std, dn_std, we_std,



def get_netcdf_variable(fname,var_name,start_time,end_time=-10000.,ttol=3600.,start_offset=0,end_offset=0):
	# Function for extracting data from a NetCDF file based on the given time value (or range).
	#
	# Inputs:
	#	fname			Complete file path and name
	#	var_name		NetCDF variable name
	#	start_time	Target or start (when end_time is specified) time value
	# Optional inputs
	#	end_time		Optional end time value
	#	ttol			Time tolelance (s) for finding averaging window
	#	start_offset	Point offset to start index (time in the NetCDF files is the save time, so typically should ignore the first points)
	#	end_offset	-||- end index
	#
	# Example:
	#	file='/ibrix/arch/ClimRes/aholaj/case_emulator_DESIGN_v1.4.0_LES_cray.dev20170324_LVL4/emul01/emul01.ts.nc'
	#	lmax=get_netcdf_variable(file,'lmax',3*3600,ttol=10)
	#	lmax=get_netcdf_variable(file,'lmax',2.5*3600,3.5*3600,ttol=10.,start_offset=1)
	import os	
	import numpy
	import netCDF4 as netcdf
	#
	# File must exist
	if not os.path.lexists(fname): raise RuntimeError(fname+' not found!')
	#
	# Open the target NetCDF file
	ncid = netcdf.Dataset(fname,'r')
	#
	if 'time' not in ncid.variables:
		raise RuntimeError('Time not found from '+fname+'!')
	elif var_name not in ncid.variables:
		raise RuntimeError('Variable '+var_name+' not found from '+fname+'!')
	elif 'time' not in ncid.variables[var_name].dimensions:
		raise RuntimeError('Time is not a dimension for '+var_name+' (file '+fname+')!')
	#
	# Time
	times = ncid.variables['time']
	#
	# Find the correct time value
	i=0
	if end_time>-9999.:
		# Time from start_time to end_time (closests matching start and end time values)
		ind_start=0
		ind_end=0
		for tt in times:
			# Closest match
			if abs(tt-start_time)<abs(times[ind_start]-start_time): ind_start=i
			if abs(tt-end_time)<abs(times[ind_end]-end_time): ind_end=i
			i+=1
		#
		if abs(times[ind_start]-start_time)>ttol or abs(times[ind_end]-end_time)>ttol:
			print 'Matching start or end time not found from '+fname+'!'
			return -999.
		#
		# Apply offset (typically the first point is ignored
		ind_start+=start_offset
		ind_end+=end_offset
		if ind_start<0 or ind_start>ind_end or ind_end>=len(times):
			print 'Invalid data range for '+fname+'!'
			return -999.
		#
		# Average over time dimension
		ind=ncid.variables[var_name].dimensions.index('time')
		#
		out=numpy.mean( ncid.variables[var_name][ind_start:ind_end,],axis=ind )
		# Could add standard deviations?
		#out_std = numpy.std( ncid.variables[var_name][ind_start:ind_end,],axis=ind )
	else:
		# Single time value (closest match)
		ind=0
		for tt in times:
			# Closest match
			if abs(tt-start_time)<abs(times[ind]-start_time): ind=i
			i=i+1
		#
		if abs(times[ind]-tstart)>ttol:
			print 'Matching time not found from '+fname+'!'
			return -999.
		#
		# Apply offset (typically the first point is ignored
		ind+=start_offset
		if ind<0 or ind>=len(times):
			print 'Invalid index for '+fname+'!'
			return -999.
		#
		out=ncid.variables[var_name][ind,]
	#
	# Close file
	ncid.close()
	return out



def extract_write_data(fname_template,specs,name_out='',nmax=200,skip_errs=False):
	# Extract and process data from one or more NetCDF files, and write it to a text file (optional)
	#
	# Inputs:
	#	fname_template	File name template with complete path
	#						e.g. '/ibrix/arch/ClimRes/aholaj/case_emulator_DESIGN_v1.4.0_LES_cray.dev20170324_LVL4/emul%02u/emul%02u.ts.nc'
	#	specs				List of variables including slizing and numpy operations
	#	name_out			Output file name (optional)
	# 	nmax				Maximum number of files (optional)
	#	skip_errs			Don't stop on errors - needed when complete data set is not available (saves just NaN)
	#
	# Examples:
	#	fname_template='/ibrix/arch/ClimRes/aholaj/case_emulator_DESIGN_v1.4.0_LES_cray.dev20170324_LVL4/emul%02u/emul%02u.ts.nc'
	#	specs=['cfrac[10]','wmax[10]','lwp_bar[10]']
	#	aa=extract_write_data(fname_template,specs,name_out='tmp.dat')
	#	specs=['cfrac']
	#	aa=extract_write_data(fname_template,specs,name_out='tmp.dat')
	#
	#	fname_template='/ibrix/arch/ClimRes/aholaj/case_emulator_DESIGN_v1.4.0_LES_cray.dev20170324_LVL4/emul%02u/emul%02u.ps.nc'
	#	specs=['numpy.amax(l[10,:])']
	#	aa=extract_write_data(fname_template,specs,name_out='tmp.dat')
	#
	#	fname_template='/ibrix/arch/ClimRes/aholaj/case_emulator_DESIGN_v1.4.0_LES_cray.dev20170324_LVL4/emul%02u/emul%02u.nc'
	#	specs=['numpy.amax(numpy.amax(numpy.amax(l[2,:],axis=0),axis=0),axis=0)']
	#	aa=extract_write_data(fname_template,specs,name_out='tmp.dat')
	#
	import os
	import netCDF4 as netcdf
	import numpy
	import sys
	#
	# Function for converting command line commands to NetCDF format
	def interpret_fun(cmd):
		# Interpret function call, e.g. 'numpy.amax(l[89,:])': just replace variable name x with "ncid.variables['x']",
		# e.g. 'numpy.amax(l[2,100,100,:])
		frst=-1
		lst=-1
		i=0
		for tt in cmd:
			if (tt=='[' or tt==')' or tt==',') and lst==-1:
				# e.g. 'numpy.amax(l[89,:])', 'numpy.amax(l)' or 'numpy.amax(P_Rwca,axis=0)'
				lst=i
			elif tt=='(':
				frst=i
			i+=1
		# Return complete command as fun
		fun=cmd[:frst+1]+'ncid.variables[\''+cmd[frst+1:lst]+'\']'+cmd[lst:]
		return fun
	#
	#
	# Any '%' in file name template will be replaced by counter i=1,2,3,...
	n=fname_template.count('%')
	if n==0: nmax=1	# Template as is
	#
	# Output to text file
	if len(name_out): fid_out=open(name_out,'w')
	#
	nerr=0
	files=0		# Count files
	values=0	# Count values extracted
	out=[] 		# Complete output
	for i in range(1,nmax):
		# File name with full path
		if n==0:
			file_name=fname_template
		elif n==1:
			file_name=fname_template % (i)
		elif n==2:
			file_name=fname_template % (i,i)
		elif n==3:
			file_name=fname_template % (i,i,i)
		else:
			# No more than three?
			file_name=fname_template % (i,i,i,i)
		#
		ncid=0
		if not os.path.isfile(file_name):
			if i==1 and n>0: print file_name+' not found!'
			if not skip_errs or i>90:
				break
			else:
				# Ignore missing file (<90)
				ncid=-999
				nerr+=1
				msg=file_name+' not found!'
				#
				#
				#row=[] # one row
				#for nam in specs:
				#	row.append(obj)
				#
				#out.append(row)
				#
				# Save data
				#if len(name_out):
				#	# Save variables as space separated strings
				#	if not hasattr(obj, '__iter__'):
				#		# Scalar
				#		fid_out.write( str(obj)+" " )
				#		values+=1
				#	else:
				#		# Vector
				#		for val in obj:
				#			fid_out.write( str(val)+" " )
				#			values+=1
				#
				#continue
		#
		# Open input file
		if ncid==0: ncid = netcdf.Dataset(file_name,'r')
		#
		# Outputs
		row=[] # one row
		for nam in specs:
			# Interpret command
			if ncid<0:
				# File not found
				obj = float('nan')	# Scalar works quite often
			elif '(' in nam:
				# There is a call to a function
				fun=interpret_fun(nam)
				try:
					obj =eval( fun )
				except:
					if not skip_errs:
						print "Unexpected error:", sys.exc_info()[0]
						raise
					#
					# Ignore errors
					obj = float('nan')	# Scalar works quite often
					nerr+=1
					msg=sys.exc_info()[0]
			elif '[' in nam:
				# Selected data range
				name=nam[:nam.index('[')]
				ind=nam[nam.index('['):]
				try:
					obj =eval( 'ncid.variables[\''+name+'\']'+ind )
				except:
					if not skip_errs:
						print "Unexpected error:", sys.exc_info()[0]
						raise
					#
					# Ignore errors
					obj = float('nan')	# Scalar works quite often
					nerr+=1
					msg=sys.exc_info()[0]
			else:
				# Data as is
				try:
					obj = ncid.variables[nam][:]
				except:
					if not skip_errs:
						print "Unexpected error:", sys.exc_info()[0]
						raise
					#
					# Ignore errors
					obj = float('nan')	# Scalar works quite often
					nerr+=1
					msg=sys.exc_info()[0]
			#
			# Append data
			row.append(obj)
			#
			# Save data
			if len(name_out):
				# Save variables as space separated strings
				if not hasattr(obj, '__iter__'):
					# Scalar
					fid_out.write( str(obj)+" " )
					values+=1
				else:
					# Vector/matrix
					for val in obj:
						if not hasattr(val, '__iter__'):
							# Scalar (1D vector)
							fid_out.write( str(val)+" " )
							values+=1
						else:
							# Vector (2D matrix)
							for val2 in val:
								fid_out.write( str(val2)+" " )
								values+=1
		#
		# New line
		if len(name_out):	fid_out.write( "\r\n" )
		#
		out.append(row)
		#
		if ncid>0: ncid.close()
		files+=1
	#
	if len(name_out):
		fid_out.close()
		print str(files)+' files examined, '+str(values)+' values saved to '+name_out
		if nerr>0: print '   '+str(nerr)+' error(s) ignored: ',msg
	#
	# Return the data
	return out


#
# LES inputs and outputs
#

def calc_cloud_base(p_surf,theta,rw):
	# Calulate cloud base heigh when liquid water potential temperature (theta [K]) and water
	# vapor mixing ratio (rw [kg/kg]) are constants. Surface pressure p_surf is given in Pa.
	# For more information, see "lifted condensation level" (LCL).
	#
	# Constants
	R=287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
	Rm=461.5	# -||- for water
	ep2=Rm/R-1.0 #M_air/M_water-1
	cp=1005.0	# Specific heat for a constant pressure
	rcp=R/cp
	cpr=cp/R
	g=9.8
	p00=1.0e+05
	#
	# Integrate to cloud base altitude
	dz=1.			# 1 m resolution
	z=0.				# The first altitude
	press=p_surf	# Start from surface
	RH=0
	while RH<100 and z<10000:
		# Temperature (K)
		tavg=theta*(press/p00)**rcp
		#
		# Current RH (%)
		RH=calc_rh(rw,tavg,press)
		if RH>100: break
		#
		# From z to z+dz
		z+=dz
		# Virtual temperature: T_virtual=T*(1+ep2*rl)
		xsi=(1+ep2*rw)
		# Pressure (Pa)
		press-=g*dz*press/(R*tavg*xsi)
	#
	# No cloud
	if RH<100: return -999
	#
	# Return cloud base altitude
	return z


def calc_lwc_altitude(p_surf,theta,rw,zz):
	# Calculate cloud water mixing ratio at a given altitude z (m) when liquid water potential 
	# temperature (theta [k]) and water vapor mixing ratio (rw [kg/kg]) are constants. 
	# Surface pressure p_surf is given in Pa.
	#
	# Constants
	R=287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
	Rm=461.5	# -||- for water
	ep2=Rm/R-1.0 #M_air/M_water-1
	cp=1005.0	# Specific heat for a constant pressure
	rcp=R/cp
	cpr=cp/R
	g=9.8
	p00=1.0e+05
	alvl = 2.5e+06 #  ! latent heat of vaporization
	#
	# a) Integrate to cloud base altitude
	dz=1.			# 1 m resolution
	z=0.				# The first altitude
	press=p_surf	# Start from surface
	RH=0
	while z<zz:
		# Temperature (K) 
		tavg=theta*(press/p00)**rcp
		#
		# Current RH (%)
		RH=calc_rh(rw,tavg,press)
		if RH>100: break
		#
		# From z to z+dz
		z+=dz
		# Virtual temperature: T_virtual=T*(1+ep2*rl)
		xsi=(1+ep2*rw)
		# Pressure (Pa)
		press-=g*dz*press/(R*tavg*xsi)
	#
	# No cloud or cloud water
	if RH<100: return 0.0
	#
	# b) Integrate up to given altitude
	while z<zz:
		# From z to z+dz
		z+=dz
		#
		# Moist adiabatic lapse rate
		q_sat=calc_sat_mixr(press,tavg)
		tavg-=g*(1+alvl*q_sat/(R*tavg))/(cp+alvl**2*q_sat/(Rm*tavg**2))*dz
		#
		# New pressure
		xsi=(1+ep2*q_sat)
		press-=g*dz*press/(R*tavg*xsi)
	#
	# Return cloud water mixing ratio = totol - vapor
	return rw-q_sat


def solve_rw(p_surf,theta,lwc,zz):
	# Solve total water mixing ratio (rw, kg/kg) from surface pressure (p_surf, Pa), liquid water potential
	# temperature (theta, K) and liquid water mixing ratio (lwc) at altitude zz (m)
	#
	# Constants
	R=287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
	Rm=461.5	# -||- for water
	ep2=Rm/R-1.0 #M_air/M_water-1
	cp=1005.0	# Specific heat for a constant pressure
	rcp=R/cp
	cpr=cp/R
	g=9.8
	p00=1.0e+05
	alvl = 2.5e+06 #  ! latent heat of vaporization
	#
	# Mimimum water vapor mixing ratio is at least lwc
	q_min=lwc
	#
	# Maximum water vapor mixing ratio is unlimited, but should be smaller
	# than that for a cloud which base is at surface
	t_surf=theta*(p_surf/p00)**rcp
	q_max=calc_sat_mixr(p_surf,t_surf)
	#
	k=0
	while k<100:
		q_new=(q_min+q_max)/2
		lwc_calc=calc_lwc_altitude(p_surf,theta,q_new,zz)
		#	
		if abs(lwc-lwc_calc)<1e-7:
			break
		elif lwc<lwc_calc:
			q_max=q_new
		else:
			q_min=q_new
		k+=1
		# Failed
		if k==50: return -999
	#
	return q_new




#
#
#
# ================ Helper functions ================
#
def ls_fit(xx,yy):
	# Simple linear least squares fit: y=a+b*x
	import numpy
	#
	# Ignore NaN's
	x=[]; y=[]
	i=0
	for val in xx:
		if not (numpy.isnan(xx[i]) or numpy.isnan(yy[i])):
			x.append(xx[i])
			y.append(yy[i])
		i+=1
	#
	if len(x)<=1:
		# Scalar
		a=0.0; a_std=0.0
		b=1.0; b_std=0.0
	else:
		# Matrix H
		H = numpy.matrix( numpy.vstack([numpy.ones(len(x)),x]).T )
		# LS solution
		th=numpy.linalg.inv( H.T*H ) * H.T *numpy.matrix(y).T
		# Outputs
		a=numpy.asscalar(th[0])
		b=numpy.asscalar(th[1])
		# Parameter uncertainty
		if len(x)>2:
			# Variance
			sv2=((numpy.matrix(y).T-H*th).T * (numpy.matrix(y).T-H*th) )/(len(x)-2)
			std=numpy.sqrt( numpy.asscalar(sv2) * numpy.diagonal( numpy.linalg.inv( H.T*H ) ) )
			# Outputs
			a_std=numpy.asscalar(std[0])
			b_std=numpy.asscalar(std[1])
		else:
			a_std=0.0
			b_std=0.0
	# 
	return a,b,a_std,b_std,

def average_scaled(x,y):
	# Calculate average of x/y so that points where y=0 are ignored
	import numpy
	sx=0.
	sx2=0.
	n=0
	i=0
	for yy in y:
		if yy>0.:
			sx+=x[i]/yy
			sx2+=(x[i]/yy)**2
			n+=1
		i+=1
	#
	if n==0:
		return -999., -999.
	elif n==1:
		return sx, -999.
	else:
		return sx/n, numpy.sqrt( sx2/n - (sx/n)**2 )

#
# Functions from the LES model
# 

def calc_psat_w(T):
	# Function calculates the saturation vapor pressure (Pa) of liquid water as a function of temperature (K)
	#
	# thrm.f90:  real function rslf(p,t)
	c0=0.6105851e+03
	c1=0.4440316e+02
	c2=0.1430341e+01
	c3=0.2641412e-01
	c4=0.2995057e-03
	c5=0.2031998e-05
	c6=0.6936113e-08
	c7=0.2564861e-11
	c8=-.3704404e-13
	#
	x=max(-80.,T-273.16)
	return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))


def calc_sat_mixr(p,T):
	# Function calculates saturation mixing ratio for water (kg/kg)
	#
	# thrm.f90: real function rslf(p,t)
	#
	# r=m_w//m_air
	# R/Rm=287.04/461.5=.622
	#
	esl=calc_psat_w(T)
	return .622*esl/(p-esl)


def calc_rh(rw,T,press):
	# Calculate RH (%) from water vapor mixing ratio rw (r=m_w/m_air [kg/kg]), temperature (K) and pressure (Pa)
	#
	# r=m_w//m_air=pw/Rm/(pair/R)=pw/(p-pw)*R/Rm => pw=p*r/(R/Rm+r)
	#
	R=287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
	Rm=461.5	# Specific gas constant for water
	ep=R/Rm
	#
	psat=calc_psat_w(T)
	return press*rw/(ep+rw)/psat*100
	# When ep>>rw => RH=press*rw/(ep*psat)*100
# ================================