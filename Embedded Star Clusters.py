###########################################################################
###########################################################################
###########################################################################
'''
Part 3.1 & 3.2
Including the potential from a gas Plummer Sphere
'''
###########################################################################
###########################################################################
###########################################################################

'''
The nbody.py code already has these code lines in the main script but they
are commented. We will want to uncomment these line code in order to
be able to run the whole simulation including the gas potential.
'''

########################## First commented line ##########################

# gas potential params (assuming Msim=1 Msol, Lsim=1) 
fgas = 0.5 # a fraction (0.0-1.0)
mplum = 1000.0 #msol total=gas+star
rplum = 1.0 #pc

########################## Second commented line ##########################

# adjust particles masses depending on gas fraction
	mass[:,0] = (1.0-fgas)*np.array(infile['m'])

########################## Third commented line ##########################

# new gas plummer accel
rad = np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
accterm = G * mplum * fgas / (rplum*rplum*rplum*(1.0+rad*rad/(rplum*rplum))**1.5)
acc[:,0] = acc[:,0] - accterm*pos[:,0]
acc[:,1] = acc[:,1] - accterm*pos[:,1]
acc[:,2] = acc[:,2] - accterm*pos[:,2]

