###########################################################################
###########################################################################
###########################################################################
'''
Part 3.2 & 3.3
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

'''
Be aware that there are two sections that needs to be uncommmented
for the gas Plummer's acceleration. Both has the same information
but are in different code lines.
'''

# new gas plummer accel
rad = np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
accterm = G * mplum * fgas / (rplum*rplum*rplum*(1.0+rad*rad/(rplum*rplum))**1.5)
acc[:,0] = acc[:,0] - accterm*pos[:,0]
acc[:,1] = acc[:,1] - accterm*pos[:,1]
acc[:,2] = acc[:,2] - accterm*pos[:,2]

###########################################################################
###########################################################################
###########################################################################
'''
Part 3.4
Simulating a supernova explosion
'''
###########################################################################
###########################################################################
###########################################################################

'''
To simulate a loss fraction of the gas in a spontaneous time denoted as
after 5 crossing times (optional) we will need to make sure that the
gas fraction is set to 0 after this condition is accomplished.
'''

if t >= tEnd*0.5:
	fgas = 0

'''
this line code needs to be held before the - # new gas plummer accel -
line code (see Part 3.3) that is INSIDE the simulated loop.
'''

###########################################################################
###########################################################################
###########################################################################
'''
Part 3.5
Simulating a supernova explosion
'''
###########################################################################
###########################################################################
###########################################################################


