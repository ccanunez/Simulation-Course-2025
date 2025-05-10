import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
from astropy.io import ascii

"""
Create Your Own N-body Simulation (With Python)
Philip Mocz (2020) Princeton Univeristy, @PMocz

Simulate orbits of stars interacting due to gravity
Code calculates pairwise forces according to Newton's Law of Gravity
"""

def getAcc( pos, mass, G, softening ):
	"""
    Calculate the acceleration on each particle due to Newton's Law 
	pos  is an N x 3 matrix of positions
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	softening is the softening length
	a is N x 3 matrix of accelerations
	"""
	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r^3 for all particle pairwise particle separations 
	inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)
	inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)

	ax = G * (dx * inv_r3) @ mass
	ay = G * (dy * inv_r3) @ mass
	az = G * (dz * inv_r3) @ mass
	
	# pack together the acceleration components
	a = np.hstack((ax,ay,az))

	return a
	
def getEnergy( pos, vel, mass, G ):
	"""
	Get kinetic energy (KE) and potential energy (PE) of simulation
	pos is N x 3 matrix of positions
	vel is N x 3 matrix of velocities
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	KE is the kinetic energy of the system
	PE is the potential energy of the system
	"""
	# Kinetic Energy:
	KE = 0.5 * np.sum(np.sum( mass * vel**2 ))


	# Potential Energy:

	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r for all particle pairwise particle separations 
	inv_r = np.sqrt(dx**2 + dy**2 + dz**2)
	inv_r[inv_r>0] = 1.0/inv_r[inv_r>0]

	# sum over upper triangle, to count each interaction only once
	PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*inv_r,1)))
	
	return KE, PE;


def main():
	""" N-body simulation """
	
	# Simulation parameters
	N         = 1000     # Number of particles
	t         = 0      # current time of the simulation
	tCross = 0.19784190099458493	# Tsim
	tEnd      = 10.0 * tCross   # time at which simulation ends
	n = 100
	dt        =  tCross / n  # timestep --> dt = t_cross / n_steps
	softening = 0.23777004980940408    # softening length
	'''
	softening = 1/n**(1/3) = (volume / N/2)**(1/3)
	'''
	G         = 1.0    # Newton's Gravitational Constant
	plotRealTime = True # switch on for plotting as the simulation goes along
	
	# open and read input file
	infile = ascii.read("/Users/ccanunez/Documents/USM/9no Semestre/Simulaciones/Project 2/plumparts.csv")
#	everyXsnaps = 10

	# gas potential params (assuming Msim=1 Msol, Lsim=1) 
	fgas = 0.5 # a fraction (0.0-1.0)
	mplum = 1000.0 #msol total=gas+star
	rplum = 1.0 #pc

	
	# open an output file
	f = open('/Users/ccanunez/Documents/USM/9no Semestre/Simulaciones/Project 2/output.dat', 'wb')
	
	# Generate Initial Conditions
	np.random.seed(17)            # set the random number generator seed
	
	mass = np.zeros((N,1))  # mass of particles
	pos = np.zeros((N,3))
	vel = np.zeros((N,3))
	pos[:,0] = np.array(infile['x'])
	pos[:,1] = np.array(infile['y'])
	pos[:,2] = np.array(infile['z'])
	vel[:,0] = np.array(infile['vx'])
	vel[:,1] = np.array(infile['vy'])
	vel[:,2] = np.array(infile['vz'])
	mass[:,0] = np.array(infile['m'])
# adjust particles masses depending on gas fraction
	mass[:,0] = (1.0-fgas)*np.array(infile['m'])

	tmat = np.zeros((N,1))
	ival = np.zeros((N,1))
	
	# Convert to Center-of-Mass frame
	vel -= np.mean(mass * vel,0) / np.mean(mass)
	
	# calculate initial gravitational accelerations
	acc = getAcc( pos, mass, G, softening )

# new gas plummer accel
	rad = np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
	accterm = G * mplum * fgas / (rplum*rplum*rplum*(1.0+rad*rad/(rplum*rplum))**1.5)
	acc[:,0] = acc[:,0] - accterm*pos[:,0]
	acc[:,1] = acc[:,1] - accterm*pos[:,1]
	acc[:,2] = acc[:,2] - accterm*pos[:,2]
	
	# calculate initial energy of system
	KE, PE  = getEnergy( pos, vel, mass, G )
	
	# number of timesteps
	Nt = int(np.ceil(tEnd/dt))
	
	# save energies, particle orbits for plotting trails
	pos_save = np.zeros((N,3,Nt+1))
	pos_save[:,:,0] = pos
	KE_save = np.zeros(Nt+1)
	KE_save[0] = KE
	PE_save = np.zeros(Nt+1)
	PE_save[0] = PE
	t_all = np.arange(Nt+1)*dt
	
	# prep figure
	fig = plt.figure(figsize=(4,5), dpi=80)
	grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
	ax1 = plt.subplot(grid[0:2,0])
	ax2 = plt.subplot(grid[2,0])
	
	for i in range(N):
		ival[i] = int(i)+1
#		print(i,ival[i])
#	print(ival[1],ival[5])	
	
	# Simulation Main Loop
	for i in range(Nt):
		# (1/2) kick
		vel += acc * dt/2.0
		
		# drift
		pos += vel * dt
		
		# update accelerations
		acc = getAcc( pos, mass, G, softening )

		# new gas plum accel
		tGone = tEnd / 2
		tSN = tEnd / 2

		if t >= tEnd/2:
			fgas = -( (t-tSN) / tGone) + 1
		rad = np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
		accterm = G * mplum * fgas / (rplum*rplum*rplum*(1.0+rad*rad/(rplum*rplum))**1.5)		
		acc[:,0] = acc[:,0] - accterm*pos[:,0]
		acc[:,1] = acc[:,1] - accterm*pos[:,1]
		acc[:,2] = acc[:,2] - accterm*pos[:,2]
		
		# (1/2) kick
		vel += acc * dt/2.0
		
		# update time
		t += dt

		tmat[:] = t
		
		# output to ascii file
		new = np.append(pos,ival,axis=1)
		new2 = np.append(new,tmat,axis=1)
		if(i==1): 
			np.savetxt(f,new2, header='x y z id snap')
#		if(i % everyXsnaps==1): 
#			np.savetxt(f,new2)
		if(i>1): 
			np.savetxt(f,new2)

		
		# get energy of system
		KE, PE  = getEnergy( pos, vel, mass, G )
		
		# save energies, positions for plotting trail
		pos_save[:,:,i+1] = pos
		KE_save[i+1] = KE
		PE_save[i+1] = PE
		
		# plot in real time
		if plotRealTime or (i == Nt-1):
			plt.sca(ax1)
			plt.cla()
			xx = pos_save[:,0,max(i-50,0):i+1]
			yy = pos_save[:,1,max(i-50,0):i+1]
			plt.scatter(xx,yy,s=1,color=[.7,.7,1])
			plt.scatter(pos[:,0],pos[:,1],s=10,color='blue')
			ax1.set(xlim=(-10, 10), ylim=(-10, 10))
			ax1.set_aspect('equal', 'box')
			ax1.set_xticks([-10,-5,0,5,10])
			ax1.set_yticks([-10,-5,0,5,10])
			
			plt.sca(ax2)
			plt.cla()
			plt.scatter(t_all,KE_save,color='red',s=1,label='KE' if i == Nt-1 else "")
			plt.scatter(t_all,PE_save,color='blue',s=1,label='PE' if i == Nt-1 else "")
			plt.scatter(t_all,KE_save+PE_save,color='black',s=1,label='Etot' if i == Nt-1 else "")
			ax2.set(xlim=(0, tEnd))#, ylim=(-0.00000000003, 0.00000000003))
			#ax2.set_aspect(0.1)
			
			plt.pause(0.001)
	    
	
	
	# add labels/legend
	plt.sca(ax2)
	plt.xlabel('time')
	plt.ylabel('energy')
	ax2.legend(loc='upper right')
	
	# Save figure
	plt.savefig('/Users/ccanunez/Documents/USM/9no Semestre/Simulaciones/Project 2/nbody.png',dpi=240)
	plt.show()
	
	    
	return 0
	


  
if __name__== "__main__":
  main()
