###########################################################################
###########################################################################
###########################################################################
'''
Part 1.1.1
Setting up a Plummer Sphere
'''
###########################################################################
###########################################################################
###########################################################################

import matplotlib.pyplot as plt
n = 17
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    'font.size': n,         #Tamaño de fuente general
    'axes.titlesize': n,    #Tamaño de fuente para títulos de ejes
    'axes.labelsize': n,    #Tamaño de fuente para etiquetas de ejes
    'xtick.labelsize': n,   #Tamaño de fuente para etiquetas del eje x
    'ytick.labelsize': n,   #Tamaño de fuente para etiquetas del eje y
    'legend.fontsize': n,   #Tamaño de fuente para la leyenda
    'figure.titlesize': n   #Tamaño de fuente para el título de la figura
})
import numpy as np

# Fixing random state for reproducibility
rng = np.random.default_rng()

# Arrays to store values
x1_vals = []
x2_vals = []
x3_vals = []

x_value = []
y_value = []
z_value = []

M_e = []
radii = []

# Plotting
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

n = 1000                    # number of galaxies
rmax = 1                    # kpc
m = '.'                     # marker symbol
r_s = 1                     # kpc
Mo = 1                      # solar mass

for i in range(n):
    x1 = rng.random()
    x2 = rng.random()
    x3 = rng.random()

    # Store the values
    x1_vals.append(x1)
    x2_vals.append(x2)
    x3_vals.append(x3)

    rs = x3 * rmax
    #rs = 1
    z = (1 - 2 * x1) * rs
    x = np.sqrt(rs**2 - z**2) * np.cos(2 * np.pi * x2)
    y = np.sqrt(rs**2 - z**2) * np.sin(2 * np.pi * x2)
    r = np.sqrt(x**2 + y**2 + z**2)

    ax.scatter(x, y, z, marker=m)

    M = (Mo*rs**3) / (rs**2 + r_s**2)**(1.5)
    
    M_e.append(M)
    radii.append(r)
    x_value.append(x)
    y_value.append(y)
    z_value.append(z)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.tight_layout()
plt.show()

# Plot histograms of x1, x2, x3
fig, axs = plt.subplots(1, 3, figsize=(15, 4))

axs[0].hist(x1_vals, bins=30, color='skyblue', edgecolor='black')
axs[0].set_title('Histogram of x1')

axs[1].hist(x2_vals, bins=30, color='salmon', edgecolor='black')
axs[1].set_title('Histogram of x2')

axs[2].hist(x3_vals, bins=30, color='lightgreen', edgecolor='black')
axs[2].set_title('Histogram of x3')

for ax in axs:
    ax.set_xlabel('Value')
    ax.set_ylabel('Frequency')

plt.tight_layout()
plt.show()

###########################################################################
###########################################################################
###########################################################################
'''
Part 1.1.2
Half total mass radius
'''
###########################################################################
###########################################################################
###########################################################################
# Necessary initial condition
r_pl = 1

# Arrays to store values
radii = []
Mass = []

# Iteration loop
for j in range(0,n,1):
    x1 = rng.random()
    r = r_pl*((x1**(-2/3) - 1))**(-1/2)
    Mass.append(x1)
    radii.append(r)

# Sorting the radii array
r_half = np.sort(radii)

'''Since there is only 1000 stars and each one
with the same solar mass, the 500th position
corresponds to the half mass radius'''
print(f'the radius containing half the total mass is {r_half[500]:.2}')

# Plotting
fig, axs = plt.subplots(1, 2, figsize=(15, 4))

axs[0].plot(radii)
axs[0].set_title('original radii arrangement')

axs[1].plot(r_half)
axs[1].set_title('sorted radii arrangement')

plt.show()

# Plotting the mass
plt.scatter(r_half, np.sort(Mass))
plt.xlabel('sorted radii')
plt.ylabel('sorted mass enclosed')
plt.show()

###########################################################################
###########################################################################
###########################################################################
'''
Part 1.2
Velocity distribution
'''
###########################################################################
###########################################################################
###########################################################################
from astropy import units as u
from astropy import constants as const

M_pl = 1e3  
G = const.G.value

# Arrays to store values
V_esc = []
V_tot = []
v = []
radius = []
vx_values = []
vy_values = []
vz_values = []

# Iterable loop
k=0
while k<=999:
    x1 = rng.random()
    x2 = rng.random()
    if 0.1*x2 < x1**2 * (1 - x1**2)**(3.5):
        vel = (2*G*M_pl/r_pl)**(0.5)*(1 + (radii[k]*2/r_pl**2))**(-1/4)
        v_tot = x1*vel
        
        # Calculation of velocity components
        vz = (1 - 2*x1)*v_tot
        vx = (v_tot**2 - vz**2)**(0.5)*np.cos(2*np.pi*x2)
        vy = (v_tot**2 - vz**2)**(0.5)*np.sin(2*np.pi*x2)
        velocity = np.sqrt(vx**2 + vy**2 + vz**2)
        
        # Storing values
        V_esc.append(vel)
        V_tot.append(v_tot)
        v.append(velocity)        
        vx_values.append(vx)
        vy_values.append(vy)
        vz_values.append(vz)
        radius.append(radii[k])       
        
        plt.scatter(x1,x2*0.1)
        
        k+=1

# Plotting the function velocity distribution
print(f'the total particles considered are {len(V_esc)} from the {n} total')
plt.xlabel('random number x1')
plt.ylabel('random number x2')
plt.title('Accepted values')
plt.hlines(y=0.1, xmin=0, xmax=1, colors='red', linestyles='dashed')
plt.grid(alpha=0.3)
plt.show()

# Plotting the velocities over radius
plt.scatter(radius, v, color='gray', label='particle´s velocity')
plt.scatter(radius, V_esc, color='red', linewidths=0.1, label='scape velocities')
plt.xlabel('radius [kpc]')
plt.ylabel('velocity [km/s]')
plt.xscale('log')
plt.grid(alpha=0.3)
plt.legend()
plt.show()
