###########################################################################
###########################################################################
###########################################################################
'''
Part 2.1 & 2.2 & 2.3
2.1 Convertion of units into G = 1 units
2.2 Calculation of crossing time
2.3 Calculation of softening
'''
###########################################################################
###########################################################################
###########################################################################

'''G = 1 units'''
G = 1
Lsim = 1    # pc
Tsim = 14.9 # Myr
Vsim = 65.8 # m/s
Msim = 1    # solar mass
N = 1e3

'''converting previous parameters into G = 1 units'''
Vel = []
for i in V_tot:
    g_unit_velocity = i/Vsim
    Vel.append(g_unit_velocity/Vsim)

M_pl = M_pl/Msim
r_pl = r_pl/Lsim

'''calculating t cross'''
t_cross = (3*np.pi/32)**(-3/2) * (r_pl**3/(G*M_pl))**(1/2)
print(f'the crossing time is {t_cross} in Tsim')
print(f'and 10 crossing times equals {10*t_cross} in Tsim')
print(f'therefore 10 crossing times in Tsim are {10*Tsim*t_cross:.3} Myr')

'''softening'''
volumn = 4/3 * np.pi * r_half[500]**2
softening = (volumn / (N/2))**(1/3)
print(f'the softening is {softening}')

###########################################################################
###########################################################################
###########################################################################
'''
Part 2.4
Storing values for velocity and position & creation of the dataFrame table
in order to use in the nbody.py code
'''
###########################################################################
###########################################################################
###########################################################################

from astropy import units as u
from astropy import constants as const

# Initial parameters and definition of G
M_pl = 1e3  
G = const.G.value

# Arrays to store values
x_value = []
y_value = []
z_value = []
M_e = []
V_esc = []
V_tot = []
v = []
vx_values = []
vy_values = []
vz_values = []


# Iterable loop for velocity
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
        
        k+=1

# Itetable loop for position
for j in radii:
        # Random sample
        x1 = rng.random()
        x2 = rng.random()

        # Creation of spatial coordinate
        rs = j
        z = (1 - 2 * x1) * rs
        x = np.sqrt(rs**2 - z**2) * np.cos(2 * np.pi * x2)
        y = np.sqrt(rs**2 - z**2) * np.sin(2 * np.pi * x2)
        r = np.sqrt(x**2 + y**2 + z**2)

        # Constant mass for each particle
        M = 1

        # Saving the data
        M_e.append(M)
        x_value.append(x)
        y_value.append(y)
        z_value.append(z)

print(f'the len of pos is {len(x_value)}')
import pandas as pd

# Create a DataFrame
df = pd.DataFrame({
    'x': x_value,
    'y': y_value,
    'z': z_value,
    'm': M_e,
    'vx': vx_values,
    'vy': vy_values,
    'vz': vz_values,
})

# Export to CSV
df.to_csv('plumparts.csv', index=False)

###########################################################################
###########################################################################
###########################################################################
'''
Part 2.5
Plotting the Lagrangian radii
'''
###########################################################################
###########################################################################
###########################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

# Initial parameters in G = 1 units
G = 1
Lsim = 1  # pc
Msim = 1  # Msun
Tsim = 14.9 * np.sqrt(Lsim**3 / Msim)
Vsim = 980.4 * (Lsim / Tsim)

# Files
input_file = 'output.dat'
output_file = 'data.csv'

# Reads the input file in ASCII format
df = pd.read_csv(input_file, delim_whitespace=True)

# Transforms and saves the file in CSV format
df.to_csv(output_file, index=False)
print(f"Converted '{input_file}' to '{output_file}'")

# Columns information extraction
x_sim = df['x'].values
y_sim = df['y'].values
z_sim = df['z'].values
part_id = df['id'].values
snap = df['snap'].values

# Compoute radius
r = np.sqrt(x_sim**2 + y_sim**2 + z_sim**2)
df['r'] = r  # Adds a new column for the radius

# Obtains unique snaps and their first index
snaps, ind_snaps = np.unique(snap, return_index=True)
print(f'the first snapshot is {snaps[0]}')            # First snapshot
print(f'the number of snapshots are {len(ind_snaps)}')      # Number of unique snapshots

# Stacks radius per snapshor
r_arr = []
for s in tqdm(snaps, desc='Processing'):
    pos = np.where(snap == s)[0]
    r_arr.append(r[pos])

# Compute percentiles
r_half, r20, r80 = [], [], []
for r_snapshot in r_arr:
    r_sort = np.sort(r_snapshot)
    N = len(r_sort)
    r_half.append(r_sort[N // 2])
    r20.append(r_sort[int(N * 0.2)])
    r80.append(r_sort[int(N * 0.8)])

# Plotting proccess of radius evolution
plt.plot(snaps * Tsim, r20, c='navy', lw=1.5, label=r'$R_{20}$')
plt.plot(snaps * Tsim, r_half, c='purple', lw=1.5, label=r'$R_{50}$')
plt.plot(snaps * Tsim, r80, c='darkred', lw=1.5, label=r'$R_{80}$')
plt.ylabel(r'$r\ [pc]$', fontsize=15)
plt.xlabel(r'$t\ [Myr]$', fontsize=15)
plt.ylim(0, 4)
plt.legend()
plt.show()
