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
