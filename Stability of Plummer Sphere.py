

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
