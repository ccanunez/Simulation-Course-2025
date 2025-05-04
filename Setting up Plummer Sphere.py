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
