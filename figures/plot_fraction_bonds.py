import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# latex options
from matplotlib import rc
rc('font', **{'family':'sans-serif', 'size' : 12})
rc('text', usetex=False)

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

# width as measured in inkscape
# width = 3.375
# height = width / 1.618
width = 8
height = 8/2

# define colors
color_red = (0.73, 0.13869999999999993, 0.)
color_orange = (1., 0.6699999999999999, 0.)
color_green = (0.14959999999999996, 0.43999999999999995, 0.12759999999999994)
color_blue = (0.06673600000000002, 0.164512, 0.776)
color_purple = (0.25091600000000003, 0.137378, 0.29800000000000004)
color_ocker = (0.6631400000000001, 0.71, 0.1491)
color_pink = (0.71, 0.1491, 0.44730000000000003)
color_brown = (0.651, 0.33331200000000005, 0.054683999999999955)
color_red2 = (0.766, 0.070, 0.183)
color_turquoise = (0., 0.684, 0.676)
color_yellow = (0.828, 0.688, 0.016)
color_grey = (0.504, 0.457, 0.410)

# Load data
input_filename = "Data-Fraction.dat"
input_file = open(input_filename, "r")
x, filled_sites, filled_sites_error, filled_J33, filled_J33_error, filled_J34, filled_J34_error, filled_J44, filled_44_error, filled_bonds_all, filled_bonds_all_error = np.genfromtxt(input_file, comments = '#', usecols = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), unpack = True)
input_file.close

# Create figure
fig = plt.figure(1,figsize = [width,height])

# Axes 1
ax1 = plt.subplot(1,2,1)
ax1.plot([0.25, 0.25], linestyle='--', color=color_orange)
ax1.plot(x, filled_J33, marker='s', color=color_green, label = '$J_{33}$')
ax1.plot(x, filled_J34, marker='s', color=color_yellow, label = '$J_{34}$')
ax1.plot(x, filled_J44, marker='s', color=color_purple, label = '$J_{44}$')
ax1.plot(x, filled_J34 + filled_J44, marker='s', color=color_blue, label = '$J_{34} + J_{44}$')
ax1.plot(x, filled_bonds_all, marker='s', color=color_pink, label = 'All')

ax1.set_xlabel('$x$ (Co$^{4+}$)')
ax1.set_ylabel('Fraction of filled bonds')
ax1.grid(linestyle=':', linewidth=1.)

# Axes 2
ax2 = plt.subplot(1,2,2)
ax2.plot([0.25, 0.25], linestyle='--', color=color_orange)
ax2.plot(x, filled_J33, marker='s', color=color_green, label = '$J_{33}$')
ax2.plot(x, filled_J34, marker='s', color=color_yellow, label = '$J_{34}$')
ax2.plot(x, filled_J44, marker='s', color=color_purple, label = '$J_{44}$')
ax2.plot(x, filled_J34 + filled_J44, marker='s', color=color_blue, label = '$J_{34} + J_{44}$')
ax2.plot(x, filled_bonds_all, marker='s', color=color_pink, label = 'All')

ax2.set_xlabel('$x$ (Co$^{4+}$)')
ax2.grid(linestyle=':', linewidth=1.)
ax2.set_xlim(0., 0.15)
ax2.set_ylim(0., 0.3)
ax2.set_xticks([0, 0.05, 0.1, 0.15])

ax2.legend()

plt.tight_layout()

fig.savefig("Plot-Fraction_Bonds.pdf", bbox_inches="tight")
plt.show()

