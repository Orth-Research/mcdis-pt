import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# latex options
from matplotlib import rc
rc('font', **{'family':'sans-serif', 'size' : 12}) 
rc('text', usetex=False)

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

# width as measured in inkscape
width = 4
height = width

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

# Load data
input_filename = "Data-Magnetization-L_10-Lz_10-J33_0.000000-J34_-0.500000-J44_-1.000000-temperature_0.100000-run_steps_1000.dat"
input_file = open(input_filename, "r")
x_pol, isolated_polarons = np.genfromtxt(input_file, comments = '#', usecols = (1, 4), unpack = True)
input_file.close

fig = plt.figure(1,figsize = [width, height])

ax1 = plt.subplot(1,1,1)
ax1.plot([-0.05, 0.15], [0.31, 0.31], linestyle='--', color=color_orange)
ax1.plot(x, filled_sites, marker='o', color=color_red)
ax1.set_xlabel('$x$ (Co$^{4+}$)')
ax1.set_ylabel('Fraction of filled sites')

ax1.set_xticks([0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6], minor = True)
ax1.grid(linestyle=':', linewidth=1., which = 'both')
ax1.set_xlim(left = -0.05, right=0.6)
ax1.set_ylim(bottom = - 0.02, top = None)

ax2 = inset_axes(ax1, width="35%", height="60%", loc='lower right', borderpad=1.5)
ax2.plot(x_pol, 7*isolated_polarons, marker='.', color=color_blue)

ax2.set_ylabel('# spins in isol. pol./N')
ax2.grid(linestyle=':', linewidth=1., which = 'both')
ax2.set_xlim(left = -0.01, right = 0.2)
ax2.set_xticks([0, 0.1, 0.2], minor=True)
ax2.set_yticks([0,0.04, 0.08])

fig.savefig("Plot-Fraction_sites.pdf", bbox_inches="tight")
plt.show()

