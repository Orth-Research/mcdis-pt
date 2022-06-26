import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np



# latex options
from matplotlib import rc
rc('font', **{'family':'sans-serif', 'size' : 12})
#,'sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('font', family='serif', serif='Times')
rc('text', usetex=False)

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)

# width as measured in inkscape
# width = 3.375
# height = width / 1.618
width = 8
height = 8/1.618
#height = width/1.618

#fig, ax = plt.subplots()
#fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)

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
input_filename = "Data-Tc-J33_0.2-Mod.dat"
input_file = open(input_filename, "r")
x, Tc = np.genfromtxt(input_file, comments = '#', usecols = (0, 1), unpack = True)
input_file.close

# Load FM susceptibility data
# x = 0.5, chi_FM
input_filename = "../..//L_08/Output/Data-AveragedObservables-L_8-Lz_8-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_8, chi_FM_x_5_L_8, chi_FM_error_x_5_L_8 = np.genfromtxt(input_file, comments = '#', usecols = (3, 12, 13), unpack = True)
input_file.close

input_filename = "../..//L_12/Output/Data-AveragedObservables-L_12-Lz_12-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_12, chi_FM_x_5_L_12, chi_FM_error_x_5_L_12 = np.genfromtxt(input_file, comments = '#', usecols = (3, 12, 13), unpack = True)
input_file.close

input_filename = "../..//L_16/Output/Data-AveragedObservables-L_16-Lz_16-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_16, chi_FM_x_5_L_16, chi_FM_error_x_5_L_16 = np.genfromtxt(input_file, comments = '#', usecols = (3, 12, 13), unpack = True)
input_file.close

# x = 0.2, chi_FM
input_filename = "../..//L_08/Output/Data-AveragedObservables-L_8-Lz_8-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_8, chi_FM_x_2_L_8, chi_FM_error_x_2_L_8 = np.genfromtxt(input_file, comments = '#', usecols = (3, 12, 13), unpack = True)
input_file.close

input_filename = "../..//L_12/Output/Data-AveragedObservables-L_12-Lz_12-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_12, chi_FM_x_2_L_12, chi_FM_error_x_2_L_12 = np.genfromtxt(input_file, comments = '#', usecols = (3, 12, 13), unpack = True)
input_file.close

input_filename = "../..//L_16/Output/Data-AveragedObservables-L_16-Lz_16-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_16, chi_FM_x_2_L_16, chi_FM_error_x_2_L_16 = np.genfromtxt(input_file, comments = '#', usecols = (3, 12, 13), unpack = True)
input_file.close

# Load FM magnetization data
# x = 0.5, m^2_FM
input_filename = "../..//L_08/Output/Data-AveragedObservables-L_8-Lz_8-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_8, m_FM_x_5_L_8, m_FM_error_x_5_L_8 = np.genfromtxt(input_file, comments = '#', usecols = (3, 8, 9), unpack = True)
input_file.close

input_filename = "../..//L_12/Output/Data-AveragedObservables-L_12-Lz_12-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_12, m_FM_x_5_L_12, m_FM_error_x_5_L_12 = np.genfromtxt(input_file, comments = '#', usecols = (3, 8, 9), unpack = True)
input_file.close

input_filename = "../..//L_16/Output/Data-AveragedObservables-L_16-Lz_16-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_16, m_FM_x_5_L_16, m_FM_error_x_5_L_16 = np.genfromtxt(input_file, comments = '#', usecols = (3, 8, 9), unpack = True)
input_file.close

# x = 0.2, m^2_FM
input_filename = "../..//L_08/Output/Data-AveragedObservables-L_8-Lz_8-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_8, m_FM_x_2_L_8, m_FM_error_x_2_L_8 = np.genfromtxt(input_file, comments = '#', usecols = (3, 8, 9), unpack = True)
input_file.close

input_filename = "../..//L_12/Output/Data-AveragedObservables-L_12-Lz_12-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_12, m_FM_x_2_L_12, m_FM_error_x_2_L_12 = np.genfromtxt(input_file, comments = '#', usecols = (3, 8, 9), unpack = True)
input_file.close

input_filename = "../..//L_16/Output/Data-AveragedObservables-L_16-Lz_16-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_16, m_FM_x_2_L_16, m_FM_error_x_2_L_16 = np.genfromtxt(input_file, comments = '#', usecols = (3, 8, 9), unpack = True)
input_file.close

# Load AF magnetization data

# x = 0.2, m^2_FM
input_filename = "../..//L_08/Output/Data-AveragedObservables-L_8-Lz_8-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_8, m_AF_x_2_L_8, m_AF_error_x_2_L_8 = np.genfromtxt(input_file, comments = '#', usecols = (3, 24, 25), unpack = True)
input_file.close

input_filename = "../..//L_12/Output/Data-AveragedObservables-L_12-Lz_12-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_12, m_AF_x_2_L_12, m_AF_error_x_2_L_12 = np.genfromtxt(input_file, comments = '#', usecols = (3, 24, 25), unpack = True)
input_file.close


# Load AF susceptibility data
# x = 0.5, chi_AF
input_filename = "../..//L_08/Output/Data-AveragedObservables-L_8-Lz_8-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_8, chi_AF_x_5_L_8, chi_AF_error_x_5_L_8 = np.genfromtxt(input_file, comments = '#', usecols = (3, 26, 27), unpack = True)
input_file.close

input_filename = "../..//L_12/Output/Data-AveragedObservables-L_12-Lz_12-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_12, chi_AF_x_5_L_12, chi_AF_error_x_5_L_12 = np.genfromtxt(input_file, comments = '#', usecols = (3, 26, 27), unpack = True)
input_file.close

input_filename = "../..//L_16/Output/Data-AveragedObservables-L_16-Lz_16-disorder_filling_0.500000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_5_L_16, chi_AF_x_5_L_16, chi_AF_error_x_5_L_16 = np.genfromtxt(input_file, comments = '#', usecols = (3, 26, 27), unpack = True)
input_file.close

# x = 0.2, chi_AF
input_filename = "../..//L_08/Output/Data-AveragedObservables-L_8-Lz_8-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_8, chi_AF_x_2_L_8, chi_AF_error_x_2_L_8 = np.genfromtxt(input_file, comments = '#', usecols = (3, 26, 27), unpack = True)
input_file.close

input_filename = "../..//L_12/Output/Data-AveragedObservables-L_12-Lz_12-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_12, chi_AF_x_2_L_12, chi_AF_error_x_2_L_12 = np.genfromtxt(input_file, comments = '#', usecols = (3, 26, 27), unpack = True)
input_file.close

input_filename = "../..//L_16/Output/Data-AveragedObservables-L_16-Lz_16-disorder_filling_0.200000-J33_0.200000-J34_-1.000000-J44_-1.000000-TMin_0.500000-TMax_6.000000-TSteps_60.dat"
input_file = open(input_filename, "r")
T_x_2_L_16, chi_AF_x_2_L_16, chi_AF_error_x_2_L_16 = np.genfromtxt(input_file, comments = '#', usecols = (3, 26, 27), unpack = True)
input_file.close


# plot(x, y, color='green', marker='^', linestyle='dashed', linewidth=2, markersize=12)

fig = plt.figure(1,figsize = [width,height])

##############################
######## Axes 1 ##############
##############################

ax1 = plt.subplot(2,2,1)
#ax1.plot([0.2, 0.2], linestyle='--', color=color_orange)
ax1.plot(x, Tc, marker='o', color=color_red)
ax1.axvline(x=0.13307573983242402, linestyle='--', color=color_brown)
ax1.axvline(x=0.07913688344848963, linestyle='--', color=color_grey)
#ax1.plot(x, filled_J34 + filled_J44, marker='^', color=color_blue, label = '$J_{34} + J_{44}$') 
#ax1.plot(x, filled_bonds_all, marker='s', color=color_pink, label = 'All')

#ax1.set_title('')
ax1.set_xlabel('$x$ (Co$^{4+}$)')
ax1.set_ylabel('$T/J_{44}$')
ax1.grid(linestyle=':', linewidth=1.)
ax1.set_xlim(left=0.)
ax1.set_ylim(bottom=0.)

#ax1.legend()
#loc='upper right')
#, bbox_to_anchor=(0.5, 1.05), ncol=1)

##############################
######## Axes 2 ##############
##############################

ax2 = plt.subplot(2,2,2)
#ax2.plot([0.2, 0.2], linestyle='--', color=color_orange)
ax2.plot(T_x_5_L_8, chi_FM_x_5_L_8, marker='.', color=color_red, label = '$x=0.5,L=8$')
ax2.plot(T_x_5_L_12, chi_FM_x_5_L_12, marker='^', color=color_red, label = '$x=0.5,L=12$')

ax2.plot(T_x_2_L_8, chi_FM_x_2_L_8, marker='.', color=color_orange, label = '$x=0.2,L=8$')
ax2.plot(T_x_2_L_12, chi_FM_x_2_L_12, marker='^', color=color_orange, label = '$x=0.2,L=12$')
#ax2.set_xlim(right=9)

# ax2.plot(x, filled_J34 + filled_J44, marker='^', color=color_blue, label = '$J_{34} + J_{44}$') 
# ax2.plot(x, filled_bonds_all, marker='s', color=color_pink, label = 'All')

# #ax2.set_title('')
ax2.set_xlabel('$T/J_{44}$')
ax2.set_ylabel('$\\chi_{\mathrm{FM}}$')
# ax2.grid(linestyle=':', linewidth=1.)
# ax2.set_xlim(0., 0.15)
# ax2.set_ylim(0., 0.3)
# ax2.set_xticks([0, 0.05, 0.1, 0.15])
#ax2.legend()


##############################
######## Axes 3 ##############
##############################

ax3 = plt.subplot(2,2,3)
#ax3.plot([0.2, 0.2], linestyle='--', color=color_orange)
ax3.set_yscale('log')
#ax3.plot(T_x_5_L_8, m_FM_x_5_L_8, marker='.', color=color_red, label = '$x=0.5,L=8$')
ax3.plot(T_x_5_L_12, m_FM_x_5_L_12, marker='^', color=color_red, label = 'FM')

#ax3.plot(T_x_2_L_8, m_FM_x_2_L_8, marker='.', color=color_orange, label = '$x=0.2,L=8$')
ax3.plot(T_x_2_L_12, m_FM_x_2_L_12, marker='^', color=color_orange, label = 'FM')

#ax3.plot(T_x_2_L_8, m_AF_x_2_L_8, marker='.', color=color_green, label = '$x=0.2,L=8$')
ax3.plot(T_x_2_L_12, m_AF_x_2_L_12, marker='v', color=color_orange, label = 'AF')
#ax3.set_xlim(right=9)

# ax3.plot(x, filled_J34 + filled_J44, marker='^', color=color_blue, label = '$J_{34} + J_{44}$') 
# ax3.plot(x, filled_bonds_all, marker='s', color=color_pink, label = 'All')

# #ax3.set_title('')
ax3.set_xlabel('$T/J_{44}$')
ax3.set_ylabel('$m^2$')
# ax3.grid(linestyle=':', linewidth=1.)
# ax3.set_xlim(0., 0.15)
# ax3.set_ylim(0., 0.3)
# ax3.set_xticks([0, 0.05, 0.1, 0.15])

ax3.legend()

##############################
######## Axes 4 ##############
##############################
ax4 = plt.subplot(2,2,4)
#ax4.plot([0.2, 0.2], linestyle='--', color=color_orange)
ax4.plot(T_x_5_L_8, chi_AF_x_5_L_8, marker='.', color=color_red, label = '$x=0.5,L=8$')
ax4.plot(T_x_5_L_12, chi_AF_x_5_L_12, marker='^', color=color_red, label = '$x=0.5,L=12$')

ax4.plot(T_x_2_L_8, chi_AF_x_2_L_8, marker='.', color=color_orange, label = '$x=0.2,L=8$')
ax4.plot(T_x_2_L_12, chi_AF_x_2_L_12, marker='^', color=color_orange, label = '$x=0.2,L=12$')
#ax4.set_xlim(right=9)

# ax4.plot(x, filled_J34 + filled_J44, marker='^', color=color_blue, label = '$J_{34} + J_{44}$') 
# ax4.plot(x, filled_bonds_all, marker='s', color=color_pink, label = 'All')

# #ax4.set_title('')
ax4.set_xlabel('$T/J_{44}$')
ax4.set_ylabel('$\\chi_{\mathrm{AF}}$')
# ax4.grid(linestyle=':', linewidth=1.)
# ax4.set_xlim(0., 0.15)
# ax4.set_ylim(0., 0.3)
# ax4.set_xticks([0, 0.05, 0.1, 0.15])
ax4.legend()


plt.tight_layout()

fig.savefig("Plot-Phase_Diagram.pdf", bbox_inches="tight")
plt.show()

