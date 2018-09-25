import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os


def maximum_kappa_CC(T, p):
    """Determine the maximum exchange of oxygen, for cocurrent flows of CO2 and H2 with 1:3 flow rates, omega=3.
    Returns, CO2_conversion =  kappa_max"""
    # define the two flows as canter solution objects using cantera's gri30 database
    flow1 = ct.Solution('gri30.cti')
    flow2 = ct.Solution('gri30.cti')
    # complete conversion of CO2 to CO gives Kappa = 0.5, kappa in range 0-0.5
    kappa = np.arange(0, 0.50, 0.002)
    for k in kappa:
        # iterates through kappa until the CO2 and CH4 oxygen partial pressures meet in CC configuration
        # equal flows of CH4 and CO2, i.e. omega=1
        flow1.X = {'CO': 1, 'O2': 0.5 - k}
        flow2.X = {'H2': 3, 'O2': k}
        flow1.TP = T, p
        flow2.TP = T, p
        flow1.equilibrate('TP')
        flow2.equilibrate('TP')
        pO2_flow1 = flow1.X[flow1.species_index('O2')]*p
        pO2_flow2 = flow2.X[flow2.species_index('O2')]*p
        if (pO2_flow2 > pO2_flow1 or k > 0.497):
            kappa_max = k
            break
    return kappa_max

def maximum_kappa_CT(T, p):
    """Determine the maximum exchange of oxygen for countercurrent flows of CO2 and H2 with 1:3 flow rates, omega=3.
    Returns, CO2_conversion =  2*Kappa_max"""
    # define the two flows as canter solution objects using cantera's gri30 database
    flow1 = ct.Solution('gri30.cti')
    flow2 = ct.Solution('gri30.cti')
    # Make arrays tof the pO2(kappa) in each flow
    PO2_flow2 = []
    PO2_flow1 = []
    # complete conversion of CO2 to CO gives Kappa = 0.5, kappa in range 0-0.5
    kappa = np.arange(0, 0.50, 0.002)
    for k in kappa:
        flow1.X = {'CO': 1, 'O2': 0.5 - k}
        flow2.X = {'H2': 3, 'O2': k}
        flow1.TP = T, p
        flow2.TP = T, p
        flow1.equilibrate('TP')
        flow2.equilibrate('TP')
        pO2_flow1 = flow1.X[flow1.species_index('O2')] * p
        pO2_flow2 = flow2.X[flow2.species_index('O2')] * p
        # add the new pO2 values to the arrays
        PO2_flow1.append(pO2_flow1)
        # for the second flow at pO2 value to the start of array to reverse kappa for countercurrent
        PO2_flow2.insert(0, pO2_flow2)
        # if the arrays meet at any point, then we have reached the maximum exchange extent kappa_max
        if np.any(np.asarray(PO2_flow1)-np.asarray(PO2_flow2) <= 0) or k > 0.497:
            kappa_max = k
            break
    return kappa_max






# Lets look at the thermodynamics as a function of temperature in the range ata pressure of 1 bar,
Trange = np.arange(573, 1183, 10)
p = 100000

# compare results to regular water gas shift with forced selectivity to avoid methanol or methane formation. Use
# modified database gri0met to do this.
RWGS = ct.Solution('gri0met.cti')

# store results in arrays
CC_co2_conversion = np.zeros(len(Trange))
CT_co2_conversion = np.zeros(len(Trange))
RWGS_co2_conversion = np.zeros(len(Trange))



for i, T in enumerate(Trange):
    RWGS.TP = T, p
    RWGS.X = {'H2': 3, 'CO2': 1}
    RWGS.equilibrate('TP')
    RWGS_co2_conversion[i] = RWGS.X[RWGS.species_index('CO')]/(RWGS.X[RWGS.species_index('CO')]+ RWGS.X[RWGS.species_index('CO2')])
    CC_co2_conversion[i] = 2*maximum_kappa_CC(T, p)
    CT_co2_conversion[i] = 2*maximum_kappa_CT(T, p)
    print("T, CC kappa_max, CT kappa_max:", T, CC_co2_conversion[i]/2, CT_co2_conversion[i]/2)


# NEXT PLOT REASULTS
# Some graphing definitions
plt.style.use('seaborn-paper')

axisfont = 16
ticklabelfont = 13
labelfont = 14

def label_line(ax, line, at_x=0, rotation_on=False, fontsize=14, offset=(0, 0), **kwargs):
    """A function for label lines in graphs"""
    x = line.get_xdata()
    y = line.get_ydata()
    # place labels using numpy interpolation
    pos = [at_x + offset[0], np.interp(at_x, x, y) + offset[1]]
    ax.text(pos[0], pos[1], line.get_label(), size=fontsize, color=line.get_color(),
            ha="center", va="center", **kwargs)


# plot the results of conversion vs. T
xaxislabel = '$T\;\; \mathrm{[^\circ C]}$'
yaxislabel = "$\mathrm{CO_2 \; conversion \; [-]}$"
labels = [" $\mathrm{CC}$", "$\mathrm{CT}$", "$\mathrm{RWGS}$"]
label_x = [500,500,500,800]
offsets = [(10,0.07),(10,-0.07),(20,-0.07),(0,-0.07)]

lines = []

fig = plt.figure(figsize=(7, 4.5), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel, fontsize=axisfont)
ax.set_ylabel(yaxislabel, fontsize=axisfont)
ax.tick_params(labelsize=ticklabelfont)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis='both', which = 'both', direction='in', bottom = 'on', top= 'on', left = 'on', right='on')

a, = ax.plot(Trange-273, CC_co2_conversion, ls='-',  lw=1.8, color = 'black', label=labels[0])
b, = ax.plot(Trange-273, CT_co2_conversion, ls='-',  lw=1.8, color = 'black', label=labels[1])
c, = ax.plot(Trange-273, RWGS_co2_conversion, ls='--',  lw=1.8, color = 'grey', label=labels[2])
lines = [a,b,c]
for i, line in enumerate(lines):
    label_line(ax, line, at_x=label_x[i], offset=offsets[i], fontsize=labelfont)
ax.set_xlim(300,800)
ax.set_ylim(0,1)

plt.savefig(os.path.join('plots', 'figure_7b' + '.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_7b' + '.pdf'), dpi=150, bbox_inches='tight')
plt.show()


# Next plot the concentration profiles for a specific case
T = 773
p = 100000

kappamax_cc = maximum_kappa_CC(T, p)
kappamax_ct = maximum_kappa_CT(T, p)

# store results in arrays for plotting
pO2_flow1 = np.arange(0, 0.502, 0.002)
pO2_flow2_ct = np.linspace(0, kappamax_ct, 200)
pO2_flow2_cc = np.linspace(0,kappamax_cc, 200)


kappa = np.arange(0, 0.502, 0.002)
kappa_ct = np.linspace(0, kappamax_ct, 200)
kappa_cc = np.linspace(0, kappamax_cc, 200)

flow1 = ct.Solution('gri30.cti')
flow2 = ct.Solution('gri30.cti')

for i, k in enumerate(kappa):
    flow1.X = {'CO': 1, 'O2': 0.5-k}
    flow1.TP = T, p
    flow1.equilibrate('TP')
    pO2_flow1[i] = flow1.X[flow1.species_index('O2')]*p


for i, k in enumerate(kappa_ct):
    flow2.X = {'H2': 3, 'O2': k}
    flow2.TP = T, p
    flow2.equilibrate('TP')
    pO2_flow2_ct[i] = flow2.X[flow2.species_index('O2')]*p

for i, k in enumerate(kappa_cc):
    flow2.X = {'H2': 3, 'O2': k}
    flow2.TP = T, p
    flow2.equilibrate('TP')
    pO2_flow2_cc[i] = flow2.X[flow2.species_index('O2')]*p

axisfont = 16
ticklabelfont = 13
labelfont = 14

# PLOT pO2 vs. kappa
xaxislabel = '$ \kappa$ $\;\; \mathrm{[-]}$'
yaxislabel = '$p_\mathrm{O_2} \;\; \mathrm{[Pa]}$'
labels = ["$\mathrm{CO2}$",
          "$\mathrm{CC}$", "$\mathrm{CT}$"]

label_x = [0.1 , 0.1, 0.1]
offsets = [(0, 0), (0, 0), (0, 0)]

fig = plt.figure(figsize=(7, 4.5), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel, fontsize=axisfont)
ax.set_ylabel(yaxislabel, fontsize=axisfont)
ax.tick_params(labelsize=ticklabelfont)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis='both', which = 'both', direction='in', bottom = 'on', top= 'on', left = 'on', right='on')

ax.plot(kappa, pO2_flow1, ls='-',  lw=1.8, color = 'black', label=labels[0])
ax.plot(kappa_cc, pO2_flow2_cc, ls=':',  lw=1.8, color = 'black', label=labels[1])
ax.plot(kappamax_ct-kappa_ct, pO2_flow2_ct, ls='--',  lw=1.8, color = 'black', label=labels[2])


ax.text(0.08, 0.88, '$\mathrm{CO_2}$', horizontalalignment='center',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont, color = 'Black')
ax.text(0.08, 0.35, '$\mathrm{CC}$', horizontalalignment='center',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont, color = 'Black')
ax.text(0.08, 0.68, '$\mathrm{CT}$', horizontalalignment='center',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont, color = 'Black')

# ktatoal markers
ax.plot([kappamax_cc, kappamax_cc], [1E-40, 8.3E-25],  ls='-',  lw=1.0, color = 'grey')
ax.plot([kappamax_ct, kappamax_ct], [1E-40, 1E-27],  ls='-',  lw=1.0, color = 'grey')
#ax.plot([0.45, 0.5], [1E-28, 1E-30],  ls='-',  lw=1.0, color = 'black')
#labels
ax.text(0.91, 0.07, '$\mathrm{CT,} \; \kappa_\mathrm{max} \\rightarrow$', horizontalalignment='right',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont-1, color = 'Grey')
ax.text(0.49, 0.07, '$\mathrm{CC,} \; \kappa_\mathrm{max} \\rightarrow$', horizontalalignment='right',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont-1, color = 'Grey')

ax.text(0.75, 0.95, '$ \, T = 500\,\mathrm{^\circ C}$ \n $p =1 \, \mathrm{bar}$ \n $\omega = 3$', horizontalalignment='left',
       verticalalignment='top', transform=ax.transAxes, fontsize=labelfont-2, color = 'BLACK')

ax.set_xlim(0,0.5)
ax.set_ylim(1E-30,1E-20)
ax.set_yscale('log')


plt.savefig(os.path.join('plots', 'figure_7a' + '.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_7a' + '.pdf'), dpi=150, bbox_inches='tight')
plt.show()

