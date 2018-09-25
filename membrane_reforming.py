import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os

def maximum_kappa_CC(T, p):
    """Determine the maximum exchange of oxygen kappa_max, for cocurrent flows of CO2 and CH4 with equal flow rates, omega=1.
    Returns, CO2_conversion =  kappa_max, CH4_conversion, and mole fractions of main species in the CH4 stream."""
    # define the two flows as canter solution objects using cantera's gri30 database
    flow1 = ct.Solution('gri30.cti')
    flow2 = ct.Solution('gri30.cti')
    # complete conversion of CO2 to CO gives Kappa = 0.5, kappa in range 0-0.5
    kappa = np.arange(0, 0.50, 0.002)
    for k in kappa:
        # iterates through kappa until the CO2 and CH4 oxygen partial pressures meet in CC configuration
        # equal flows of CH4 and CO2, i.e. omega=1
        flow1.X = {'CO': 1, 'O2': 0.5 - k}
        flow2.X = {'CH4': 1, 'O2': k}
        flow1.TP = T, p
        flow2.TP = T, p
        flow1.equilibrate('TP')
        flow2.equilibrate('TP')
        pO2_flow1 = flow1.X[flow1.species_index('O2')]*p
        pO2_flow2 = flow2.X[flow2.species_index('O2')]*p
        if (pO2_flow2 > pO2_flow1 or k > 0.497):
            kappa_max =  k
            # get all mole fractions
            xCH4 = flow2.X[flow2.species_index('CH4')]
            xH2O = flow2.X[flow2.species_index('H2O')]
            xCO2 = flow2.X[flow2.species_index('CO2')]
            xCO = flow2.X[flow2.species_index('CO')]
            xH2 = flow2.X[flow2.species_index('H2')]
            CH4_conversion = 1-xCH4/(xCH4+xCO2+xCO)
            mole_fractions_flow2 = [xCH4, xCO2, xH2O, xCO, xH2]
            break
    return(kappa_max, CH4_conversion, mole_fractions_flow2)

def maximum_kappa_CT(T, p):
    """Determine the maximum conversion of oxygen for countercurrent flows of CO2 and CH4 with equal flow rates, omega=1.
    Returns, CO2_conversion =  Kappa_max, CH4_conversion, and mole fractions of main species in the CH4 stream."""
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
        flow2.X = {'CH4': 1, 'O2': k}
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
            kappa_max =  k
            # get all mole fractions
            xCH4 = flow2.X[flow2.species_index('CH4')]
            xH2O = flow2.X[flow2.species_index('H2O')]
            xCO2 = flow2.X[flow2.species_index('CO2')]
            xCO = flow2.X[flow2.species_index('CO')]
            xH2 = flow2.X[flow2.species_index('H2')]
            CH4_conversion = 1-xCH4/(xCH4+xCO2+xCO)
            mole_fractions_flow2 = [xCH4, xCO2, xH2O, xCO, xH2]
            break
    return(kappa_max, CH4_conversion, mole_fractions_flow2)



# Lets look at the thermodynamics as a function of temperature in the range,
Trange = np.arange(673, 1383, 10)
p=100000

# store results in arrays
CC_co2_conversion = np.zeros(len(Trange))
CT_co2_conversion = np.zeros(len(Trange))
CC_ch4_conversion = np.zeros(len(Trange))
CT_ch4_conversion = np.zeros(len(Trange))
CC_molefractions = np.zeros((len(Trange), 5))
CT_molefractions = np.zeros((len(Trange), 5))

for i, T in enumerate(Trange):
    # Caculate maximum kappa  for CC case
    cc_results = maximum_kappa_CC(T, p)
    CC_co2_conversion[i] = 2*cc_results[0]
    CC_ch4_conversion[i] = cc_results[1]
    CC_molefractions[i] = cc_results[2]
    # Caculate maximum kappa for CT case
    ct_results = maximum_kappa_CT(T, p)
    CT_co2_conversion[i] = 2*ct_results[0]
    CT_ch4_conversion[i] = ct_results[1]
    CT_molefractions[i] = ct_results[2]
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


# Plot the results of CO2 and CH4 conversion vs. T
# Figure 6b in the manuscript
xaxislabel = '$T\;\; \mathrm{[^\circ C]}$'
yaxislabel = "$\mathrm{Conversion \; extent\;\; [-]}$"
labels = [" ", "$\mathrm{\\leftarrow  CO_2 \\rightarrow}$", " ", "$\mathrm{CH_4}$"]


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
b, = ax.plot(Trange-273, CT_co2_conversion, ls='--',  lw=1.8, color = 'black', label=labels[1])
c, = ax.plot(Trange-273, CC_ch4_conversion, ls='-',  lw=1.8, color = 'black', label=labels[0])
d, = ax.plot(Trange-273, CT_ch4_conversion, ls='--',  lw=1.8, color = 'black', label=labels[3])

lines = [a,b,c,d]
label_x = [500,500,720,800]
offsets = [(10,0.07),(35,-0.07),(20,-0.07),(0,-0.07)]
for i, line in enumerate(lines):
    label_line(ax, line, at_x=label_x[i], offset=offsets[i], fontsize=labelfont)
ax.set_xlim(400,900)
ax.set_ylim(0,1)

plt.savefig(os.path.join('plots', 'figure_6b' + '.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_6b' + '.pdf'), dpi=150, bbox_inches='tight')
plt.show()


# Plot mole fractions for the CH4 stream in the CC case, not in the manuscript
lines =[]
xaxislabel = '$T\;\; \mathrm{[^\circ C]}$'
yaxislabel = '$x \;\; [-]}$'
labels = ["$\mathrm{CH_4}$","$\mathrm{CO_2}$","$\mathrm{H_2O}$",
          "$\mathrm{CO}$", "$\mathrm{H_2}$"]

label_x = [470,500, 550, 900, 900]
offsets = [(10,0.06),(0,0.04),(0,-0.04),(0,0.04),(0,0.04)]
fig = plt.figure(figsize=(7, 4.5), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel, fontsize=axisfont)
ax.set_ylabel(yaxislabel, fontsize=axisfont)
ax.tick_params(labelsize=ticklabelfont)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis='both', which = 'both', direction='in', bottom = 'on', top= 'on', left = 'on', right='on')
for i, col in enumerate(CC_molefractions.T):
    a, = ax.plot(Trange-273, col, ls='-',  lw=1.8, color = 'black', label=labels[i])
    lines.append(a)
for i, line in enumerate(lines):
    label_line(ax, line, at_x=label_x[i], offset=offsets[i], fontsize=labelfont)
ax.set_xlim(400,1100)
ax.set_ylim(0,0.8)

plt.savefig(os.path.join('plots', 'CC_ch4_conversion_x' + '.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'CC_ch4_conversion_x' + '.pdf'), dpi=150, bbox_inches='tight')
plt.show()

# Plot mole fractions for the CH4 stream in the CT case, not in the manuscript
lines =[]
xaxislabel = '$T\;\; \mathrm{[^\circ C]}$'
yaxislabel = '$x \;\; [-]}$'
labels = ["$\mathrm{CH_4}$","$\mathrm{CO_2}$","$\mathrm{H_2O}$",
          "$\mathrm{CO}$", "$\mathrm{H_2}$"]

label_x = [450,525, 600, 900, 900]
offsets = [(20,0.06),(0,0.04),(0,-0.045),(0,0.04),(0,0.04)]
fig = plt.figure(figsize=(7, 4.5), facecolor='white')
ax = plt.subplot()
ax.set_xlabel(xaxislabel, fontsize=axisfont)
ax.set_ylabel(yaxislabel, fontsize=axisfont)
ax.tick_params(labelsize=ticklabelfont)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis='both', which = 'both', direction='in', bottom = 'on', top= 'on', left = 'on', right='on')
for i, col in enumerate(CT_molefractions.T):
    a, = ax.plot(Trange-273, col, ls='-',  lw=1.8, color = 'black', label=labels[i])
    lines.append(a)
for i, line in enumerate(lines):
    label_line(ax, line, at_x=label_x[i], offset=offsets[i], fontsize=labelfont)
ax.set_xlim(400,1000)
ax.set_ylim(0,0.8)

plt.savefig(os.path.join('plots', 'CT_ch4_conversion_x' + '.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'CT_ch4_conversion_x' + '.pdf'), dpi=150, bbox_inches='tight')
plt.show()

plt.show()



# Show the partial pressures at the limiting case for a given temperature.
# figure 6a in Manuscript

T = 873
p = 100000

# In this case countercurrent allows for complete CO2 conversion
cc_results = maximum_kappa_CC(T, p)
kappa_max_cc = cc_results[0]

pO2_flow1 = np.arange(0, 0.50, 0.001)
pO2_flow2_CT = np.arange(0, 0.50, 0.001)
pO2_flow2_CC = np.linspace(0, kappa_max_cc, 200)


kappa = np.arange(0, 0.50, 0.001)
kappa_ct = np.arange(0.5, 0, -0.001)
kappa_cc = np.linspace(0,kappa_max_cc, 200)


# Define our two flows as cantera solutions
flow1 = ct.Solution('gri30.cti')
flow2 = ct.Solution('gri30.cti')
for i, k in enumerate(kappa):
    flow1.X = {'CO': 1, 'O2': 0.5-k}
    flow2.X = {'CH4': 1, 'O2': k}
    flow1.TP = T, p
    flow2.TP = T, p
    flow1.equilibrate('TP')
    flow2.equilibrate('TP')
    pO2_flow1[i] = flow1.X[flow1.species_index('O2')]*p
    pO2_flow2_CT[i] = flow2.X[flow2.species_index('O2')]*p


for i, k in enumerate(kappa_cc):
    flow2.X = {'CH4': 1, 'O2': k}
    flow2.TP = T, p
    flow2.equilibrate('TP')
    pO2_flow2_CC[i] = flow2.X[flow2.species_index('O2')]*p

axisfont = 16
ticklabelfont = 13
labelfont = 14

# PLOT pO2 vs. kappa
xaxislabel = '$ \kappa$ $\;\; \mathrm{[-]}$'
yaxislabel = '$p_\mathrm{O_2} \;\; \mathrm{[Pa]}$'

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
a, = ax.plot(kappa, pO2_flow1, ls='-',  lw=1.8, color = 'black')
b, = ax.plot(kappa_cc, pO2_flow2_CC, ls=':',  lw=1.8, color = 'black')
c, = ax.plot(kappa_ct, pO2_flow2_CT, ls='--',  lw=1.8, color = 'black')

# labels
ax.text(0.1, 0.9, '$\mathrm{CO_2}$', horizontalalignment='center',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont, color = 'Black')
ax.text(0.1, 0.35, '$\mathrm{CC}$', horizontalalignment='center',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont, color = 'Black')
ax.text(0.1, 0.7, '$\mathrm{CT}$', horizontalalignment='center',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont, color = 'Black')

# kappa_max markers
ax.plot([kappa_max_cc, kappa_max_cc], [1E-40, 3.46E-21],  ls='-',  lw=1.0, color = 'grey')
#ax.plot([0.45, 0.5], [1E-28, 1E-30],  ls='-',  lw=1.0, color = 'black')
#labels
ax.text(0.99, 0.07, '$\mathrm{CT,} \; \kappa_\mathrm{max} \\rightarrow$', horizontalalignment='right',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont-1, color = 'GREY')
ax.text(0.69, 0.07, '$\mathrm{CC,} \; \kappa_\mathrm{max} \\rightarrow$', horizontalalignment='right',
       verticalalignment='center', transform=ax.transAxes, fontsize=labelfont-1, color = 'GREY')
ax.text(0.75, 0.95, '$ \, T = 600\,\mathrm{^\circ C}$ \n $p =1 \, \mathrm{bar}$ \n $\omega = 1$',
        horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=labelfont-2,
        color = 'BLACK')

ax.set_xlim(0,0.5)
ax.set_ylim(1E-30,1E-15)
ax.set_yscale('log')

plt.savefig(os.path.join('plots', 'figure_6a' + '.png'), dpi=150, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_6a' + '.pdf'), dpi=150, bbox_inches='tight')
plt.show()