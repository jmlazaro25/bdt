import numpy as np
from os import getcwd
import matplotlib.pyplot as plt
import mods.ROOTmanager as manager

# First assumed to be bkg
cwd = getcwd()
procs = {
        'Inclusive': cwd+'/testout/inc_unsorted.root',
        'ECal PN': cwd+'/testout/pn_unsorted.root',
        r'Late $\gamma$ conv': cwd+'/testout/conv_unsorted.root',
        r"$m_{A'}$= 5 MeV": cwd+'/testout/0.005_unsorted.root',
        r"$m_{A'}$= 10 MeV": cwd+'/testout/0.01_unsorted.root',
        r"$m_{A'}$= 50 MeV": cwd+'/testout/0.05_unsorted.root'
        }
bkg_name = next(iter( procs.keys() ))

# Energies set up
marker_en = 3100
marker_ind = (marker_en - 100)//100
#marker_ind = 59 - 5 # 54 => 5.5 GeV
minenergy = 3_000
maxenergy = 10_000 # MeV
energy_space = 100 # MeV
energies = np.arange(
                energy_space, # Don't start at 0
                maxenergy + energy_space, # Include the last one
                energy_space
            )

effs = { proc: np.zeros(len(energies)) for proc in procs}
# Get efficiencies
for proc, file in procs.items():

    tree = manager.load( [file] , 'EcalVeto' )

    # Loop over tree for each en
    for en_ind, en in enumerate( energies ):
        
        # Calc frac of events passing this en cut
        npass_evs = 0
        for ev in tree:
            if ev.energy_20 < en: npass_evs += 1
        effs[proc][en_ind] = npass_evs / tree.GetEntries()
    print(effs[proc][marker_ind])


# Print effs so we know points and fractions
print(f'energies:')
print(energies)
print(f'{bkg_name}:')
print(effs[ bkg_name ])
for sig in tuple( procs.keys() )[1:]:
    print(f'{sig}:')
    print(effs[sig])
    print(effs[sig]/effs[bkg_name])

# Plot
for sig in tuple( procs.keys() )[1:]:
    plt.plot(effs[ bkg_name ], effs[sig], label=sig)
    plt.scatter(
            effs[bkg_name][marker_ind],
                effs[sig][marker_ind],
                marker='*',
                s=100,
                alpha=1,
                color='red'
                )

plt.xlim(0,1)
plt.ylim(0,1)
plt.margins(x=0)
plt.margins(y=0)
plt.xlabel(r'$\epsilon$'+f'({bkg_name})')
plt.ylabel(r'$\epsilon$(vissible signal)')
plt.legend()

plt.show()
