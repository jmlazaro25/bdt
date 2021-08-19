import os
import math
import ROOT as r
import numpy as np
from mods import ROOTmanager as manager
from mods import physTools, mipTracking
cellMap = np.loadtxt('mods/cellmodule.txt')
r.gSystem.Load('libFramework.so')

maxLayer = 20

# TreeModel to build here
branches_info = {
        f'energy_{llayer + 1}': {'rtype': float, 'default': 0.} \
                for llayer in range(maxLayer)
        }
branches_info['energy_tot'] = {'rtype': float, 'default': 0.}

def main():

    # Inputs and their trees and stuff
    pdict = manager.parse()
    batch_mode = pdict['batch']
    separate = pdict['separate']
    inlist = pdict['inlist']
    outlist = pdict['outlist']
    group_labels = pdict['groupls']
    startEvent = pdict['startEvent']
    maxEvents = pdict['maxEvents']
    # Should maybe put in parsing eventually and make event_process *arg

    # Construct tree processes
    procs = []
    for gl, group in zip(group_labels,inlist):
        procs.append( manager.TreeProcess(event_process, group,
                                          ID=gl, batch=batch_mode, pfreq=1000) )

    # Process jobs
    for proc in procs:

        # Move into appropriate scratch dir
        os.chdir(proc.tmp_dir)

        # Branches needed
        proc.ecalRecHits  = proc.addBranch('EcalHit', 'EcalRecHits_v12')

        # Tree/Files(s) to make
        print('\nRunning %s'%(proc.ID))

        proc.tfMakers = {'unsorted': None}
        for tfMaker in proc.tfMakers:
            proc.tfMakers[tfMaker] = manager.TreeMaker(group_labels[procs.index(proc)]+\
                                        '_{}.root'.format(tfMaker),\
                                        "EcalVeto",\
                                        branches_info,\
                                        outlist[procs.index(proc)]
                                        )

        # Gets executed at the end of run()
        proc.extrafs = [ proc.tfMakers[tfMaker].wq for tfMaker in proc.tfMakers ]

        # RUN
        proc.run(strEvent=startEvent, maxEvents=maxEvents)

    # Remove scratch directory if there is one
    if not batch_mode:     # Don't want to break other batch jobs when one finishes
        manager.rmScratch()

    print('\nDone!\n')

# Process an event
def event_process(self):

    # Initialize BDT input variables w/ defaults
    feats = next(iter(self.tfMakers.values())).resetFeats()

    # Add energies
    for hit in self.ecalRecHits:
        if hit.getEnergy() < 0: continue

        feats['energy_tot'] += hit.getEnergy() 
        for llayer in range(maxLayer):

            if physTools.ecal_layer(hit) <= llayer:
                feats[f'energy_{llayer + 1}'] += hit.getEnergy()

    self.tfMakers['unsorted'].fillEvent(feats)

if __name__ == "__main__":
    main()
