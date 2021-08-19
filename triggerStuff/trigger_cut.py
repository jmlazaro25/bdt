import ROOT as r
from glob import glob
from os import path, makedirs
from argparse import ArgumentParser
r.gSystem.Load('libFramework.so')

def main():

    # Parse
    parser = ArgumentParser()
    parser.add_argument('-i', dest='infile')
    parser.add_argument('-o', dest='outdir')
    args = parser.parse_args()

    # Skim
    if not path.exists(args.outdir): makedirs(args.outdir)
    trigger_file(args.infile, args.outdir)

def trigger_file(rfile, outdir):

    # Tha two constants we need
    cut = 3_000 # MeV
    maxLayer = 20

    # Load original
    ogTree = load( rfile )
    ecalRecHits = addBranch(ogTree, 'EcalHit', 'EcalRecHits_v12')

    # Tskim file and tree
    outfile = outdir + '/' + rfile.split('/')[-1][:-5] + '_tskim.root'
    tskimF = r.TFile( outfile, 'RECREATE')
    tskimT = ogTree.CloneTree(0)

    # Loop over events
    for entry in range(cut):
        ogTree.GetEntry(entry)

        # Check trigger pass
        energy_tot = 0
        for hit in ecalRecHits:
            if hit.getEnergy() < 0: continue

            if layer(hit) < maxLayer:
                energy_tot += hit.getEnergy()

        # If passes, fill new tree
        if energy_tot < cut:
            tskimT.Fill()

    # Wrap up
    tskimF.cd()
    tskimT.Write()
    tskimF.Close()

def load(fil,treeName='LDMX_Events'):

    # Load ROOT tree

    twee = r.TChain(treeName)
    twee.Add(fil)

    return twee

def addBranch(tree, ldmx_class, branch_name):

    """ Add a new branch to read from """

    if ldmx_class == 'EventHeader': branch = r.ldmx.EventHeader()
    elif ldmx_class == 'EcalVetoResult': branch = r.ldmx.EcalVetoResult()
    elif ldmx_class == 'HcalVetoResult': branch = r.ldmx.HcalVetoResult()
    elif ldmx_class == 'TriggerResult': branch = r.ldmx.TriggerResult()
    elif ldmx_class == 'SimParticle': branch = r.map(int, 'ldmx::'+ldmx_class)()
    else: branch = r.std.vector('ldmx::'+ldmx_class)()

    tree.SetBranchAddress(branch_name,r.AddressOf(branch))

    return branch

def layer(hit):

    """ Get layerID from ecal hit """

    LAYER_MASK = 0x3F  # space for up to 64 layers
    LAYER_SHIFT = 17

    return (hit.getID() >> LAYER_SHIFT) & LAYER_MASK

if __name__ == '__main__': main()
