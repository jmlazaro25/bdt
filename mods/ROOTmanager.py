import os
import sys
import ROOT as r
import numpy as np
from mods import configuration as config

##################################################
# Constants
##################################################

# ROOT colors
colors = {
        'kBlue': 600,
        'kGreen': 417,
        'kMagenta': 616,
        'kOrange': 807, # +7
        'kBlack': 1,
        'kYellow': 400,
        'kViolet': 880,
        'kRed': 632,
        'kCyan': 432
        }

# For easier loops
color_list = [colors[key] for key in colors]

# ROOT 1D/2D line styles
lineStyles = {
        'kSolid': 1,
        'kDashed': 2,
        'kDotted': 3,
        'kDashDotted': 4
        }

# For easier loops
lineStyle_list = [i for i in range(1,11)]


##################################################
# Classes
##################################################

class TreeProcess:

    """ For analysing .root samples """

    def __init__(self,
            pConfig, 
            branches=None,
            event_process=None,
            startfs=None,
            endfs=None,
            funcs=None,
            strEvent=0,
            maxEvents=-1,
            pfreq=1000,
            color=1
            ):
        
        self.branches = branches
        self.event_process = event_process
        self.startfs = startfs
        self.endfs = endfs
        self.funcs = funcs
        self.strEvent = strEvent
        self.maxEvents = maxEvents
        self.pfreq = pfreq
        self.color = color
        self.batch =pConfig.batch
        self.infiles = pConfig.infiles
        self.tree = pConfig.tree; self.tree_name = pConfig.tree_name
        self.iD = pConfig.iD
        self.cwd = os.getcwd()

        print('\nPreparing: {}'.format(self.iD))

        # Build tree amd move operations to a scratch directory
        # if providing group_files instead of a tree

        self.mvd = False
        if self.tree == None:
            self.mvd = True

            # Give warning if config passed no valid paths
            if self.infiles == []: sys.exit('No valid paths')

            # Create the scratch directory if it doesn't already exist
            scratch_dir = self.cwd + '/scratch'
            print( 'Using scratch path %s' % scratch_dir )
            if not os.path.exists(scratch_dir):
                os.makedirs(scratch_dir)

            # Get tmp num
            num=0; check = True
            while check:
                if os.path.exists( scratch_dir+'/tmp_'+str(num) ):
                    num += 1
                else:
                    check = False 

            # Create and mv into tmp directory that can be used to copy files into
            if self.batch:
                self.tmp_dir='%s/%s' % (scratch_dir, os.environ['LSB_JOBID'])
            else:
                self.tmp_dir = '%s/%s' % (scratch_dir, 'tmp_'+str(num))
            if not os.path.exists(self.tmp_dir):
                print( 'Creating tmp directory %s' % self.tmp_dir )
            os.makedirs(self.tmp_dir)
            os.chdir(self.tmp_dir)

            # Copy input files to the tmp directory
            print( 'Copying input files into tmp directory' )
            for rfilename in self.infiles:
                os.system("cp %s ." % rfilename )
            os.system("ls .")
    
            # Just get the file names without the full path
            tmpfiles = [f.split('/')[-1] for f in self.infiles]
    
            # Load'em
            if self.tree_name != None:
                self.tree = load(tmpfiles, self.tree_name)
            else:
                self.tree = load(tmpfiles)

            # Add input branaches as attributes
            if self.branches != None:
                for btup in self.branches:
                    setattr(self, btup[0], self.addBranch(btup[1], btup[0]))    

            # Move back to cwd in case running multiple procs
            os.chdir(self.cwd)

    def addBranch(self, ldmx_class, branch_name):

        """  Add a new branch to read from """

        if self.tree == None:
            sys.exit('Set tree')

        if ldmx_class == 'EventHeader': branch = r.ldmx.EventHeader()
        elif ldmx_class == 'EcalVetoResult': branch = r.ldmx.EcalVetoResult()
        elif ldmx_class == 'HcalVetoResult': branch = r.ldmx.HcalVetoResult()
        elif ldmx_class == 'TriggerResult': branch = r.ldmx.TriggerResult()
        elif ldmx_class == 'SimParticle': branch = r.map(int, 'ldmx::'+ldmx_class)()
        else: branch = r.std.vector('ldmx::'+ldmx_class)()

        self.tree.SetBranchAddress(branch_name,r.AddressOf(branch))

        return branch
 
    def run(self, strEvent=0, maxEvents=-1, pfreq=1000):
  
        """ Run event_process on each event + any start/end funcitons """

        # Move into appropriate scratch dir
        os.chdir(self.tmp_dir)

        # Init
        if strEvent != 0: self.strEvent = strEvent
        if maxEvents != -1: self.maxEvents = maxEvents
        if self.maxEvents == -1 or self.strEvent + self.maxEvents >\
                self.tree.GetEntries():
            self.maxEvents = self.tree.GetEntries() - self.strEvent
        maxEvent = self.strEvent + self.maxEvents
        if pfreq != 1000: self.pfreq = pfreq

        # Execute any opening function(s)
        # (might impliment *args, **kwargs later)
        if self.startfs != None:
            for startf in self.startfs:
                startf()

        # Loop over events
        print( '\nRunning: {}\n'.format(self.iD) )
        self.event_count = self.strEvent
        while self.event_count < maxEvent:
            self.tree.GetEntry(self.event_count)
            if self.event_count%self.pfreq == 0:
                print('Processing Event: %s'%(self.event_count))
            self.event_process()
            self.event_count += 1

        # Move back to cwd in case running multiple procs
        os.chdir(self.cwd)

        # Execute any closing function(s)
        # (might impliment *args, **kwargs later)
        if self.endfs != None:
            for endf in self.endfs:
                endf()

        # Remove tmp directory if created in move
        """ Giving Exited message in batch even when jobs run fine
        if self.mvd:
            print( 'Removing tmp directory %s' % self.tmp_dir )
            os.system('rm -rf %s' % self.tmp_dir)
        """

class TreeMaker:

    """ To write a tree in an analysis process """

    import numpy as np

    def __init__(self, tConfig):

        self.outfile = tConfig.outfile
        self.tree_name = tConfig.tree_name
        self.branches_info = tConfig.branches_info
        self.branches = {}
        self.outdir = tConfig.outdir

        # Create output directory
        if self.outdir != '':
            if not os.path.exists(self.outdir):
                print( 'Creating %s' % (self.outdir) )
                os.makedirs(self.outdir)

        # Create output file and tree
        self.tfout = r.TFile(self.outdir+'/'+self.outfile,"RECREATE")
        self.tree = r.TTree(self.tree_name, self.tree_name)

        # Set up new tree branches if given branches_info
        if self.branches_info != {}:
            for branch_name in self.branches_info:
                self.addBranch(
                                self.branches_info[branch_name]['rtype'],
                                branch_name
                                )

    def addBranch(self, rtype, branch_name):

        """ Add a new branch to write to """

        # Basic types
        if str(rtype) == "<type 'float'>" or str(rtype) == "<class 'float'>":
            rtypeL = '/D'
        elif str(rtype) == "<type 'int'>" or str(rtype) == "<class 'int'>":
            rtypeL = '/I'
        elif str(rtype) == "<type 'bool'>" or str(rtype) == "<class 'bool'>":
            rtypeL = '/O'

        if type(rtype) != str:
            self.branches[branch_name] = np.zeros(1, dtype=rtype)
            self.tree.Branch(
                    branch_name,
                    self.branches[branch_name],
                    branch_name + rtypeL
                    )

        # Others
        elif type(rtype) == str and rtype[:7] == 'vectorOF': # Untested
            self.branches[branch_name] = r.std.vector( rtype[7:] )()
            self.tree.Branch(branch_name, self.branches[branch_name])

    def resetFeats(self):

        """
        Reset variables to defaults for new event
        Return because feats['feat'] looks nicer than
        self.tfMaker.feats['feat']
        """

        feats = {}
        for branch_name in self.branches_info:
            feats[branch_name] = self.branches_info[branch_name]['default']

        return feats

    def fillEvent(self, feats):

        """ Fill the tree with new feature values """

        for feat in feats:
            try:
                self.branches[feat][0] = feats[feat]
            except:
                sys.exit( str(feat) + ': ' + str(feats[feat]) )
        self.tree.Fill()

    def wq(self):

        """ Save the tree and close the file """

        self.tfout.Write(self.tree_name)
        self.tfout.Close()

class Histogram:

    """ Just to hold histogram-related stuff and make other py code nicer """

    def __init__(self, hist, title='', xlabel='x', ylabel='y',\
            color=1, lineStyle=1, fillStyle=1):
        self.hist = hist
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.color = color
        self.lineStyle = lineStyle
        self.fillStyle = fillStyle

##################################################
# Functions
##################################################

def parse(nolist = False):

    from glob import glob
    from argparse import ArgumentParser

    # Arguments
    parser = ArgumentParser()
    parser.add_argument('action', action='store',
            default=None,
            help='BDT actor string (trees, train, or eval)')
    parser.add_argument('config', action='store',
            help='BDT configuration file')
    parser.add_argument('-p', nargs='+', action='store', dest='mprocs',
            default=None,
            help='Process names to be run from config e.g. pn, 1.0 (for sig)')
    parser.add_argument('--batch', action='store_true', dest='batch',
            default=False,
            help='Run in batch mode [Default: False]')
    parser.add_argument('--sep', action='store_true', dest='separate',
            default = False,
            help='separate events into different files [Default: False]')
    parser.add_argument('-i', nargs='+', action='store', dest='infiles',
            default=[],
            help='input file(s)')
    parser.add_argument('--indirs', nargs='+', action='store', dest='indirs',
            default=[],
            help='Director(y/ies) of input files')
    parser.add_argument('-g','-groupls', nargs='+', action='store',
            dest='group_labels',
            default='',
            help='Human readable sample labels e.g. for legends')
    parser.add_argument('-o','--out', nargs='+', action='store', dest='out',
            default=[],
            help='output files or director(y/ies) of output files')
            # if inputting directories, it's best to make a system
            # for naming files in main() of main script 
    parser.add_argument('--notlist', action='store_true', dest='nolist',
            help="return things without lists")
    parser.add_argument('-s','--start', type=int, action='store',
            dest='startEvent',
            default=0, help='event to start at')
    parser.add_argument('-m', type=int, action='store', dest='maxEvents',
            default=-1, help='max events to run over for EACH group')
    args = parser.parse_args()

    # Input
    if args.infiles != []:
        inlist = [[f] for f in args.infiles] # Makes general loading easier
        if nolist or args.nolist == True:
            inlist = inlist[0]
    elif args.indirs != []:
        inlist = [glob(indir + '/*.root') for indir in args.indirs]
        if nolist or args.nolist == True:
            inlist = inlist[0]
    else: inlist = []

    # Output
    if args.out != []:
        outlist = args.out
        if nolist or args.nolist == True:
            outlist = outlist[0]
    else: outlist = []
    
    pdict = {
            'action': args.action,
            'config': args.config,
            'mprocs': args.mprocs,
            'batch': args.batch,
            'separate': args.separate,
            'inlist': inlist,
            'groupls': args.group_labels,
            'outlist': outlist,
            'startEvent': args.startEvent,
            'maxEvents': args.maxEvents
            }

    return pdict

def load(group,treeName='LDMX_Events'):

    """ Load a group of files into a readable tree """

    tree = r.TChain(treeName)
    for f in group:
        tree.Add(f)

    return tree

def rmScratch():

    """ Remove scratch dir """

    if os.path.exists('./scratch'):
        print( '\nRemoving scratch directory' )
        os.system('rm -rf ./scratch')
