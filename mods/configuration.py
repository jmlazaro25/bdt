import sys
from glob import glob
from mods import feats
from configparser import ConfigParser

# NOTE: Add handling for batch outfile names in tConf

def print_dict(d, prefix='', skip=set()):
    
    """ Print dictionary neatly (esp. useful for configurations) """

    if prefix != '':
        print(prefix)

    for k,v in d.items():
        if k not in skip:
            print('{}: {}'.format(k,v))

def tree_building_info(configFile):

    """ Parse branches for tree making action """

    config = ConfigParser()
    config.read(configFile)
    
    tracker = config.get('trees', 'tracker')
    ecal    = config.get('trees', 'ecal')
    hcal    = config.get('trees', 'hcal')
    analysis = config.get('trees', 'analysis')

    # Used to name trees
    dict_names = tuple(
            [ d_name for d_name in (tracker, ecal, hcal, analysis) \
                    if d_name != 'None']
            )

    # Build branches_info and a list of the funcs needed to calculations
    branches = {}
    fs_dict = {}
    for d_name in dict_names:

        # Get dictionaries from feats
        for dn in d_name.split('AND'):
            if 'trees_info_{}'.format(dn) in feats.__dict__:
                branches.update( getattr(feats, 'trees_info_{}'.format(dn)) )
                fs_dict.update( getattr(feats, '{}_funcs'.format(dn)) )
            else: sys.exit('\nInvalid dictionary name: {}'.format(dn))

    return branches, fs_dict, '_'.join(dict_names) + '_Veto'

def feature_info(configFile):

    """ Parse branches used in bdt training and evaluation """

    config = ConfigParser()
    config.read(configFile)
    
    tracker = config.get('setup', 'tracker')
    ecal    = config.get('setup', 'ecal')
    hcal    = config.get('setup', 'hcal')

    # Used to name eval tree and discValue
    d_names = tuple( [l_name for l_name \
                        in (tracker, ecal, hcal) if l_name != ''] )

    # Build feature dictionary
    f_branches = {}
    for d_name in d_names:
        if 'feats_{}'.format(d_name) in feats.__dict__:
            f_branches.update( getattr(feats, 'feats_{}'.format(d_name)) )
        else: sys.exit('\nInvalid feature dictionary: {}'.format(d_name))

    return f_branches, '_'.join(d_names)

def parse_bdt_config(action_str, configFile, clargs = {}):

    """ Main parsing function outputing config objects below """

    config = ConfigParser()
    config.read(configFile)

    procs   = config.get('setup', 'processes').replace(' ','').split(',')
    labels  = config.get('setup', 'labels').replace(' ','').split(',')
    branches_info, calc_funcs, tree_nem = tree_building_info(configFile)

    # Overwrite processes list from config if specified in command line
    proc_configs = []
    if clargs['mprocs'] != None:

        procs = clargs['mprocs']

        # Give overriding message
        print('\nOnly running over the following processes:')
        for proc in procs: print(proc)

    # Action dependent items
    if action_str == 'trees':
        
        train_or_test = config.get('trees', 'train_or_test')
        if train_or_test != '':
            print('\nOnly running over {}ing set!!!'.format(train_or_test))
        sets = {'test', 'train'} if train_or_test == '' else {train_or_test} 

        for proc in procs:
            for st in sets:

                pConf = ProcessConfig(
                        infiles = glob(
                            config.get( 'trees', '{}_indir'.format(st) )+\
                                    '/{}/*.root'.format(proc)
                                ),
                        tree_name = 'LDMX_Events',
                        iD = '{}_{}'.format(proc, st),
                        strEvent = clargs['startEvent'],
                        maxEvents = clargs['maxEvents'],
                        pfreq = 1000,
                        batch = clargs['batch']
                        )
                
                tConf = TreeConfig(
                        outfile = '{}_{}.root'.format(proc, st),
                        tree_name = tree_nem,
                        branches_info = branches_info,
                        funcs = calc_funcs,
                        outdir = \
                            config.get('trees', '{}_outdir'.format(st))+\
                                    '/{}'.format(proc)
                        )

                proc_configs.append( BdtConfig(pConf, tConf) )

        # Only print branches_info once
        #proc_configs[0].tConfig.print_branches()

        return proc_configs            

    # Won't need this anymore
    del calc_funcs

    if action_str == 'train':
        
        # Overwrite branches info
        branches_info = feature_info(configFile)[0]

        return ({
                'tree_name': tree_nem,
                'branches': (*branches_info,), # Only need names for training
                'bkg_indir': config.get('train', 'bkg_indir'),
                'sig_indir': config.get('train', 'sig_indir'),
                'outdir': config.get('train', 'outdir'),
                'seed': config.get('train', 'seed'),
                'eta': config.get('train', 'eta'),
                'tree_number': config.get('train', 'tree_number'),
                'tree_depth': config.get('train', 'tree_depth'),
                'batch': clargs['batch']
                },) # Tuple wrapper for process init loop in mainframe

    if action_str == 'eval':

        # Overwrite branches_info and Add discValue
        branches_info, name = feature_info(configFile)
        branches_info[ 'discValue_'+name ] = { 'rtype': float, 'default': 0.5 }

        for proc in procs:

                pConf = ProcessConfig(
                        infiles = glob(
                            config.get( 'eval', 'indir' ) \
                                    + '/{}/*.root'.format(proc)
                                ),
                        tree_name = tree_nem,
                        iD = proc + '_eval',
                        strEvent = clargs['startEvent'],
                        maxEvents = clargs['maxEvents'],
                        pfreq = 1000,
                        batch = clargs['batch'],
                        bdt = config.get('eval', 'bdt')
                        )

                tConf = TreeConfig(
                        outfile = '{}_eval.root'.format(proc),
                        tree_name = name + '_Veto',
                        branches_info = branches_info,
                        outdir =  config.get('eval', 'outdir')
                        )

                proc_configs.append( BdtConfig(pConf, tConf) )

        # Only print branches_info once
        #proc_configs[0].tConfig.print_branches()

        return proc_configs 

def parse_batch_config(clargs):

    """ Build single BdtConfig object for input (ideally from batch) """

    config = ConfigParser()
    config.read(clargs['config'])

    infile = clargs['inlist'][0][0]
    fname = infile.split('/')[-1].split('.root')[0]
    branches_info, calc_funcs, tree_nem = tree_building_info(clargs['config'])

    print( '\nUsing {} on {}'.format(clargs['action'], infile) )

    # Action dependent items
    if clargs['action'] == 'trees':

        pConf = ProcessConfig(
                infiles = [infile],
                tree_name = 'LDMX_Events',
                iD = '{}_tree'.format(fname),
                strEvent = clargs['startEvent'],
                maxEvents = clargs['maxEvents'],
                pfreq = 1000,
                batch = True
                )

        tConf = TreeConfig(
                outfile = fname + '_tree.root',
                tree_name = tree_nem,
                branches_info = branches_info,
                funcs = calc_funcs,
                outdir = clargs['outlist'][0]
                )

        return BdtConfig(pConf, tConf, printbrs=True )

    # Won't need this anymore
    del calc_funcs

    if clargs['action'] == 'eval':

        # Overwrite branches_info and add discValue
        branches_info, name = feature_info(clargs['config'])
        branches_info[ 'discValue_'+name ] = { 'rtype': float, 'default': 0.5 }

        pConf = ProcessConfig(
                infiles = [infile],
                tree_name = tree_nem,
                iD = '{}_eval'.format(fname),
                strEvent = clargs['startEvent'],
                maxEvents = clargs['maxEvents'],
                pfreq = 1000,
                batch = clargs['batch'],
                bdt = config.get('eval', 'bdt')
                )

        tConf = TreeConfig(
                outfile = '{}_eval.root'.format(fname),
                tree_name = name + 'Veto',
                branches_info = branches_info,
                outdir = clargs['outlist'][0] 
                )

        return BdtConfig(pConf, tConf, printbrs=True )

class BdtConfig:

    """
    Has all ProcessConfigs and TreeConfigs for
    neatly an uniformly passing to all BDT actors.
    """

    def __init__(self, pConfig=None, tConfig=None, printbrs=False):
        self.pConfig = pConfig
        self.tConfig = tConfig

        self.print(printbrs)

    def print(self, printbs):

        """ Print parameters neatly """

        print() # Some space for clearity

        if self.pConfig != None: self.pConfig.print()

        if self.tConfig != None: self.tConfig.print( print_branches=printbs )

class ProcessConfig:

    """
    A container for TreeProcess confugrations
    Might maka a method to parse config file here
    """

    def __init__(
            self,
            event_process = None,
            infiles = [],
            tree=None, tree_name = None, # If using both, you're doing it wrong
            startfs = None,
            endfs = None,
            iD = '',
            strEvent = 0,
            maxEvents = -1,
            pfreq = 1000,
            batch = False,
            color = 1,
            bdt = None # For eval processes
            ):

        self.event_process = event_process
        self.infiles = infiles
        self.tree = tree; self.tree_name = tree_name
        self.startfs = startfs
        self.endfs = endfs
        self.iD = iD
        self.strEvent = strEvent
        self.maxEvents = maxEvents
        self.pfreq = pfreq
        self.batch = batch
        self.color = color
        self.bdt = bdt

    def print(self):

        """ Print parameters neatly """
       
        print_dict(self.__dict__, '\n{} configuration:'.format(self.iD))

class TreeConfig:

    """ A container for TreeMaker confugrations """
    
    def __init__(
            self,
            outfile = 'outfile.root',
            tree_name = 'tree_name',
            branches_info = {},
            funcs = {}, # For building event_process in tree making
            outdir = 'outdir'
            ):

        self.outfile = outfile
        self.tree_name = tree_name
        self.branches_info = branches_info
        self.funcs = funcs
        self.outdir = outdir

    def print(self, print_branches=False):

        """ Make branches_info (and other things) human readable """

        print_dict(self.__dict__,
                prefix = '\nTreeMaker config:',
                skip = {'branches_info'}
                )

        if print_branches: self.print_branches()

    def print_branches(self):

        """ Extra care for branches_info + it takes up a lot of space """

        print_dict( self.branches_info, prefix='\nFeatures:')
