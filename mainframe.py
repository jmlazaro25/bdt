import sys
import ROOT
import mods.ROOTmanager as manager
import mods.configuration as config
#ROOT.gSystem.Load('/nfs/slac/g/ldmx/users/${USER}/ldmx-sw/install/lib/libFramework.so')
ROOT.gSystem.Load('libFramework.so')

def main():
    
    # Parse command line args and init most important ones
    pdict = manager.parse()
    action_str = pdict['action']
    configFile = pdict['config']

    if action_str == 'trees': actor = BdtTreeMaker
    elif action_str == 'train': actor = BdtTrainer
    elif action_str == 'eval': actor = BdtEval
    else: sys.exit('\nProvide a valid action')

    print('\nUsing {} action from conf file: {}'.format(action_str,configFile))

    # If batch, just do the thing and exit
    if pdict['batch'] and action_str != 'train':
        actor( config.parse_batch_config(pdict) ).run()
        sys.exit('\nDONE!')

    # Parse an print Config and (Overrides parse options if provided in both)
    # Maybe give an overriding message in these cases
    proc_confs = config.parse_bdt_config(
                                        action_str,
                                        configFile,
                                        clargs = pdict
                                        )

    # Wait for confirmation
    if not pdict['batch']:
        input('\nPress Enter to continue... ( Also hi, have a good day :D )')

    # Construct processes
    # Mayhaps multiprocessing for this?
    procs = []
    for proc_config in proc_confs:
        procs.append( actor( proc_config ) )

    # Process processes
    for proc in procs:

        # RUN
        proc.run(
            strEvent = pdict['startEvent'],
            maxEvents = pdict['maxEvents'],
            pfreq = 1000
            )

        del proc # Done with process, remove from memory

    print('\nU DONE GOOD UWU !\n')

##################################################
# "Actors"
##################################################
class BdtTreeMaker(manager.TreeProcess):

    """
    Make flat trees to train on
    Consider multiple inheritance later
    """

    def __init__(self, proc_conf):
        
        """ Set up TreeMaker, define event_proces, and init TreeProccess """

        # Set up tfMaker to have branches_info
        self.tfMaker = manager.TreeMaker(proc_conf.tConfig)

        # Get set of branches
        branches = set()
        for fun_info in proc_conf.tConfig.funcs.values():
            for br_tup in fun_info['brs']:
                branches.add(br_tup)

        # Init parent TreeProcess (sets attrs used in next funcs loop)
        super().__init__(
                        proc_conf.pConfig,
                        branches = branches,
                        endfs=[ self.tfMaker.wq ]
                        )

        # Get lists of dicts containing funcs, args, and priorities
        # for use in appropriate event_process sections

        func_groups = {
                'init':    ( 0, 2),
                'closing': (40,50)
                }

        dets = ('tracker', 'ecal', 'hcal')
        steps = ('_init', '_l1', '_med', '_l2', '_closing')

        for det in dets:
            for step in steps:

                func_groups[ det + step ] =\
                        (10*(dets.index(det)+1) + steps.index(step),
                         10*(dets.index(det)+1) + steps.index(step) + 1)

        # Begin with list containin dummy to simplify appending condition
        for g in func_groups:
            setattr(self, g, [{'func': None}])

        # Check each function
        for fun, info in proc_conf.tConfig.funcs.items():

            # Against func_groups and determine which it belongs in
            for lab, lims in func_groups.items():
                if lims[0] <= info['priority'] < lims[1]:
                    g = lab
                    break

            # If it's new, add its info to its group
            if fun not in [ f for f in getattr(self,g) ]:

                getattr(self,g).append( {
                                'func': fun,
                                'args': { tup[0]: getattr(self, tup[0]) \
                                    for tup in info['brs'] },
                                'priority': info['priority']
                            }
                        )

        # Sort function groups in case of an internal hierarchy, then pop prio
        for g in func_groups:
            getattr(self,g).pop(0)
            if getattr(self,g) != []:
                getattr(self,g).sort(key = lambda d: d['priority'])
                for fund in getattr(self,g):
                    fund.pop('priority')

        def doFeatFuncs(f_dict, funcs_list, store_dict={}, lq=None):

            """ Short function for looping over others in funcs lists """

            for fun in funcs_list:
                fun['func']( f_dict, fun['args'], store_dict, lq)

        def doDetector(f_dict, det, lb_prefix):

            """ Do all doFeatFuncs steps for a given detector """

            store = {}

            # Init loop items
            doFeatFuncs(f_dict, getattr(self, det + '_init'), store)

            if getattr(self, det + '_l1') != []:

                # Get input branch name
                for d in getattr(self, det + '_l1'):
                    for br in d['args']:
                        if br[:len(lb_prefix)] == lb_prefix: loop_branch = br

                # Loop over said branch while doing lfs
                for hit in getattr(self, loop_branch):
                    doFeatFuncs(f_dict, getattr(self, det + '_l1'), store, hit)

                # Do any intermediary functions
                doFeatFuncs(f_dict, getattr(self, det + '_med'), store)

                # Loop again if needed
                if getattr(self, det + '_l2') != []:
                    for hit in getattr(self, loop_branch):
                        doFeatFuncs(f_dict, getattr(self, det + '_l2'), store,
                                    hit)

            # Any further det functions
            doFeatFuncs(f_dict, getattr(self, det + '_closing'), store)

        # Main event algorithm
        def event_process():
            
            """ Algorithm for computing and storing all features """

            # Initialize BDT input variables w/ defaults
            feats = self.tfMaker.resetFeats()

            # Copy from input and other basiic assignments
            doFeatFuncs(feats, self.init)

            # Tracker
            doDetector(feats, 'tracker', 'TrackerRecHits') # ?

            # Ecal
            doDetector(feats, 'ecal', 'EcalRecHits')

            # Hcal
            doDetector(feats, 'hcal', 'HcalRecHits')

            # Any final closing functions
            doFeatFuncs(feats, self.closing)

            self.tfMaker.fillEvent(feats)

        # Tell self about event_process
        setattr(self, 'event_process', event_process)

        # Light-hearted attempt to save on memory
        del branches
        del dets
        del steps
        del func_groups

class BdtTrainer():

    """ Train a BDT """

    def __init__(self, conf_dict):

        """ Set up BDT parameters """

        # Print config
        config.print_dict(conf_dict, prefix='\nBDT configuration:')
        
        # Set config items as attrs
        for k,v in conf_dict.items():
            setattr(self, k, v)

        # Yet another conf dictionary
        self.params_dict = {
                "objective": "binary:logistic",
                "eta": self.eta,
                "max_depth": self.tree_depth,
                "min_child_weight": 20,
                "silent": 1,
                "subsample": .9,
                "colsample_bytree": .85,
                "eval_metric": 'error',
                "seed": 1,
                "nthread": 1,
                "verbosity": 1,
                "early_stopping_rounds" : 10
                }

    def run(self, strEvent=None, maxEvents=1.25e6 , pfreq=None):

        """ Run the traning - startEvent and pfreq are placeholders """

        # import some stuff
        import os
        import logging
        import numpy as np
        import pickle as pkl
        import xgboost as xgb
        import matplotlib as plt

        # Seed and logging 
        np.random.seed( int(self.seed) )
        ml_logger = logging.getLogger('matplotlib')
        ml_logger.setLevel(logging.WARNING)
        plt.use('Agg')

        # Organize data for training
        sets = ('bkg', 'sig')
        for st in sets:

            # Load tree
            tree = manager.load(
                    [getattr(self, '{}_indir'.format(st))],
                    self.tree_name
                    )
 
            events = []
            for event in tree:
                if len(events) == maxEvents: break
                events.append(
                        [ getattr(event, feat) for feat in self.branches ]
                    )

            events = np.array(events)
            new_idx = np.random.permutation(
                    np.arange( np.shape(events)[0] )
                    )
            np.take(events, new_idx, axis = 0, out=events)

            setattr(self, '{}_train_x'.format(st), events)
            setattr(self,
                    '{}_train_y'.format(st),
                    np.zeros(
                        len(
                            getattr( self, '{}_train_x'.format(st) )
                            )
                        ) + (st == 'sig')
                    )

        # Combine data
        train_x = np.vstack(( self.sig_train_x, self.bkg_train_x ))
        train_y = np.append( self.sig_train_y, self.bkg_train_y )
        train_x[ np.isnan( train_x ) ] = 0.000
        train_y[ np.isnan( train_y ) ] = 0.000
        training_matrix = xgb.DMatrix(train_x, train_y)

        # Actual training
        gbm = xgb.train(
                self.params_dict,
                training_matrix,
                int(self.tree_number)
                )

        # Store BDT
        outname = self.outdir.split('/')[-1]
        if not os.path.exists(self.outdir):
            print( 'Creating %s' % (self.outdir) )
            os.makedirs(self.outdir)
        output = open('{}/{}_weights.pkl'.format(self.outdir, outname), 'wb')
        pkl.dump(gbm, output)

        # Plot feature importances
        xgb.plot_importance(gbm)
        plt.pyplot.savefig(
            '{}/{}_fimportance.png'.format(self.outdir, outname), # png name
            dpi=500, bbox_inches='tight', pad_inches=0.5 # png parameters
            )
        
        # Anounce save location
        print('Files saved in: {}'.format(self.outdir))

class BdtEval(manager.TreeProcess):

    """ Evaluate BDT on reserved flat trees """

    def __init__(self, proc_conf):

        # Set up tfMaker to have branches_info
        self.tfMaker = manager.TreeMaker(proc_conf.tConfig)

        # Init parent TreeProcess (sets attrs used in next funcs loop)
        super().__init__(
                        proc_conf.pConfig,
                        endfs=[ self.tfMaker.wq ]
                        )

        # import some stuff
        import numpy as np
        import pickle as pkl
        import xgboost as xgb

        # Set bdt
        self.bdt = pkl.load( open( proc_conf.pConfig.bdt, 'rb' ) )

        # Store discValue name
        for k in self.tfMaker.branches_info.keys():
            if k[:9] == 'discValue': discValue_name = k

        # Main event algorithm
        def event_process():

            # Collect features from this event
            feats = []
            for feat in self.tfMaker.branches_info:
                if feat != discValue_name:
                    feats.append( getattr( self.tree, feat ) )

            # Copy features to new tree
            for f_name, f_value in zip(self.tfMaker.branches_info, feats):
                self.tfMaker.branches[f_name][0] = f_value

            # Add prediction to new tree
            evtarray = np.array([feats])
            self.tfMaker.branches[discValue_name][0] =\
                    float( self.bdt.predict( xgb.DMatrix(evtarray) )[0] )

            # Fill new tree with current event values
            self.tfMaker.tree.Fill()

        # Tell self about event_process
        setattr(self, 'event_process', event_process)

if __name__ == '__main__': main()
