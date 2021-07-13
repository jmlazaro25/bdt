"""
# Standard tree-building dictionaries

## Base stuff for testing

## Visible decay dictionaries next because thesis
    ### Bases
        - back_v1
        - side_v1 (coming soon)
        - tracker_v1 (coming soon)

    ### Combos
        - hcal_vx (back_vx + side_vx)
        - hcal_tracker_vx (hcal_vx + tracker_vx)
        - full_vxG (tracker_vx + ecal(Gabrielle) + hcal_vx)
        - full_vxSy (tracker_vx + ecal(Segmipy) + hcal_vx)

## Gabrielle because standardization

## SegCont

## MipTracking

## Analysis quantities

Priorities (floats allowed for finer precision):

    0 : Needed for event classification (All)
    1 : Easy to compute, no prereqs, just do (All)

    if Tracker: det = 1
    if Ecal:    det = 2
    if Hcal:    det = 3

    if Must be done before det loops: step = 0
    if Done in first det loop:        step = 1
    if Done between det loops:        step = 2
    if Done in second det loop:       step = 3
    if Finish det:                    step = 4

    10*det + step: "main" (e.g. Done in loop 1 of Hcal: 31)

    40: Finish all

# Standard feature dictioanries (similar organization)

"""

# NOTE: Maybe split feats in feature
# NOTE: If needed, see git-history for giving each feat its own funcs dict
#       removed for brievity

from mods import hcal
from mods.ecal import segcont
from mods.ecal import mipTracking
from mods.ecal import main as emain

# Production process label
pp = 'v12'

##################################################
# Tree-building dictionaries
##################################################

# Base stuff (mostly) for testing
##################################################

# Base tracker (findableTracks)

# Base ecal
base_ecal_funcs = {
                    emain.base_ecal: {
                        'priority': 1,
                        'brs': (
                            ('EcalVeto_'+pp, 'EcalVetoResult'),
                            )
                        }
                }

trees_info_base_ecal = {
        'nReadoutHits':    {'rtype': int,   'default': 0 },
        'summedDet':       {'rtype': float, 'default': 0.},
        'summedTightIso':  {'rtype': float, 'default': 0.},
        'maxCellDep':      {'rtype': float, 'default': 0.},
        'showerRMS':       {'rtype': float, 'default': 0.},
        'xStd':            {'rtype': float, 'default': 0.},
        'yStd':            {'rtype': float, 'default': 0.},
        'avgLayerHit':     {'rtype': float, 'default': 0.},
        'stdLayerHit':     {'rtype': float, 'default': 0.},
        'deepestLayerHit': {'rtype': int,   'default': 0 },
        'ecalBackEnergy':  {'rtype': float, 'default': 0.},
        }

# Base hcal (maxPE)
hcal_base_funcs = {
                    hcal.maxPE: {
                        'priority': 1,
                        'brs': (
                            ('HcalVeto_'+pp, 'HcalVetoResult'),
                            ('HcalRecHits_'+pp, 'HcalHit')
                            )
                        }
                }

trees_info_maxPE = {
        'maxPE': {'rtype': int, 'default': 0 },
        }

# Intemded for Back HCal (section 0) back_v1
# "Full" list
# avg_x/y commented out in case of pT concerns, but may not be relevant in
# visible decay scenarios anyway
# ("Intended", because easily extendable to other sections)
##################################################

# All of back section
hcal_back_v1_all_funcs = {
                            hcal.prep_hcal_lfs: {
                                'priority': 30,
                                'brs': ()
                                },
                            # HCalRecHits not technically used in this function
                            # but it's got to be somewhere and this seems like
                            # the best place
                            hcal.collect: {
                                'priority': 31,
                                'brs': (
                                    ('HcalRecHits_'+pp, 'HcalHit'),
                                    )
                                },
                            hcal.back_v1_all: {
                                'priority': 32,
                                'brs': ()
                                }
                        }

trees_info_back_v1 = {                                                               
        'back_nHits':   {'rtype': int,   'default': 0.},
        'back_totE':    {'rtype': float, 'default': 0.},
        'back_totPE':   {'rtype': float, 'default': 0.},
        'back_maxE':    {'rtype': float, 'default': 0.},
        'back_maxPE':   {'rtype': int,   'default': 0.},
        'back_avgE':    {'rtype': float, 'default': 0.},
        'back_avgPE':   {'rtype': float, 'default': 0.},
        'back_avgR':    {'rtype': float, 'default': 0.},
        'back_std_e_x': {'rtype': float, 'default': 0.},
        'back_std_e_y': {'rtype': float, 'default': 0.},
        'back_std_e_z': {'rtype': float, 'default': 0.},
        'back_dz_e':    {'rtype': float, 'default': 0.},
        }

# Global (insludes sides)
trees_info_back_v1.update( trees_info_maxPE )

# Back section segments
hcal_back_v1_seg_funcs = {
                            hcal.prep_hcal_lfs: {
                                'priority': 30,
                                'brs': ()
                                },
                            hcal.collect: {
                                'priority': 31,
                                'brs': (
                                    ('HcalRecHits_'+pp, 'HcalHit'),
                                    )
                                },
                            hcal.back_v1_seg: {
                                'priority': 33,
                                'brs': ()
                                }
                        }

for i in range(1, hcal.back_segments + 1):

    trees_info_back_v1['back_nHits_{}e'.format(i)]  = {'rtype': int,   'default': 0.}
    trees_info_back_v1['back_tot_{}e_e'.format(i)]  = {'rtype': float, 'default': 0.}
    trees_info_back_v1['back_tot_{}e_pe'.format(i)] = {'rtype': int,   'default': 0.}
    trees_info_back_v1['back_max_{}e_e'.format(i)]  = {'rtype': float, 'default': 0.}
    trees_info_back_v1['back_max_{}e_pe'.format(i)] = {'rtype': int,   'default': 0.}
    trees_info_back_v1['back_avg_{}e_e'.format(i)]  = {'rtype': float, 'default': 0.}
    #trees_info_back_v1['back_avg_{}e_x'.format(i)]  = {'rtype': float, 'default': 0.} # Maybe pT bias
    #trees_info_back_v1['back_avg_{}e_y'.format(i)]  = {'rtype': float, 'default': 0.} # Don't train
    trees_info_back_v1['back_avg_{}e_pe'.format(i)] = {'rtype': float, 'default': 0.}
    trees_info_back_v1['back_avg_{}e_r'.format(i)]  = {'rtype': float, 'default': 0.}
    trees_info_back_v1['back_std_{}e_e'.format(i)]  = {'rtype': float, 'default': 0.}
    trees_info_back_v1['back_std_{}e_x'.format(i)]  = {'rtype': float, 'default': 0.}
    trees_info_back_v1['back_std_{}e_y'.format(i)]  = {'rtype': float, 'default': 0.}
    trees_info_back_v1['back_std_{}e_pe'.format(i)] = {'rtype': float, 'default': 0.}

back_v1_funcs = {
            **hcal_base_funcs,
            **hcal_back_v1_all_funcs,
            **hcal_back_v1_seg_funcs
        }

# Gabrielle
##################################################
trees_info_gabrielle = trees_info_base_ecal.copy()

gabrielle_funcs = {
                    **base_ecal_funcs,
                    emain.gabrielle_containment: {
                        'priority': 1,
                        'brs': (
                            ('EcalVeto_'+pp, 'EcalVetoResult'),
                            )
                        }
                }

for i in range(1, emain.contRegions + 1):

    trees_info_gabrielle['electronContainmentEnergy_x{}'.format(i)] = {'rtype': float, 'default': 0 }
    trees_info_gabrielle['photonContainmentEnergy_x{}'.format(i)]   = {'rtype': float, 'default': 0.}
    trees_info_gabrielle['outsideContainmentEnergy_x{}'.format(i)]  = {'rtype': float, 'default': 0.}
    trees_info_gabrielle['outsideContainmentNHits_x{}'.format(i)]   = {'rtype': int,   'default': 0 }
    trees_info_gabrielle['outsideContainmentXStd_x{}'.format(i)]    = {'rtype': float, 'default': 0.}
    trees_info_gabrielle['outsideContainmentYStd_x{}'.format(i)]    = {'rtype': float, 'default': 0.}


# SegCont
##################################################
trees_info_segcont = trees_info_base_ecal.copy()

segcont_funcs = {
            emain.ecal_init: {
                'priority': 20,
                'brs': (
                    ('TargetScoringPlaneHits_'+pp, 'SimTrackerHit'),
                    ('EcalScoringPlaneHits_'+pp, 'SimTrackerHit')
                    )
                },
            segcont.segcont_init: {
                'priority': 20.1,
                'brs': ()
                },
            segcont.segcont_means: {
                'priority': 21,
                'brs': (
                    ('EcalRecHits_'+pp, 'EcalHit'),
                    ),
                },
            segcont.segcont_means_norm: {
                'priority': 22,
                'brs': ()
                },
            segcont.segcont_stds: {
                'priority': 23,
                'brs': (
                    ('EcalRecHits_'+pp, 'EcalHit'),
                    ),
                },
            segcont.segcont_stds_norm: {
                'priority': 24,
                'brs': ()
                }
        }

for s in range(1, segcont.nSegments + 1):
    
    # Longitudinal segment features
    trees_info_segcont['energy_s{}'.format(s)]    = {'rtype': float, 'default': 0.}
    trees_info_segcont['nHits_s{}'.format(s)]     = {'rtype': int,   'default': 0 }
    trees_info_segcont['xMean_s{}'.format(s)]     = {'rtype': float, 'default': 0.}
    trees_info_segcont['yMean_s{}'.format(s)]     = {'rtype': float, 'default': 0.}
    trees_info_segcont['layerMean_s{}'.format(s)] = {'rtype': float, 'default': 0.}
    trees_info_segcont['xStd_s{}'.format(s)]      = {'rtype': float, 'default': 0.}
    trees_info_segcont['yStd_s{}'.format(s)]      = {'rtype': float, 'default': 0.}
    trees_info_segcont['layerStd_s{}'.format(s)]  = {'rtype': float, 'default': 0.}

    for r in range(1, emain.contRegions + 1):

        # Electron RoC features
        trees_info_segcont['eContEnergy_x{}_s{}'.format(r,s)]    = {'rtype': float, 'default': 0.}
        trees_info_segcont['eContNHits_x{}_s{}'.format(r,s)]     = {'rtype': int,   'default': 0 }
        trees_info_segcont['eContXMean_x{}_s{}'.format(r,s)]     = {'rtype': float, 'default': 0.}
        trees_info_segcont['eContYMean_x{}_s{}'.format(r,s)]     = {'rtype': float, 'default': 0.}
        trees_info_segcont['eContLayerMean_x{}_s{}'.format(r,s)] = {'rtype': float, 'default': 0.}
        trees_info_segcont['eContXStd_x{}_s{}'.format(r,s)]      = {'rtype': float, 'default': 0.}
        trees_info_segcont['eContYStd_x{}_s{}'.format(r,s)]      = {'rtype': float, 'default': 0.}
        trees_info_segcont['eContLayerStd_x{}_s{}'.format(r,s)]  = {'rtype': float, 'default': 0.}

        # Photon RoC features
        trees_info_segcont['gContEnergy_x{}_s{}'.format(r,s)]    = {'rtype': float, 'default': 0.}
        trees_info_segcont['gContNHits_x{}_s{}'.format(r,s)]     = {'rtype': int,   'default': 0 }
        trees_info_segcont['gContXMean_x{}_s{}'.format(r,s)]     = {'rtype': float, 'default': 0.}
        trees_info_segcont['gContYMean_x{}_s{}'.format(r,s)]     = {'rtype': float, 'default': 0.}
        trees_info_segcont['gContLayerMean_x{}_s{}'.format(r,s)] = {'rtype': float, 'default': 0.}
        trees_info_segcont['gContXStd_x{}_s{}'.format(r,s)]      = {'rtype': float, 'default': 0.}
        trees_info_segcont['gContYStd_x{}_s{}'.format(r,s)]      = {'rtype': float, 'default': 0.}
        trees_info_segcont['gContLayerStd_x{}_s{}'.format(r,s)]  = {'rtype': float, 'default': 0.}

        # Outside RoC features
        trees_info_segcont['oContEnergy_x{}_s{}'.format(r,s)]    = {'rtype': float, 'default': 0.}
        trees_info_segcont['oContNHits_x{}_s{}'.format(r,s)]     = {'rtype': int,   'default': 0 }
        trees_info_segcont['oContXMean_x{}_s{}'.format(r,s)]     = {'rtype': float, 'default': 0.}
        trees_info_segcont['oContYMean_x{}_s{}'.format(r,s)]     = {'rtype': float, 'default': 0.}
        trees_info_segcont['oContLayerMean_x{}_s{}'.format(r,s)] = {'rtype': float, 'default': 0.}
        trees_info_segcont['oContXStd_x{}_s{}'.format(r,s)]      = {'rtype': float, 'default': 0.}
        trees_info_segcont['oContYStd_x{}_s{}'.format(r,s)]      = {'rtype': float, 'default': 0.}
        trees_info_segcont['oContLayerStd_x{}_s{}'.format(r,s)]  = {'rtype': float, 'default': 0.}


# mipTracking
##################################################
mipTracking_funcs = {
            emain.ecal_init: {
                'priority': 20,
                'brs': (
                    ('TargetScoringPlaneHits_'+pp, 'SimTrackerHit'),
                    ('EcalScoringPlaneHits_'+pp, 'SimTrackerHit')
                    )
                },
            mipTracking.mipTracking_init: {
                'priority': 20.1,
                'brs': ()
                },
            mipTracking.epSepAndDot: {
                'priority': 20.2,
                'brs': ()
                },
            mipTracking.mipTracking_collect: {
                'priority': 21,
                'brs': (
                    ('EcalRecHits_'+pp, 'EcalHit'),
                    )
                },
            mipTracking.fullTerritories: {
                'priority': 21,
                'brs': (
                    ('EcalRecHits_'+pp, 'EcalHit'),
                    )
                },
            mipTracking.territories: {
                'priority': 22,
                'brs': ()
                },
            mipTracking.nearPhotonInfo: {
                'priority': 22,
                'brs': ()
                },
            mipTracking.tracks: {
                'priority': 22,
                'brs': ()
                },
        }

trees_info_mipTracking = {
            'straight4':                      {'rtype': int,   'default': 0 },
            'firstNearPhLayer':               {'rtype': int,   'default': 33},
            'nNearPhHits':                    {'rtype': int,   'default': 0 },
            'fullElectronTerritoryHits':      {'rtype': int,   'default': 0 },
            'fullPhotonTerritoryHits':        {'rtype': int,   'default': 0 },
            'fullTerritoryRatio':             {'rtype': float, 'default': 1.},
            'electronTerritoryHits':          {'rtype': int,   'default': 0 },
            'photonTerritoryHits':            {'rtype': int,   'default': 0 },
            'TerritoryRatio':                 {'rtype': float, 'default': 1.},
            'epSep':                          {'rtype': float, 'default': 0.},
            'epDot':                          {'rtype': float, 'default': 0.},
        }

# Analysis quantities (NOT FOR TRAINING)
# Should be added at end of flat trees to minimize impact on training feats
# Consider working in sample_analysis stuff
##################################################
"""
branches_info_analysis = {
        'trigPass':
            {
                'rtype': int,
                'default': 1 ,
                'func': 'copy_from_iput'
            },
        'recoilPT':
            {
                'rtype': float,
                'default': 0.,
                'func': 'recoilPT'
            }
        # Noise next
        }
"""
##################################################
# Feature dictionaries
##################################################

# Bases
##################################################

# Tracker

# Hcal
feats_maxPE = trees_info_maxPE 

feats_back_v1 = trees_info_back_v1

# Ecal
feats_base_ecal = trees_info_base_ecal

feats_gabrielle = trees_info_gabrielle

feats_mipTracking = trees_info_mipTracking

feats_gabmip = { **feats_gabrielle, **feats_mipTracking}

feats_segmipv3 = {

            # Base features
            **feats_base_ecal,

            # MIP tracking features
            'straight4':                 {'rtype': int,   'default': 0 },
            'firstNearPhLayer':          {'rtype': int,   'default': 33},
            'nNearPhHits':               {'rtype': int,   'default': 0 },
            'photonTerritoryHits':       {'rtype': int,   'default': 0 },
            'epSep':                     {'rtype': float, 'default': 0.},
            'epDot':                     {'rtype': float, 'default': 0.},

            # Longitudinal segment features
            'energy_s1':                 {'rtype': float, 'default': 0.},
            'xMean_s1':                  {'rtype': float, 'default': 0.},
            'yMean_s1':                  {'rtype': float, 'default': 0.},
            'layerMean_s1':              {'rtype': float, 'default': 0.},
            'xStd_s1':                   {'rtype': float, 'default': 0.},
            'yStd_s1':                   {'rtype': float, 'default': 0.},
            'layerStd_s1':               {'rtype': float, 'default': 0.},
            'energy_s2':                 {'rtype': float, 'default': 0.},
            'nHits_s2':                  {'rtype': int,   'default': 0 },
            'yMean_s2':                  {'rtype': float, 'default': 0.},
            'energy_s3':                 {'rtype': float, 'default': 0.},
            'yMean_s3':                  {'rtype': float, 'default': 0.},

            # Electron RoC features
            'eContEnergy_x1_s1':         {'rtype': float, 'default': 0.},
            'eContEnergy_x2_s1':         {'rtype': float, 'default': 0.},
            'eContYMean_x1_s1':          {'rtype': float, 'default': 0.},
            'eContEnergy_x1_s2':         {'rtype': float, 'default': 0.},
            'eContEnergy_x2_s2':         {'rtype': float, 'default': 0.},
            'eContNHits_x1_s2':          {'rtype': int,   'default': 0 },
            'eContYMean_x1_s2':          {'rtype': float, 'default': 0.},
            'eContEnergy_x1_s3':         {'rtype': float, 'default': 0.},

            # Photon RoC features
            'gContNHits_x1_s1':          {'rtype': int,   'default': 0 },
            'gContYMean_x1_s1':          {'rtype': float, 'default': 0.},
            'gContEnergy_x1_s2':         {'rtype': float, 'default': 0.},
            'gContEnergy_x2_s2':         {'rtype': float, 'default': 0.},
            'gContEnergy_x3_s2':         {'rtype': float, 'default': 0.},
            'gContNHits_x1_s2':          {'rtype': int,   'default': 0 },
            'gContYMean_x1_s2':          {'rtype': float, 'default': 0.},

            # Outside RoC features
            'oContEnergy_x1_s1':         {'rtype': float, 'default': 0.},
            'oContEnergy_x2_s1':         {'rtype': float, 'default': 0.},
            'oContEnergy_x3_s1':         {'rtype': float, 'default': 0.},
            'oContNHits_x1_s1':          {'rtype': int,   'default': 0 },
            'oContXMean_x1_s1':          {'rtype': float, 'default': 0.},
            'oContYMean_x1_s1':          {'rtype': float, 'default': 0.},
            'oContYMean_x2_s1':          {'rtype': float, 'default': 0.},
            'oContLayerMean_x1_s1':      {'rtype': float, 'default': 0.},
            'oContYStd_x1_s1':           {'rtype': float, 'default': 0.},
            'oContYStd_x2_s1':           {'rtype': float, 'default': 0.},
            'oContEnergy_x1_s2':         {'rtype': float, 'default': 0.},
            'oContEnergy_x2_s2':         {'rtype': float, 'default': 0.},
            'oContEnergy_x3_s2':         {'rtype': float, 'default': 0.},
            'oContNHits_x2_s2':          {'rtype': int,   'default': 0 },
            'oContXMean_x1_s2':          {'rtype': float, 'default': 0.},
            'oContYMean_x1_s2':          {'rtype': float, 'default': 0.},
            'oContYMean_x2_s2':          {'rtype': float, 'default': 0.},
            'oContYMean_x3_s2':          {'rtype': float, 'default': 0.},
            'oContLayerMean_x1_s2':      {'rtype': float, 'default': 0.},
            'oContYStd_x1_s2':           {'rtype': float, 'default': 0.},
            'oContYStd_x2_s2':           {'rtype': float, 'default': 0.},
            'oContLayerStd_x1_s2':       {'rtype': float, 'default': 0.},
            'oContEnergy_x1_s3':         {'rtype': float, 'default': 0.},
            'oContLayerMean_x1_s3':      {'rtype': float, 'default': 0.},
    }

feats_segmipv4 = {

            # Base features
            **feats_base_ecal,

            # MIP Tracking features
            **feats_mipTracking,

            # Longitudinal segment features
            'energy_s1':                      {'rtype': float, 'default': 0.},
            'nHits_s1':                       {'rtype': int,   'default': 0 },
            'xMean_s1':                       {'rtype': float, 'default': 0.},
            'yMean_s1':                       {'rtype': float, 'default': 0.},
            'layerMean_s1':                   {'rtype': float, 'default': 0.},
            'energy_s2':                      {'rtype': float, 'default': 0.},
            'xStd_s2':                        {'rtype': float, 'default': 0.},
            'yStd_s2':                        {'rtype': float, 'default': 0.},
            'layerStd_s2':                    {'rtype': float, 'default': 0.},
            'nHits_s3':                       {'rtype': int,   'default': 0 },
            'xMean_s3':                       {'rtype': float, 'default': 0.},
            'yMean_s3':                       {'rtype': float, 'default': 0.},
            'layerMean_s3':                   {'rtype': float, 'default': 0.},
            'xStd_s3':                        {'rtype': float, 'default': 0.},
            'yStd_s3':                        {'rtype': float, 'default': 0.},
            'layerStd_s3':                    {'rtype': float, 'default': 0.},

            # Electron RoC features
            'eContEnergy_x1_s1':              {'rtype': float, 'default': 0.},
            'eContEnergy_x2_s1':              {'rtype': float, 'default': 0.},
            'eContNHits_x1_s1':               {'rtype': int,   'default': 0 },
            'eContNHits_x2_s1':               {'rtype': int,   'default': 0 },
            'eContYMean_x1_s1':               {'rtype': float, 'default': 0.},
            'eContYMean_x2_s1':               {'rtype': float, 'default': 0.},
            'eContLayerMean_x1_s1':           {'rtype': float, 'default': 0.},
            'eContLayerMean_x2_s1':           {'rtype': float, 'default': 0.},
            'eContEnergy_x1_s2':              {'rtype': float, 'default': 0.},
            'eContEnergy_x2_s2':              {'rtype': float, 'default': 0.},
            'eContXMean_x1_s2':               {'rtype': float, 'default': 0.},
            'eContXMean_x2_s2':               {'rtype': float, 'default': 0.},
            'eContYMean_x1_s2':               {'rtype': float, 'default': 0.},
            'eContYMean_x2_s2':               {'rtype': float, 'default': 0.},

            # Photon RoC features
            'gContEnergy_x1_s1':              {'rtype': float, 'default': 0.},
            'gContEnergy_x2_s1':              {'rtype': float, 'default': 0.},
            'gContNHits_x1_s1':               {'rtype': int,   'default': 0 },
            'gContNHits_x2_s1':               {'rtype': int,   'default': 0 },
            'gContYMean_x1_s1':               {'rtype': float, 'default': 0.},
            'gContYMean_x2_s1':               {'rtype': float, 'default': 0.},
            'gContLayerMean_x1_s1':           {'rtype': float, 'default': 0.},
            'gContLayerMean_x2_s1':           {'rtype': float, 'default': 0.},
            'gContNHits_x1_s2':               {'rtype': int,   'default': 0 },
            'gContNHits_x2_s2':               {'rtype': int,   'default': 0 },
            'gContXMean_x1_s2':               {'rtype': float, 'default': 0.},
            'gContXMean_x2_s2':               {'rtype': float, 'default': 0.},
            'gContLayerMean_x1_s2':           {'rtype': float, 'default': 0.},
            'gContLayerMean_x2_s2':           {'rtype': float, 'default': 0.},
            'gContEnergy_x1_s3':              {'rtype': float, 'default': 0.},
            'gContEnergy_x2_s3':              {'rtype': float, 'default': 0.},

            # Outside RoC features
            'oContEnergy_x1_s1':              {'rtype': float, 'default': 0.},
            'oContEnergy_x2_s1':              {'rtype': float, 'default': 0.},
            'oContEnergy_x3_s1':              {'rtype': float, 'default': 0.},
            'oContNHits_x1_s1':               {'rtype': int,   'default': 0 },
            'oContNHits_x2_s1':               {'rtype': int,   'default': 0 },
            'oContNHits_x3_s1':               {'rtype': int,   'default': 0 },
            'oContXMean_x1_s1':               {'rtype': float, 'default': 0.},
            'oContXMean_x2_s1':               {'rtype': float, 'default': 0.},
            'oContXMean_x3_s1':               {'rtype': float, 'default': 0.},
            'oContYMean_x1_s1':               {'rtype': float, 'default': 0.},
            'oContYMean_x2_s1':               {'rtype': float, 'default': 0.},
            'oContYMean_x3_s1':               {'rtype': float, 'default': 0.},
            'oContXStd_x1_s1':                {'rtype': float, 'default': 0.},
            'oContXStd_x2_s1':                {'rtype': float, 'default': 0.},
            'oContXStd_x3_s1':                {'rtype': float, 'default': 0.},
            'oContYStd_x1_s1':                {'rtype': float, 'default': 0.},
            'oContYStd_x2_s1':                {'rtype': float, 'default': 0.},
            'oContYStd_x3_s1':                {'rtype': float, 'default': 0.},
            'oContLayerStd_x1_s1':            {'rtype': float, 'default': 0.},
            'oContLayerStd_x2_s1':            {'rtype': float, 'default': 0.},
            'oContLayerStd_x3_s1':            {'rtype': float, 'default': 0.},
            'oContEnergy_x1_s2':              {'rtype': float, 'default': 0.},
            'oContEnergy_x2_s2':              {'rtype': float, 'default': 0.},
            'oContEnergy_x3_s2':              {'rtype': float, 'default': 0.},
            'oContLayerMean_x1_s2':           {'rtype': float, 'default': 0.},
            'oContLayerMean_x2_s2':           {'rtype': float, 'default': 0.},
            'oContLayerMean_x3_s2':           {'rtype': float, 'default': 0.},
            'oContLayerStd_x1_s2':            {'rtype': float, 'default': 0.},
            'oContLayerStd_x2_s2':            {'rtype': float, 'default': 0.},
            'oContLayerStd_x3_s2':            {'rtype': float, 'default': 0.},
            'oContEnergy_x1_s3':              {'rtype': float, 'default': 0.},
            'oContEnergy_x2_s3':              {'rtype': float, 'default': 0.},
            'oContEnergy_x3_s3':              {'rtype': float, 'default': 0.},
            'oContLayerMean_x1_s3':           {'rtype': float, 'default': 0.},
            'oContLayerMean_x2_s3':           {'rtype': float, 'default': 0.},
            'oContLayerMean_x3_s3':           {'rtype': float, 'default': 0.},
    }

# segmipx is the intersection of segmipv3 and segmipv4
feats_segmipx = { 
                    feat: info for feat, info in feats_segmipv3.items() \
                        if feat in feats_segmipv4
        }
