"""
# Standard tree-building dictionaries

## Analysis quantities

## Base stuff for testing

## Visible decay dictionaries next because thesis
    ### Bases
        - backv1
        - sidev1 (coming soon)
        - trackerv1 (coming soon)

    ### Combos
        - hcal_vx (back_vx + side_vx)
        - hcal_tracker_vx (hcal_vx + tracker_vx)
        - full_vxG (tracker_vx + ecal(Gabrielle) + hcal_vx)
        - full_vxSy (tracker_vx + ecal(Segmipy) + hcal_vx)

## Gabrielle because standardization

## SegCont

## MipTracking

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
from mods import analysis
from mods.ecal import segcont
from mods.ecal import mipTracking
from mods.ecal import main as emain

# Production process label (definatley going into config ... eventually)
pp = 'v12'

##################################################
# Tree-building dictionaries
##################################################

# Analysis quantities (NOT FOR TRAINING)
##################################################
analysis_funcs = {
                    #analysis.trigger: {
                    #        'priority': 20.1,
                    #        'brs:' (
                    #            ('Trigger_'+pp, 'TriggerResult'),
                    #            )
                    #    },
                    #analysis.recoil_tracks: {
                    #    'priority': 0,
                    #    'brs': (
                    #        ('EcalScoringPlaneHits_'+pp, 'SimTrackerHit'),
                    #        )
                    #    },
                    emain.ecal_init: {
                        'priority': 20,
                        'brs': (
                            ('TargetScoringPlaneHits_'+pp, 'SimTrackerHit'),
                            ('EcalScoringPlaneHits_'+pp, 'SimTrackerHit')
                            )
                        },
                    analysis.fidcats: {
                        'priority': 20.1,
                        'brs': ()
                        },
                    analysis.recoilPT: {
                        'priority': 20.1,
                        'brs': ()
                        },
                    analysis.decay_verts: {
                        'priority':  1,
                        'brs': (
                            ('SimParticles_'+pp, 'SimParticle'),
                            )
                        }
        }

trees_info_analysis = {
            #'n_recoil_tracks': {'rtype': int,  'default': None },
            #'max_track_pman':  {'rtype': int,  'default': None },
            'e_fid':           {'rtype': bool, 'default': False},
            'g_fid':           {'rtype': bool, 'default': False},
            'recoilPT':        {'rtype': float,'default': None },
            'Ap_dt':           {'rtype': float,'default': None },
            'Ap_dx':           {'rtype': float,'default': None },
            'Ap_dy':           {'rtype': float,'default': None },
            'Ap_dz':           {'rtype': float,'default': None },
        }

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
maxPE_funcs = {
                    hcal.maxPE: {
                        'priority': 30,
                        'brs': (
                            ('HcalVeto_'+pp, 'HcalVetoResult'),
                            ('HcalRecHits_'+pp, 'HcalHit')
                            )
                        },
                    hcal.maxPEsearch: {
                        'priority': 31,
                        'brs': (
                            ('HcalRecHits_'+pp, 'HcalHit'),
                            )
                        }
                }

trees_info_maxPE = {
        'maxPE': {'rtype': int, 'default': 0 },
        }

# Intemded for Back HCal (section 0) backv1
# "Full" list
# avg_x/y commented out in case of pT concerns, but may not be relevant in
# visible decay scenarios anyway
# ("Intended", because easily extendable to other sections)
##################################################

# All of back section
backv1_all_funcs = {
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
                        hcal.backv1_all: {
                            'priority': 32,
                            'brs': ()
                            }
                    }

trees_info_backv1_all = {
        'back_nHits':   {'rtype': int,   'default': 0.},
        'back_totE':    {'rtype': float, 'default': 0.},
        'back_totPE':   {'rtype': float, 'default': 0.},
        'back_maxE':    {'rtype': float, 'default': 0.},
        'back_maxPE':   {'rtype': int,   'default': 0.},
        'back_avgE':    {'rtype': float, 'default': 0.},
        'back_avgPE':   {'rtype': float, 'default': 0.},
        'back_std_e_x': {'rtype': float, 'default': 0.},
        'back_std_e_y': {'rtype': float, 'default': 0.},
        'back_std_e_z': {'rtype': float, 'default': 0.},
        'back_dz_e':    {'rtype': float, 'default': 0.},
        }

# Back section segments
backv1_seg_funcs = {
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
                        hcal.backv1_seg: {
                            'priority': 32,
                            'brs': ()
                            }
                    }

trees_info_backv1_seg = {}
for s in range(1, hcal.back_segments + 1):

    trees_info_backv1_seg[f'back_nHits_{s}e']  = {'rtype': int,   'default': 0 }
    trees_info_backv1_seg[f'back_tot_{s}e_e']  = {'rtype': float, 'default': 0.}
    trees_info_backv1_seg[f'back_tot_{s}e_pe'] = {'rtype': int,   'default': 0 }
    trees_info_backv1_seg[f'back_max_{s}e_e']  = {'rtype': float, 'default': 0.}
    trees_info_backv1_seg[f'back_max_{s}e_pe'] = {'rtype': int,   'default': 0 }
    trees_info_backv1_seg[f'back_avg_{s}e_e']  = {'rtype': float, 'default': 0.}
    trees_info_backv1_seg[f'back_avg_{s}e_pe'] = {'rtype': float, 'default': 0.}
    trees_info_backv1_seg[f'back_std_{s}e_e']  = {'rtype': float, 'default': 0.}
    trees_info_backv1_seg[f'back_std_{s}e_x']  = {'rtype': float, 'default': 0.}
    trees_info_backv1_seg[f'back_std_{s}e_y']  = {'rtype': float, 'default': 0.}
    trees_info_backv1_seg[f'back_std_{s}e_pe'] = {'rtype': float, 'default': 0.}

backv1_funcs = {
            **maxPE_funcs,
            **backv1_all_funcs,
            **backv1_seg_funcs
        }

trees_info_backv1 = {
                    **trees_info_maxPE,
                    **trees_info_backv1_all,
                    **trees_info_backv1_seg
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

for r in range(1, emain.contRegions + 1):

    trees_info_gabrielle[f'electronContainmentEnergy_x{r}'] = {'rtype': float, 'default': 0 }
    trees_info_gabrielle[f'photonContainmentEnergy_x{r}']   = {'rtype': float, 'default': 0.}
    trees_info_gabrielle[f'outsideContainmentEnergy_x{r}']  = {'rtype': float, 'default': 0.}
    trees_info_gabrielle[f'outsideContainmentNHits_x{r}']   = {'rtype': int,   'default': 0 }
    trees_info_gabrielle[f'outsideContainmentXStd_x{r}']    = {'rtype': float, 'default': 0.}
    trees_info_gabrielle[f'outsideContainmentYStd_x{r}']    = {'rtype': float, 'default': 0.}

# SegCont
##################################################
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

trees_info_segcont = trees_info_base_ecal.copy()
for s in range(1, segcont.nSegments + 1):

    # Longitudinal segment features
    trees_info_segcont[f'energy_s{s}']    = {'rtype': float, 'default': 0.}
    trees_info_segcont[f'nHits_s{s}']     = {'rtype': int,   'default': 0 }
    trees_info_segcont[f'xMean_s{s}']     = {'rtype': float, 'default': 0.}
    trees_info_segcont[f'yMean_s{s}']     = {'rtype': float, 'default': 0.}
    trees_info_segcont[f'layerMean_s{s}'] = {'rtype': float, 'default': 0.}
    trees_info_segcont[f'xStd_s{s}']      = {'rtype': float, 'default': 0.}
    trees_info_segcont[f'yStd_s{s}']      = {'rtype': float, 'default': 0.}
    trees_info_segcont[f'layerStd_s{s}']  = {'rtype': float, 'default': 0.}

    for r in range(1, emain.contRegions + 1):

        # Electron RoC features
        trees_info_segcont[f'eContEnergy_x{r}_s{s}']    = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'eContNHits_x{r}_s{s}']     = {'rtype': int,   'default': 0 }
        trees_info_segcont[f'eContXMean_x{r}_s{s}']     = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'eContYMean_x{r}_s{s}']     = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'eContLayerMean_x{r}_s{s}'] = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'eContXStd_x{r}_s{s}']      = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'eContYStd_x{r}_s{s}']      = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'eContLayerStd_x{r}_s{s}']  = {'rtype': float, 'default': 0.}

        # Photon RoC features
        trees_info_segcont[f'gContEnergy_x{r}_s{s}']    = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'gContNHits_x{r}_s{s}']     = {'rtype': int,   'default': 0 }
        trees_info_segcont[f'gContXMean_x{r}_s{s}']     = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'gContYMean_x{r}_s{s}']     = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'gContLayerMean_x{r}_s{s}'] = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'gContXStd_x{r}_s{s}']      = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'gContYStd_x{r}_s{s}']      = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'gContLayerStd_x{r}_s{s}']  = {'rtype': float, 'default': 0.}

        # Outside RoC features
        trees_info_segcont[f'oContEnergy_x{r}_s{s}']    = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'oContNHits_x{r}_s{s}']     = {'rtype': int,   'default': 0 }
        trees_info_segcont[f'oContXMean_x{r}_s{s}']     = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'oContYMean_x{r}_s{s}']     = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'oContLayerMean_x{r}_s{s}'] = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'oContXStd_x{r}_s{s}']      = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'oContYStd_x{r}_s{s}']      = {'rtype': float, 'default': 0.}
        trees_info_segcont[f'oContLayerStd_x{r}_s{s}']  = {'rtype': float, 'default': 0.}

# RsegCont
##################################################
rsegcont_funcs = {
            emain.base_rsegcont: {
                'priority': 20,
                'brs': (
                    ('EcalVeto_'+pp, 'EcalVetoResult'),
                    )
                },
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
            segcont.rsegcont_means: {
                'priority': 21,
                'brs': (
                    ('EcalRecHits_'+pp, 'EcalHit'),
                    ),
                },
            segcont.rsegcont_means_norm: {
                'priority': 22,
                'brs': ()
                },
            segcont.rsegcont_stds: {
                'priority': 23,
                'brs': (
                    ('EcalRecHits_'+pp, 'EcalHit'),
                    ),
                },
            segcont.rsegcont_stds_norm: {
                'priority': 24,
                'brs': ()
                }
        }
"""
segcont.prep_ecal_lfs: {
    'priority': 20.2,
    'brs': ()
    },
segcont.collect: {
    'priority': 21,
    'brs': (
        ('EcalRecHits_'+pp, 'EcalHit'),
        )
    },
segcont.rsegcont: {
    'priority': 22,
    'brs': ()
    }
"""

trees_info_rsegcont = {
        'nReadoutHits':    {'rtype': int,   'default': 0 },
        'summedDet':       {'rtype': float, 'default': 0.},
        'summedTightIso':  {'rtype': float, 'default': 0.},
        'maxCellDep':      {'rtype': float, 'default': 0.},
        'showerRMS':       {'rtype': float, 'default': 0.},
        'xStd':            {'rtype': float, 'default': 0.},
        'yStd':            {'rtype': float, 'default': 0.},
        'stdLayerHit':     {'rtype': float, 'default': 0.},
        }

for s in range(1, segcont.nSegments + 1):

    # Longitudinal segment features Should probably remove x/yMeans for pT bias
    trees_info_rsegcont[f'energy_rs{s}']    = {'rtype': float, 'default': 0.}
    trees_info_rsegcont[f'nHits_rs{s}']     = {'rtype': int,   'default': 0 }
    trees_info_rsegcont[f'xMean_rs{s}']     = {'rtype': float, 'default': 0.}
    trees_info_rsegcont[f'yMean_rs{s}']     = {'rtype': float, 'default': 0.}
    trees_info_rsegcont[f'xStd_rs{s}']      = {'rtype': float, 'default': 0.}
    trees_info_rsegcont[f'yStd_rs{s}']      = {'rtype': float, 'default': 0.}
    trees_info_rsegcont[f'layerStd_rs{s}']  = {'rtype': float, 'default': 0.}

    for r in range(1, emain.contRegions + 1):

        # Electron RoC features
        trees_info_rsegcont[f'eContEnergy_x{r}_rs{s}']    = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'eContNHits_x{r}_rs{s}']     = {'rtype': int,   'default': 0 }
        trees_info_rsegcont[f'eContXMean_x{r}_rs{s}']     = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'eContYMean_x{r}_rs{s}']     = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'eContXStd_x{r}_rs{s}']      = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'eContYStd_x{r}_rs{s}']      = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'eContLayerStd_x{r}_rs{s}']  = {'rtype': float, 'default': 0.}

        # Photon RoC features
        trees_info_rsegcont[f'gContEnergy_x{r}_rs{s}']    = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'gContNHits_x{r}_rs{s}']     = {'rtype': int,   'default': 0 }
        trees_info_rsegcont[f'gContXMean_x{r}_rs{s}']     = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'gContYMean_x{r}_rs{s}']     = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'gContXStd_x{r}_rs{s}']      = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'gContYStd_x{r}_rs{s}']      = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'gContLayerStd_x{r}_rs{s}']  = {'rtype': float, 'default': 0.}

        # Outside RoC features
        trees_info_rsegcont[f'oContEnergy_x{r}_rs{s}']    = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'oContNHits_x{r}_rs{s}']     = {'rtype': int,   'default': 0 }
        trees_info_rsegcont[f'oContXMean_x{r}_rs{s}']     = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'oContYMean_x{r}_rs{s}']     = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'oContXStd_x{r}_rs{s}']      = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'oContYStd_x{r}_rs{s}']      = {'rtype': float, 'default': 0.}
        trees_info_rsegcont[f'oContLayerStd_x{r}_rs{s}']  = {'rtype': float, 'default': 0.}

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

##################################################
# Feature dictionaries
##################################################

# Bases
##################################################

# Analysis (Not actually feats; ignored in training and eval actors)
feats_analysis = trees_info_analysis

# Tracker

# Hcal
feats_maxPE = trees_info_maxPE

feats_backv1_all = trees_info_backv1_all

feats_backv1_seg = trees_info_backv1_seg

feats_backv1 = trees_info_backv1

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

# Make rsegmip dicts smartly
feats_rsegmipv3 = {

            # Base features
            'nReadoutHits':    {'rtype': int,   'default': 0 },
            'summedDet':       {'rtype': float, 'default': 0.},
            'summedTightIso':  {'rtype': float, 'default': 0.},
            'maxCellDep':      {'rtype': float, 'default': 0.},
            'showerRMS':       {'rtype': float, 'default': 0.},
            'xStd':            {'rtype': float, 'default': 0.},
            'yStd':            {'rtype': float, 'default': 0.},
            'stdLayerHit':     {'rtype': float, 'default': 0.},

            # MIP tracking features
            'straight4':                 {'rtype': int,   'default': 0 },
            'firstNearPhLayer':          {'rtype': int,   'default': 33},
            'nNearPhHits':               {'rtype': int,   'default': 0 },
            'photonTerritoryHits':       {'rtype': int,   'default': 0 },
            'epSep':                     {'rtype': float, 'default': 0.},
            'epDot':                     {'rtype': float, 'default': 0.},
    }
for k,v in feats_segmipv3.items():
    if '_s' in k:
        candidate_k = k.replace('_s','_rs')
        if candidate_k in trees_info_rsegcont:
            feats_rsegmipv3[ candidate_k ] = v

feats_rsegmipv4 = {

            # Base features
            'nReadoutHits':    {'rtype': int,   'default': 0 },
            'summedDet':       {'rtype': float, 'default': 0.},
            'summedTightIso':  {'rtype': float, 'default': 0.},
            'maxCellDep':      {'rtype': float, 'default': 0.},
            'showerRMS':       {'rtype': float, 'default': 0.},
            'xStd':            {'rtype': float, 'default': 0.},
            'yStd':            {'rtype': float, 'default': 0.},
            'stdLayerHit':     {'rtype': float, 'default': 0.},

            # MIP Tracking features
            **feats_mipTracking,
    }
for k,v in feats_segmipv4.items():
    if '_s' in k:
        candidate_k = k.replace('_s','_rs')
        if candidate_k in trees_info_rsegcont:
            feats_rsegmipv4[ candidate_k ] = v

feats_rsegmipx = {
                    feat: info for feat, info in feats_rsegmipv3.items() \
                        if feat in feats_rsegmipv4
    }
