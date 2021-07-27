# Python Implimentation of BDT

## Purpose: Test potential of ML using multiple subdetectors, especially in the case of visible signal decays
## Requirments: Working install of `ldmx-sw-v2.3.0` or greater and `v12` samples.
* If using v2.X.X, change `libFramework.so` to `libEvent.so` in `mainframe.py`.
* `confs/gabreille_back_v1.ini` is used here as an example config. See its comments for further explanation.

In the directory containing `ldmx-sw`, `$LDMX_BASE`, enter `source ldmx-sw/scripts/ldmx-env.sh`.
* You may clone `ldmx-sw` with a different name (e.g. `ldmx-sw-v3.0.0`), just chnage it in the above command and in `ROOT.gSystem.Load()` in `mainframe.py`

### Interactive
To make flat trees from ldmx event files:
```
ldmx python3 mainframe.py trees confs/test.ini
```
* `-m` can be used to give maximum number of events to run over 
* `-p` can be used to selected a subset of processeses listed in the config file
* Other options can be found in `mods/ROOTmanager.py` and `mods/configuration.py`

To train BDT:
`ldmx hadd` flat trees into the paths given in the config file
```
ldmx python3 mainfraim.py train confs/test.ini
```

To evaluate trained BDT on test samples:
```
ldmx python3 mainframe.py eval confs/test.ini
```

### Batch
```cp batch.py $LDMX_BASE; cd $LDMX_BASE``` \
Remove `ldmx` and replace `maineframe.py` in corresponding interactive commands with `batch.py`.\
(Some options are removed at this stage to avoid making big mistakes; add them at your own risk.)
