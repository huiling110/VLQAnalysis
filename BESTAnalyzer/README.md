# BESTAnalyzer: An Ntuplizer that uses BEST

This portion of the code makes ntuples using a trained BEST network. In `src/BESTAnalyzer.cc` is the main code.
The code for loading the neural network using tensorflow is in `interface` and the tools for making the BEST input 
variables are in `src/BESTToolbox.cpp`

## Instructions

To run this code locally, use one of the python run files:

```bash
cmsRun test/run_ZZ.py
```


