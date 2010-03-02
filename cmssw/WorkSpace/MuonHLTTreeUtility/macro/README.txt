The following contains some examples of macros used for getting results from the MHTU.  These macros were developed originally on the fly for quick analysis and will hopefully be updated and refined.

The following scripts are used for
 - L1 and L2 efficiencies, resolutions, and distributions
 - Examining L2 muons with zero valid hits
 - Obtaining trigger efficiencies for all 3 levels
 - Obtaining reconstruction efficiencies for all 3 levels
 - Examining trigger rates by parent type expected in early running.
 - Examining isolation at L2 and L3

In addition to these scripts, the "LibDict" directory, which containes the files that create libTreeLib.so.  This library file allows the macros to process some of the data types used by the MHTU that root otherwise struggles with.  This directory also contains a README file detailing how to create and load this library.

These scripts are all run via the following set of commands

TChain f("MuTrigMC")
f->Add("<your files>")
.L <one of the macros>
ScanChain(MuTrigMC,"<output file name>") // ScanChain(MuTrigMC,"<output file name>",<process type>) for isolation macro

All of these scripts may involve some hands-on fine-tuning for best results.  A partial list of the parameters the user might need to play with is below.

All macros:
 - Including the command gSystem->Load("libSmatrix") might be necessary for running some of these macros.

Reco & trigger efficiencies:
 - Binning for efficiency vs. pT
 - Reversion from "poor man's association by Delta R" to association variables.  This was due to a "feature" in 3.X, where the official CMSSW associators did not account for the new feature of invalid hits at L2.

Trigger rates: 
 - Luminosity (10^30 by default), cross section, filter efficiency

L1 and L2 quantities:
 - Binning for pT-dependent variables
