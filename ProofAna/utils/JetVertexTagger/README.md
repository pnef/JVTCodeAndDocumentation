JetVertex Tagger RootCore packate
=======================
This is a RootCore package I've written to use JVT in run-I physics analysis as a root core package with NTUP_COMMON samples. 

The jet-vertex tagger JVT is as a 2-dimensional likelihood to tag and suppress pileup jets. Similarly to JVF, it used tracking and vertexing information to compute a variable to discriminate hard-scatter from pileup jets. JVT has the following advantages over JVF:
Performance: JVT is more effective than JVF at suppressing pileup. For a large range of cut values resulting hard-scatter jet efficiencies, the corresponding fraction of pileup jets passing such cuts (pileup fake rate) is significantly lower for JVT.
NPV independent hard-scatter jet efficiencies: The input variables of JVT, namely corrJVF and RpT are constructed to be insensitive to the pileup activity in the event. Consequently, cutting on JVT to suppress pileup jets results in hard-scatter jet efficiencies that are pileup independent.

Note: This RootCore implementation of JVT is only approximate, based on the information available in the NTUP_COMMON D3PD branches, leading to slight performance differences w.r.t. the CONF Note version of JVT. Please see details under Package Validation below.

The complete JVT documentation can be found under ATLAS-CONF-2014-018.

The ATLAS internal twiki page is here: [JVT Twiki](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetVertexTagger)

# Usage Instructions

## Package and Tag For standalone, Rootcore and Athena analyses:
```
svn co $SVNOFF/Reconstruction/Jet/JetAnalysisTools/JetVertexTagger/tags/JetVertexTagger-00-00-01 JetVertexTagger
```

The code is Standalone/athena/RootCore compatible.

## To compile it Standalone:
```
cd Reconstruction/Jet/JetAnalysisTools/JVFUncertaintyTool/cmt
cmt make -f Makefile.Standalone
```

It is then possible to link the compiled .so to your standalone executable.

## How to use the tool?

### In the header file
Include the header file, initializing the JetVertexTagger

```
#include "JetVertexTagger/JetVertexTagger.h"

JetVertexTagger* jvt;
```

### Before the event loop
Call the constructor: specify the JVT cut value and point to the JVT likelihood file

```
jvt = new JetVertexTagger(0.2, [path-to]/JetVertexTagger/data/JVTlikelihood_20140805.root");
```

### In the event loop
Initialize the tool with the list of tracks (pt, z0_wrtPV) and track-to-vertex association (vxp_trk_index). These vectors should be taken directly from the NTUP COMMON D3PD branches, but the track pT need to be converted to GeV!

```
jvt->init(trk_pt, trk_z0_wrtPV, vxp_trk_index);
```

### For each jet
To check if a jet passes the JVT cut, pass the jet pT and the associated tracks to the tool. Here, assoc_trk_indices is a vector of indices corresponding to the indices of the associated tracks in the vectors trk_pt and trk_z0_wrtPV . It should be taken directly from the D3PD branch.

WARNING: The indices of vxp_trk_index and assoc_trk_indices has to correspond to the ones from trk_pt and trk_z0_wrtPV!

```
bool pass = (*jvt)(myjet_pt, assoc_trk_indices);
```

where myjet_pt is the fully calibrated jet pT, in GeV!
The values of corrJVF, RpT and JVT for the jet can then be accessed as (not necessary for tagging with a specified cut):

```
float corrJVF = jvt->corrJVF();
float RpT = jvt->RpT();
float JVT = jvt->JVT();
```

## Test your installation!
We prepared a simple example explicitly implementing the steps described above based on a simple MakeClass generated on a sample D3PD file containing a few events. We recommend this example is used as a guideline of how to use the JetVertexTagger tool.

### Running the example
To run the example on all 100 events in the file mc12_8TeV.147771.Sherpa_CT10_Zmumu.merge.NTUP_COMMON.e1434_s1499_s1504_r3658_r3549_p1562.NTUP_COMMON.01313919._000880.skimmed100Events.root , do the following.

```
root -b -q JVTTest.C 
```

This will only run after you're successfully compiled the JetVertexTagger package and created the shared library StandAlone/libJetVertexTagger.so.
The script will print to the stdout a list of event number, jet pTs, JVT, corrJVF and RpT values in the following form:

```
JVTTest: >>> EventNumber 4141202 jet pt       24.3 eta       2.31 jvt       0.73 corrJVF      0.868 RpT       0.15 <<< 
JVTTest: >>> EventNumber 4141202 jet pt       22.9 eta     -0.162 jvt          1 corrJVF      0.965 RpT       1.81 <<< 
JVTTest: >>> EventNumber 4141202 jet pt       20.1 eta       2.14 jvt       0.05 corrJVF      0.376 RpT      0.101 <<< 
JVTTest: >>> EventNumber 4141204 jet pt       22.5 eta    -0.0479 jvt       0.97 corrJVF      0.873 RpT      0.267 <<< 
JVTTest: >>> EventNumber 4141206 jet pt       27.5 eta     -0.355 jvt          1 corrJVF      0.948 RpT      0.666 <<< 
JVTTest: >>> EventNumber 4141208 jet pt       37.8 eta      0.271 jvt          1 corrJVF      0.943 RpT      0.768 <<< 
JVTTest: >>> EventNumber 4141214 jet pt       43.6 eta      -1.88 jvt       0.99 corrJVF      0.959 RpT      0.908 <<< 
JVTTest: >>> EventNumber 4141214 jet pt       27.9 eta      -1.12 jvt          1 corrJVF       0.93 RpT       2.38 <<< 
JVTTest: >>> EventNumber 4141216 jet pt         36 eta      -1.13 jvt       0.01 corrJVF      0.156 RpT     0.0555 <<< 
JVTTest: >>> EventNumber 4141216 jet pt         33 eta      0.205 jvt       0.95 corrJVF      0.699 RpT      0.414 <<< 
JVTTest: >>> EventNumber 4141216 jet pt       31.4 eta      0.551 jvt       0.22 corrJVF      0.395 RpT      0.205 <<< 
JVTTest: >>> EventNumber 4141217 jet pt       21.6 eta      -1.73 jvt       0.96 corrJVF      0.822 RpT      0.257 <<< 
JVTTest: >>> EventNumber 4141219 jet pt       29.6 eta       1.01 jvt          1 corrJVF      0.976 RpT       1.33 <<< 
JVTTest: >>> EventNumber 4141219 jet pt         23 eta       1.88 jvt          0 corrJVF      0.191 RpT     0.0568 <<< 
```

In this script, the jet pT is taken from the D3PD at the constituent scale. This is so that the results of this script can easily be used as a reference, independently of the version of the JES that is applied. In an actual implementation of the tool in a physics analysis, the fully calibrated jet pT should be used as input!
We strongly suggest you test your own implementation of the JetVertexTagger tool against this script, running your analysis on the same D3PD input file.



