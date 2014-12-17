JetVertex Tagger RootCore packate
=======================
This is a RootCore package I've written to use JVT in run-I physics analysis as a root core package with NTUP_COMMON samples. 

The jet-vertex tagger JVT is as a 2-dimensional likelihood to tag and suppress pileup jets. Similarly to JVF, it used tracking and vertexing information to compute a variable to discriminate hard-scatter from pileup jets. JVT has the following advantages over JVF:
Performance: JVT is more effective than JVF at suppressing pileup. For a large range of cut values resulting hard-scatter jet efficiencies, the corresponding fraction of pileup jets passing such cuts (pileup fake rate) is significantly lower for JVT.
NPV independent hard-scatter jet efficiencies: The input variables of JVT, namely corrJVF and RpT are constructed to be insensitive to the pileup activity in the event. Consequently, cutting on JVT to suppress pileup jets results in hard-scatter jet efficiencies that are pileup independent.

Note: This RootCore implementation of JVT is only approximate, based on the information available in the NTUP_COMMON D3PD branches, leading to slight performance differences w.r.t. the CONF Note version of JVT. Please see details under Package Validation below.

The complete JVT documentation can be found under ATLAS-CONF-2014-018.



