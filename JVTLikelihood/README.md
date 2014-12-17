JVTLikelihood
=======================

Here you find instructions on how to rederive the JVT likelihood. 


## Step 1:

Construct a root tree containing the RpT and corrJVT information for each jet. The tree should be filled on a jet-by-jet basis (rather than event-by-event) so that each entry in the tree corresponds to a new jet. (This makes it easier to feed the tree to TMVA.). 

As an example -- take a look at the Analysis [Analysis_PileUpStudiesTreeFiller.cxx](../ProofAna/PileUpStudies/analyses/Analysis_PileUpStudiesTreeFiller.cxx):

 * corrJVT and RpT can be computed using the JVT RootCore package. See [here](..//ProofAna/utils/JetVertexTagger) for information and usage instructions. 

## Step 2 
