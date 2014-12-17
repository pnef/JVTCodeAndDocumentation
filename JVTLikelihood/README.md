JVTLikelihood
=======================

Here you find instructions on how to rederive the JVT likelihood. 


## Step 1: Construct a root tree containing the RpT and corrJVT information for each jet

The tree should be filled on a jet-by-jet basis (rather than event-by-event) so that each entry in the tree corresponds to a new jet. (This makes it easier to feed the tree to TMVA.). 

As an example -- take a look at the Analysis [Analysis_PileUpStudiesTreeFiller.cxx](../ProofAna/PileUpStudies/analyses/Analysis_PileUpStudiesTreeFiller.cxx):

 * corrJVT and RpT can be computed using the JVT RootCore package. See [here](..//ProofAna/utils/JetVertexTagger) for information and usage instructions. 

## Step 2: run TMVA

The file [TMVA_Train_JVT_likelihood.py](TMVA_Train_JVT_likelihood.py) is a python script I used to construct the kNN discriminator using TMVA. 
The (self-explaining) arguments / input files / variables are defined in the header of the file, but can also be passed as command line arguments. 

```
DEFAULT_OUTFNAME         = "skimmedJpt20to50Eta2p4.20140606.13.00_PileUpStudies_rev360973.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN100trimmed.nPUtrkCorrJVF_RpT.JHStrkPtSumGE0.root"
DEFAULT_INFNAME          = "/atlas/local/pnef//Pileup/Data/skimmedJpt20to50Eta2p4.20140606.13.00_PileUpStudies_rev360973.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.root"
DEFAULT_TREESIG          = "JetTree0"
DEFAULT_TREEBKG          = "JetTree0"
DEFAULT_METHODS          = "kNN100trim,likelihood"
DEFAULT_NEVENTS_TEST_S   =0
DEFAULT_NEVENTS_TEST_B   =0
DEFAULT_NEVENTS_TRAIN_S  =0
DEFAULT_NEVENTS_TRAIN_B  =0
 
VARIABLES             = ["corrJVF_RootCore", "RpT_RootCore"]
SELECTION             = "fabs(Jeta)<2.4 && Jpt>20 && Jpt<50 && RpT_RootCore<2 && VtxDzTruth<0.1 && RpT_RootCore>0"
SPECTATORS            = ["NPV", "Mu" ,"Jpt", "Jeta", "JisPU", "JisHS", "Weight", "VtxDzTruth", "Jtruthpt",
                          "JJVF", "JCorrPUtrkPtSumOverPt", "nPUTracks"]
```
Run the script as 
```
python TMVA_Train_JVT_likelihood.py
```
and this will produce an output tree containing the specified variables are branches (VARIABLES and SPECTATORS) including the new discriminator. Also, a xml file is produced, containing the JVT information for each (RpT, corrJVF) pair of the training sample of jets. In the above example config, this file will be
```
skimmedJpt20to50Eta2p4.20140606.13.00_PileUpStudies_rev360973.PythJ1and2mc12aJETMET.jetmet2012pileupcustom.kNN100trimmed.nPUtrkCorrJVF_RpT.JHStrkPtSumGE0_KNN100trim.weights.xml
```

## Step 3: Construct the JVT Likelihood histogram

The script [JVT_makeLikelihood.py](JVT_makeLikelihood.py) serves two purposes:
 * produce a 2D histogram, showing the JVT value as a function of corrJVT and RpT
 * read RpT and corrJVT from a jet-by-jet tree, compute JVT based on the previously produced histogram, and add a JVT branch to the tree. 
 

#### producing the histogram
Running ``` JVT_makeLikelihood.py --makeLikelihood ``` will take the specified TMVA-xml file, load the discriminator,  scan over all possible corrJVF and RpT values and compute the resulting JVT value. This is then filled in a 2D histogram, producing a plot of this kind:
![Fig 5a] (https://raw.githubusercontent.com/pnef/JVTCodeAndDocumentation/master/JVTConfNotePlots/plots/fig_05a.png)





