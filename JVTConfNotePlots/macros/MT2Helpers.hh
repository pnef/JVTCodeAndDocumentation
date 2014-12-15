#ifndef MT2Helpers_HH
#define MT2Helpers_HH

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TROOT.h"

// ---------------- MT2RatioCanvas ---------------------------------------------
class MT2RatioCanvas  {

public:
	MT2RatioCanvas();
	~MT2RatioCanvas();
	TCanvas *canv;
	TPad    *pad1;
	TPad    *pad2;
	void UseMainPad();
	void UseRatioPad();
};

MT2RatioCanvas::MT2RatioCanvas(){
	canv = new TCanvas("c1","Histogram Drawing Options",0,0,600,800);
	pad1 = new TPad("pad1", "",0.0,0.25,1,1   );
	pad2 = new TPad("pad2", "",0.0,0.0, 1,0.25);
	pad1->Draw();
	pad2->Draw();
}

MT2RatioCanvas::~MT2RatioCanvas(){
	if (pad1!=NULL) delete pad1;
	if (pad2!=NULL) delete pad2;
	if (canv!=NULL) delete canv;
	gROOT->SetStyle("MT2-Style");
	gROOT->ForceStyle();
	gROOT->UseCurrentStyle();
}

void MT2RatioCanvas::UseMainPad(){
	gROOT->SetStyle("MT2-Style");
	pad1->cd();
	pad1 ->UseCurrentStyle();
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle();
}

void MT2RatioCanvas::UseRatioPad(){
	gROOT->SetStyle("MT2-ratio-Style");
	pad2 ->cd();
	pad2 ->UseCurrentStyle();
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle();
}

#endif
