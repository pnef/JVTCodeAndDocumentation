/**************************************************************************
 **
 **   File:         Analysis_WritePRW.cxx
 **
 **   Description:  Writes out events in an Event tree
 **                 
 **
 **   Author:       B. Butler
 **
 **   Created:      1-20-12
 **   Modified:
 **
 **************************************************************************/

#define Analysis_WritePRW_cxx

#include "Analysis_WritePRW.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <iostream>
#include <cstdlib>
#include "PileupReweighting/TPileupReweighting.h"


void Analysis_WritePRW::WorkerBegin()
{
  if(Debug()) cout << " working begin for writePRW" << endl;
	prwConfig = 0;
	  if(Debug()) cout << " working begin ending for writePRW" << endl;
}

bool Analysis_WritePRW::ProcessEvent()
{
  if(isMC() && !prwConfig)
		prwConfig = (Root::TPileupReweighting*) InputList()->FindObject("PRWConfigure");



	static const MomKey EventWeight("EventWeight"), averageIntPerXing("averageIntPerXing");

	if(isMC()){
		//prwConfig->Fill(RunNumber(),ChannelNumber(),Float(EventWeight),Float(averageIntPerXing));
	}
	return kTRUE;
	
}

void Analysis_WritePRW::WorkerTerminate()
{
  if(Debug()) cout << "worker terminate for write prw begin " << endl;
	//if(isMC() && prwConfig){
		//prwConfig->WriteToFile();
	//}
  if(Debug()) cout << "worker terminate for write prw end" << endl;
}

