/**************************************************************************
 **
 **   File:         Analysis_WritePRW.h
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

#ifndef Analysis_WritePRW_h
#define Analysis_WritePRW_h

#include "AnalysisBase.h"

namespace Root {class TPileupReweighting;}

using std::cout;
using std::endl;

class Analysis_WritePRW : public AnalysisBase {
public :

	Analysis_WritePRW(TTree* /*tree*/ =0) { } 
	virtual ~Analysis_WritePRW() { }

	ClassDef(Analysis_WritePRW,0);

private :
	
	bool    ProcessEvent();
	void    WorkerBegin(); 
	void    WorkerTerminate();
	
	Root::TPileupReweighting* prwConfig;
	
};

#endif

