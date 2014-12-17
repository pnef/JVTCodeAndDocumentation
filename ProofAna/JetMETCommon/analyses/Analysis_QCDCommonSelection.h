/**************************************************************************
 **
 **   File:         Analysis_QCDCommonSelection.h
 **
 **   Description:  Analysis class for common selection of top-like events
 **                 
 ** 
 **   Authors:      M. Swiatlowski, B. Nachman
 **
 **   Created:      2012-01-30
 **
 **************************************************************************/

#ifndef Analysis_QCDCommonSelection_h
#define Analysis_QCDCommonSelection_h

#include "Analysis_JetMET_Base.h"
 
using std::cout;
using std::endl;

namespace Root{
  class TTileTripReader;
}



class Analysis_QCDCommonSelection : public Analysis_JetMET_Base {

 public :
  
  Analysis_QCDCommonSelection(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }

  virtual ~Analysis_QCDCommonSelection() { }
  
  ClassDef(Analysis_QCDCommonSelection, 0);
  
  Bool_t  fDetail;

  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();
  void    Book(const TString& prefix);
  void    BookBasic();
  void    BookDetail();
   
 private :			  

  Root::TTileTripReader* m_treader;

};

#endif

