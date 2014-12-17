#ifndef __TrigRoILinksFillerTool_h__
#define __TrigRoILinksFillerTool_h__
/*
  TrigRoILinksFillerTool.h
*/
#include "D3PDMakerUtils/BlockFillerTool.h"

namespace D3PD {
  
  class TrigRoILinksFillerTool : public BlockFillerTool<RoILinks> {
  public:
    TrigRoILinksFillerTool(const std::string& type, 
			   const std::string& name, 
			   const IInterface* parent);
    ~TrigRoILinksFillerTool();
    
    StatusCode book();
    StatusCode fill(const RoILinks& x);
    
  private:
    int* mRoIType;
    std::vector<int*> mObjIndexVec;
    std::vector<std::vector<int>* > mObjIndicesVec;

    std::vector<std::string> mBranchNames;
    std::vector<std::string> mBranchNames2;
  };

}

#endif // __TrigRoILinks[FillerTool_h__
