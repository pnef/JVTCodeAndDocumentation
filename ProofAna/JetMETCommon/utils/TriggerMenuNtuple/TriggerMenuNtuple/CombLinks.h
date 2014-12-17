#ifndef __CombLinks_h__
#define __CombLinks_h__
/*
  CombLinks.h
*/
#include "TriggerMenuNtuple/FeatureIndex.h"
#include <vector>
#include <string>
#include <map>

class CombLinks {
public:
  static int featureId(const std::string& feature);
  static std::string featureName(int feature_id);
  static int addFeature(const std::string& feature);
  static int addRoIType(const std::string& feature);
  static const std::map<int, std::string>& featureIdMap();
private:
  static std::map<int, std::string> sFeatureIdMap;
  friend std::ostream& operator<<(std::ostream& o, const CombLinks& x);

public:
  typedef std::vector<FeatureIndex> FeatureIndexVec_t;
  
public:
  CombLinks();
  CombLinks(int roi_type);
  ~CombLinks();

  bool hasFeature(const std::string& feature) const;
  std::vector<std::string> allFeatureNames() const;
  std::vector<std::string> FeatureNames() const;
  std::vector<std::string> FeatureVecNames() const;
  int lastStep() const { return mLastStep; }

  int RoIType() const { return mRoIType; }
  const FeatureIndex* index(const std::string& feature) const;
  const std::vector<FeatureIndex>* indexVec(const std::string& feature) const;

  void setLastStep(int i) { mLastStep = i; }
  void setRoIType(int i) { mRoIType = i; }
  void addIndex(const std::string& feature, const FeatureIndex& i);
  void addIndexVec(const std::string& feature, const FeatureIndexVec_t& iv);
  
  bool isValid() const;
  void dump() const;
  bool isSameRoI(const CombLinks& x) const;
  bool operator==(const CombLinks& x) const;

  const FeatureIndex* index(int feature_id) const;
  const std::vector<FeatureIndex>* indexVec(int feature_id) const;

private:
  int mRoIType;
  int mLastStep;
  std::map<int, FeatureIndex> mIndexMap;
  std::map<int, FeatureIndexVec_t> mIndexVecMap;
};

std::ostream& operator<<(std::ostream& o, const CombLinks& x);

#endif // __CombLinks_h__
