#ifndef EventBuilder_jetmet2012pileupcustom_h
#define EventBuilder_jetmet2012pileupcustom_h

#include "EventBuilder_jetmet2012.h"
#include "Particle.h"

class EventBuilder_jetmet2012pileupcustom : public EventBuilder_jetmet2012{

public:
  bool DoCustomSetup(); // everything done on top of EventBuilder_jetmet2012 is called from here
  
  bool fRecoverBtracks;

  // Selector Utils
  Bool_t TrackSelector(); 
  Bool_t  TrackVertexLinks();

  // vertex calculations
  void VertexCalculations();

  // JVF related methods
  void CalculateJVF(const MomKey jetType, const MomKey trackType);

  // truth Matching
  void AddTruthMatch(const MomKey JetType, const MomKey TruthJetType, const bool doLabel, bool highestPt);
  void DoTrackTruthMatching();

  // jet helpers
  void JetSelector(const MomKey JetType);
  bool JetIsIsolated(const MomKey JetType, Particle* jet, float minDR);
  void PartonFlavorMatch(const MomKey JetType);



private:

};



#endif
