#ifndef __JVT__
#define __JVT__

#include <vector>
#include <TString.h>
#include "../JetVertexTagger/JetVertexTagger.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class JetVertexTagger+;
#pragma link C++ class JVTTrack+;
#pragma link C++ class std::vector<JVTTrack>+;

#endif

#endif
