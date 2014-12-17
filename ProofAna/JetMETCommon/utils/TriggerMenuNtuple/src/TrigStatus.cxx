/*
  TrigStatus.cxx
*/
#include "TriggerMenuNtuple/TrigStatus.h"

TrigStatus::TrigStatus(const std::string& name, int status) : 
  mName(name), mStatus(status) {
}

TrigStatus::~TrigStatus() {
}
