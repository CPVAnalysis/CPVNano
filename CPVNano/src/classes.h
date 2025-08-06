#ifndef BPHNANO_CLASSES_H
#define BPHNANO_CLASSES_H

#include "DataFormats/Common/interface/Wrapper.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "PhysicsTools/CPVNano/interface/KinVtxFitter.h"
#include <vector>


namespace {
  struct dictionary {
      edm::Wrapper<std::vector<KinVtxFitter> > wkv;
  };
}

#endif  
