#include "TChain.h"
#include <iostream>
#include "TProof.h"

void my_starter_data(){

  Bool_t isMC = false;
  
  TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/bparknano_30Dec21_Chunk0.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/Chunk0_n500/bparknano_30Dec21_nj1.root";

  TString outFileName = "flat_bparknano.root";
  if(isMC) outFileName += "_isMC";

  TChain * c = new TChain("Events");
  c->Add(inFileName);

  c->Process("TrackMatcher.C+", outFileName);
}
