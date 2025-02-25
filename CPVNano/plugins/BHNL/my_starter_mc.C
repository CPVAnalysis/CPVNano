#include "TChain.h"
#include <iostream>
#include "TProof.h"

void my_starter_mc(){

  Bool_t isMC = true;

  TString outFileName = "flat_bparknano.root";
  if(isMC) outFileName += "_isMC";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root";
  TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root";

  TChain * c = new TChain("Events");
  //c->Add(inFileName);

  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");

  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");

  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN2p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");

  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");

  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");
  c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN4p5_ctau1p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano_30Dec21.root");

  c->Process("TrackMatcher.C+", outFileName);
}
