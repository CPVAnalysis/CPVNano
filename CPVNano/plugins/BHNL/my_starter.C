#include "TChain.h"
#include <iostream>
#include "TProof.h"
//#include "utils.C"

#include <iomanip>
#include <fstream>
#include <algorithm>

using namespace std;

void my_starter(){

  Bool_t isMC = false;
  Bool_t isSignalMC = true;

  //TString inFileName = "../../test/bparknano_mc_new.root";
  TString inFileName = "../../test/bparknano.root";
  //TString inFileName = "../../test/bparknano_trigger_1file_2.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/bparknano.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018B/Chunk0_n750/bparknano_nj95.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_full/mass3.0_ctau184.256851021/nanoFiles/Chunk0_n750/bparknano_nj95.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_full/mass3.0_ctau184.256851021/nanoFiles/merged/bparknano.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/Chunk0_n500/bparknano_selected_stdtrgmu_full_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_selected.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_selected.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V03/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/Chunk0_n454/bparknano_selected_nj1.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/Chunk0_n500/bparknano_selected_stdtrgmu_full_v1_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH4_Run2018B/Chunk0_n20/merged/bparknano_for_triggermuon_matching_study.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V02/ParkingBPH1_Run2018A/Chunk0_n500/merged/bparknano_selected.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_selected_muononly_alltrgmu.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/Chunk0_n500/bparknano_selected_updatedgenmatching_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_selected_updatedgenmatching.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection_dsamuons_dr0p25.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/Chunk0_n500/bparknano_29Jun21_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V05_29Jun21/ParkingBPH1_Run2018A/Chunk0_n500/bparknano_29Jun21_nj1.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/bparknano_looselection_updatedmatching_fullgen.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/Chunk0_n500/bparknano_looseselection_mupi_v1_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/bparknano_looseselection.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/Chunk0_n129/bparknano_29Jun21_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V06_29Jun21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/Chunk0_n129/bparknano_29Jun21_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/Chunk1_n500/merged/bparknano_selectedslimmed_looseselectiondsa_loosedeltaptrel.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V0/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n500/bparknano_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V07_18Aug21/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/Chunk0_n129/bparknano_18Aug21_nj1.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V09_06Nov21/ParkingBPH1_Run2018A/merged/bparknano_data_1file_looseselection.root";
  
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/bparknano_30Dec21_Chunk0.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n86/bparknano_30Dec21_nj1.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V10_30Dec21/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/Chunk0_n500/bparknano_30Dec21_nj1.root";
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V39_Bc/mass4.5_ctau1.0/nanoFiles/merged/bparknano_24Apr22.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/V13_06Feb23/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/Chunk0_n500/bparknano_06Feb23_nj97.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V12_08Aug22/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n32/bparknano_08Aug22_nj1.root";
  //
  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/bparknano_08Aug22_nj1.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/Chunk0_n500/bparknano_08Aug22_nj1.root";

  //TString inFileName = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V12_08Aug22/ParkingBPH1_Run2018D/merged/bparknano_data_loosedr_1file.root";

  TString outFileName = "flat_bparknano.root";
  if(isSignalMC){
    outFileName += "_isSignalMC";
  }
  if(isMC && !isSignalMC){
    outFileName += "_isMC";

   // string filename = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V10_30Dec21/BToNMuX_NToEMuPi_SoftQCD_b_mN3p0_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n86/bparknano_30Dec21_nj1.root";
   // //string filename = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/mc_central/BToJPsiKstar_V0/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/Chunk0_n500/bparknano_nj1.root";
   // std::size_t index_ini = filename.find("ctau");
   // if((int)index_ini != -1){
   //   std::size_t index_fin = filename.find("mm", index_ini+4);
   //   string ctau = filename.substr(index_ini+4, index_fin-index_ini-4);
   //   TString ctau_sample = ctau.c_str();//.Replace("p", ".");
   //   ctau_sample = ctau_sample.ReplaceAll("p", ".");
   //   //Ssiz_t index_ini = inFileName.Index("ctau")+4;
   //   //std::cout << "ctau " << index_ini << std::endl;
   //   outFileName += "_isMC_" + ctau_sample;
   // }
   // else{
   //   outFileName += "_isMC";
   // }
  }
      
  TChain * c = new TChain("Events");
  c->Add(inFileName);

  //c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/T6_updatedPU_n4200000_njt200/mass3.0_ctau1013.41268062/nanoFiles/Chunk0_n2/bparknano_test_run_nj1.root");
  //c->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/T6_updatedPU_n4200000_njt200/mass3.0_ctau1013.41268062/nanoFiles/Chunk0_n2/bparknano_test_run_nj2.root");
  

  //c->Process("NanoDumper.C+", outFileName);
  c->Process("BToMuMuPiDumper.C+", outFileName);
  //c->Process("BToKMuMuDumper.C+", outFileName);
  //c->Process("TagAndProbeDumper.C+", outFileName);
  //c->Process("HNLToMuPiDumper.C+", outFileName);

  if(isMC){
    TChain* c_run = new TChain("Runs");
    c_run->Add(inFileName);
    //c_run->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/T6_updatedPU_n4200000_njt200/mass3.0_ctau1013.41268062/nanoFiles/Chunk0_n2/bparknano_test_run_nj1.root");
    //c_run->Add("/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/T6_updatedPU_n4200000_njt200/mass3.0_ctau1013.41268062/nanoFiles/Chunk0_n2/bparknano_test_run_nj2.root");

    c_run->Process("NanoRunDumper.C", outFileName);
  }
}
