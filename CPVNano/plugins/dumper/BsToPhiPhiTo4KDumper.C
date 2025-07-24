#define BsToPhiPhiTo4KDumper_cxx
// The class definition in BsToPhiPhiTo4KDumper.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("BsToPhiPhiTo4KDumper.C")
// root> T->Process("BsToPhiPhiTo4KDumper.C","some options")
// root> T->Process("BsToPhiPhiTo4KDumper.C+")
//


#include "BsToPhiPhiTo4KDumper.h"
#include <TMath.h>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "utils.C"

void BsToPhiPhiTo4KDumper::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   cout << " -------------------------------" << endl;
   cout << "       Bs->PhiPhi->4K Dumper    " << endl;
   cout << " -------------------------------" << endl;
}

void BsToPhiPhiTo4KDumper::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   TString outFileName = option;

   isSignalMC = false;
   isMC = false;

   if(outFileName.Contains("isSignalMC")){
      isMC = true;
      isSignalMC = true;
      outFileName.Resize(outFileName.Length()-11);
   }
   else if(outFileName.Contains("isMC")){
      isMC = true;
      outFileName.Resize(outFileName.Length()-5);
   }

   // check if outputfile exists
   if(gSystem->AccessPathName(outFileName)){
      my_file = new TFile(outFileName, "RECREATE");  
   }
   else{
      my_file = new TFile(outFileName, "UPDATE");  
   }
   my_file->cd();

   // if MC, get the correct content of the GenPart branches
   if(isMC){
      nGenPart = {fReader, "nGenPart"};
      GenPart_eta = {fReader, "GenPart_eta"};
      GenPart_mass = {fReader, "GenPart_mass"};
      GenPart_phi = {fReader, "GenPart_phi"};
      GenPart_pt = {fReader, "GenPart_pt"};
      GenPart_vx = {fReader, "GenPart_vx"};
      GenPart_vy = {fReader, "GenPart_vy"};
      GenPart_vz = {fReader, "GenPart_vz"};
      GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
      GenPart_pdgId = {fReader, "GenPart_pdgId"};
      GenPart_status = {fReader, "GenPart_status"};
      GenPart_statusFlags = {fReader, "GenPart_statusFlags"};
      Muon_genPartIdx = {fReader, "Muon_genPartIdx"};
      Pileup_nPU = {fReader, "Pileup_nPU"};
      Pileup_nTrueInt = {fReader, "Pileup_nTrueInt"};
   }

   signal_tree = new TTree("signal_tree", "signal_tree");

   signal_tree->Branch("event", &the_event);
   signal_tree->Branch("run", &the_run);
   signal_tree->Branch("lumi", &the_lumi);
   signal_tree->Branch("pu_ntrueint", &the_pu_ntrueint);
   signal_tree->Branch("pv_npvs", &the_pv_npvs);

   signal_tree->Branch("ntriggermuon", &the_ntriggermuon);

   signal_tree->Branch("bs_beta", &the_bs_beta);
   signal_tree->Branch("bs_eta", &the_bs_eta);
   signal_tree->Branch("bs_charge", &the_bs_charge);
   signal_tree->Branch("bs_cos2d", &the_bs_cos2d);
   signal_tree->Branch("bs_cxx", &the_bs_cxx);
   signal_tree->Branch("bs_cyx", &the_bs_cyx);
   signal_tree->Branch("bs_cyy", &the_bs_cyy);
   signal_tree->Branch("bs_czx", &the_bs_czx);
   signal_tree->Branch("bs_czy", &the_bs_czy);
   signal_tree->Branch("bs_czz", &the_bs_czz);
   signal_tree->Branch("bs_lxysig", &the_bs_lxysig);
   signal_tree->Branch("bs_mass", &the_bs_mass);
   signal_tree->Branch("bs_mass_err", &the_bs_mass_err);
   signal_tree->Branch("bs_mass_corr", &the_bs_mass_corr);
   signal_tree->Branch("bs_phi", &the_bs_phi);
   signal_tree->Branch("bs_pt", &the_bs_pt);
   signal_tree->Branch("bs_sv_prob", &the_bs_sv_prob);
   signal_tree->Branch("bs_sv_chi2", &the_bs_sv_chi2);
   signal_tree->Branch("bs_sv_ndof", &the_bs_sv_ndof);
   signal_tree->Branch("bs_vx", &the_bs_vx);
   signal_tree->Branch("bs_vy", &the_bs_vy);
   signal_tree->Branch("bs_vz", &the_bs_vz);

   signal_tree->Branch("beamspot_x", &the_beamspot_x);
   signal_tree->Branch("beamspot_y", &the_beamspot_y);
   signal_tree->Branch("beamspot_z", &the_beamspot_z);

   signal_tree->Branch("cos_theta_star_phi1", &the_cos_theta_star_phi1);
   signal_tree->Branch("cos_theta_star_phi2", &the_cos_theta_star_phi2);

   signal_tree->Branch("deltar_k1k2", &the_deltar_k1k2);
   signal_tree->Branch("deltar_k1k3", &the_deltar_k1k3);
   signal_tree->Branch("deltar_k1k4", &the_deltar_k1k4);
   signal_tree->Branch("deltar_k2k3", &the_deltar_k2k3);
   signal_tree->Branch("deltar_k2k4", &the_deltar_k2k4);
   signal_tree->Branch("deltar_k3k4", &the_deltar_k3k4);
   signal_tree->Branch("deltar_min", &the_deltar_min);
   signal_tree->Branch("deltar_max", &the_deltar_max);
   signal_tree->Branch("deltar_phi1phi2", &the_deltar_phi1phi2);

   signal_tree->Branch("k1_eta", &the_k1_eta);
   signal_tree->Branch("k1_phi", &the_k1_phi);
   signal_tree->Branch("k1_pt", &the_k1_pt);
   signal_tree->Branch("k2_eta", &the_k2_eta);
   signal_tree->Branch("k2_phi", &the_k2_phi);
   signal_tree->Branch("k2_pt", &the_k2_pt);
   signal_tree->Branch("k3_eta", &the_k3_eta);
   signal_tree->Branch("k3_phi", &the_k3_phi);
   signal_tree->Branch("k3_pt", &the_k3_pt);
   signal_tree->Branch("k4_eta", &the_k4_eta);
   signal_tree->Branch("k4_phi", &the_k4_phi);
   signal_tree->Branch("k4_pt", &the_k4_pt);

   signal_tree->Branch("k1k3_mass", &the_k1k3_mass);
   signal_tree->Branch("k1k3_pt", &the_k1k3_pt);
   signal_tree->Branch("k1k4_mass", &the_k1k4_mass);
   signal_tree->Branch("k1k4_pt", &the_k1k4_pt);
   signal_tree->Branch("k2k3_mass", &the_k2k3_mass);
   signal_tree->Branch("k2k3_pt", &the_k2k3_pt);
   signal_tree->Branch("k2k4_mass", &the_k2k4_mass);
   signal_tree->Branch("k2k4_pt", &the_k2k4_pt);

   signal_tree->Branch("phi1_idx", &phi1_idx);
   signal_tree->Branch("phi2_idx", &phi2_idx);
   signal_tree->Branch("k1_idx", &k1_idx);
   signal_tree->Branch("k2_idx", &k2_idx);
   signal_tree->Branch("k3_idx", &k3_idx);
   signal_tree->Branch("k4_idx", &k4_idx);

   signal_tree->Branch("ismatched", &ismatched);

   signal_tree->Branch("phi1_eta", &the_phi1_eta);
   signal_tree->Branch("phi1_phi", &the_phi1_phi);
   signal_tree->Branch("phi1_pt", &the_phi1_pt);
   signal_tree->Branch("phi1_mass", &the_phi1_mass);
   signal_tree->Branch("phi1_cos_theta_star_k1", &the_phi1_cos_theta_star_k1);
   signal_tree->Branch("phi1_cos_theta_star_k2", &the_phi1_cos_theta_star_k2);
   signal_tree->Branch("phi1_deltar", &the_phi1_deltar);
   signal_tree->Branch("phi1_cos2d", &the_phi1_cos2d);
   signal_tree->Branch("phi1_cxx", &the_phi1_cxx);
   signal_tree->Branch("phi1_cyx", &the_phi1_cyx);
   signal_tree->Branch("phi1_cyy", &the_phi1_cyy);
   signal_tree->Branch("phi1_czx", &the_phi1_czx);
   signal_tree->Branch("phi1_czy", &the_phi1_czy);
   signal_tree->Branch("phi1_czz", &the_phi1_czz);
   signal_tree->Branch("phi1_sv_chi2", &the_phi1_sv_chi2);
   signal_tree->Branch("phi1_sv_ndof", &the_phi1_sv_ndof);
   signal_tree->Branch("phi1_sv_prob", &the_phi1_sv_prob);
   signal_tree->Branch("phi1_lxy", &the_phi1_lxy);
   signal_tree->Branch("phi1_lxysig", &the_phi1_lxysig);
   signal_tree->Branch("phi1_mass_err", &the_phi1_mass_err);
   signal_tree->Branch("phi1_vx", &the_phi1_vx);
   signal_tree->Branch("phi1_vy", &the_phi1_vy);
   signal_tree->Branch("phi1_vz", &the_phi1_vz);
   signal_tree->Branch("phi1_is_matched", &the_phi1_is_matched);

   signal_tree->Branch("phi2_eta", &the_phi2_eta);
   signal_tree->Branch("phi2_phi", &the_phi2_phi);
   signal_tree->Branch("phi2_pt", &the_phi2_pt);
   signal_tree->Branch("phi2_mass", &the_phi2_mass);
   signal_tree->Branch("phi2_cos_theta_star_k1", &the_phi2_cos_theta_star_k1);
   signal_tree->Branch("phi2_cos_theta_star_k2", &the_phi2_cos_theta_star_k2);
   signal_tree->Branch("phi2_deltar", &the_phi2_deltar);
   signal_tree->Branch("phi2_cos2d", &the_phi2_cos2d);
   signal_tree->Branch("phi2_cxx", &the_phi2_cxx);
   signal_tree->Branch("phi2_cyx", &the_phi2_cyx);
   signal_tree->Branch("phi2_cyy", &the_phi2_cyy);
   signal_tree->Branch("phi2_czx", &the_phi2_czx);
   signal_tree->Branch("phi2_czy", &the_phi2_czy);
   signal_tree->Branch("phi2_czz", &the_phi2_czz);
   signal_tree->Branch("phi2_sv_chi2", &the_phi2_sv_chi2);
   signal_tree->Branch("phi2_sv_ndof", &the_phi2_sv_ndof);
   signal_tree->Branch("phi2_sv_prob", &the_phi2_sv_prob);
   signal_tree->Branch("phi2_lxy", &the_phi2_lxy);
   signal_tree->Branch("phi2_lxysig", &the_phi2_lxysig);
   signal_tree->Branch("phi2_mass_err", &the_phi2_mass_err);
   signal_tree->Branch("phi2_vx", &the_phi2_vx);
   signal_tree->Branch("phi2_vy", &the_phi2_vy);
   signal_tree->Branch("phi2_vz", &the_phi2_vz);
   signal_tree->Branch("phi2_is_matched", &the_phi2_is_matched);

   signal_tree->Branch("phi1_pt_times_phi2_pt", &the_phi1_pt_times_phi2_pt);

   signal_tree->Branch("k1_dcasig", &the_k1_dcasig);
   signal_tree->Branch("k1_charge", &the_k1_charge);
   signal_tree->Branch("k1_covlamlam", &the_k1_covlamlam);
   signal_tree->Branch("k1_covlamphi", &the_k1_covlamphi);
   signal_tree->Branch("k1_covphiphi", &the_k1_covphiphi);
   signal_tree->Branch("k1_covqoplam", &the_k1_covqoplam);
   signal_tree->Branch("k1_covqopphi", &the_k1_covqopphi);
   signal_tree->Branch("k1_covqopqop", &the_k1_covqopqop);
   signal_tree->Branch("k1_dxy", &the_k1_dxy);
   signal_tree->Branch("k1_dxysig", &the_k1_dxysig);
   signal_tree->Branch("k1_dz", &the_k1_dz);
   signal_tree->Branch("k1_dzsig", &the_k1_dzsig);
   signal_tree->Branch("k1_normalisedchi2", &the_k1_normalisedchi2);
   signal_tree->Branch("k1_pt_err", &the_k1_pt_err);
   signal_tree->Branch("k1_validfraction", &the_k1_validfraction);
   signal_tree->Branch("k1_vx", &the_k1_vx);
   signal_tree->Branch("k1_vy", &the_k1_vy);
   signal_tree->Branch("k1_vz", &the_k1_vz);
   signal_tree->Branch("k1_highpurityflag", &the_k1_highpurityflag);
   signal_tree->Branch("k1_islost", &the_k1_islost);
   signal_tree->Branch("k1_ispacked", &the_k1_ispacked);
   signal_tree->Branch("k1_ndof", &the_k1_ndof);
   signal_tree->Branch("k1_numberlosthits", &the_k1_numberlosthits);
   signal_tree->Branch("k1_numberpixelhits", &the_k1_numberpixelhits);
   signal_tree->Branch("k1_numberpixellayers", &the_k1_numberpixellayers);
   signal_tree->Branch("k1_numbertrackerlayers", &the_k1_numbertrackerlayers);
   signal_tree->Branch("k1_numberofvalidhits", &the_k1_numberofvalidhits);
   signal_tree->Branch("k1_numberofvalidpixelhits", &the_k1_numberofvalidpixelhits);
   signal_tree->Branch("k1_qualityindex", &the_k1_qualityindex);

   signal_tree->Branch("k2_dcasig", &the_k2_dcasig);
   signal_tree->Branch("k2_charge", &the_k2_charge);
   signal_tree->Branch("k2_covlamlam", &the_k2_covlamlam);
   signal_tree->Branch("k2_covlamphi", &the_k2_covlamphi);
   signal_tree->Branch("k2_covphiphi", &the_k2_covphiphi);
   signal_tree->Branch("k2_covqoplam", &the_k2_covqoplam);
   signal_tree->Branch("k2_covqopphi", &the_k2_covqopphi);
   signal_tree->Branch("k2_covqopqop", &the_k2_covqopqop);
   signal_tree->Branch("k2_dxy", &the_k2_dxy);
   signal_tree->Branch("k2_dxysig", &the_k2_dxysig);
   signal_tree->Branch("k2_dz", &the_k2_dz);
   signal_tree->Branch("k2_dzsig", &the_k2_dzsig);
   signal_tree->Branch("k2_normalisedchi2", &the_k2_normalisedchi2);
   signal_tree->Branch("k2_pt_err", &the_k2_pt_err);
   signal_tree->Branch("k2_validfraction", &the_k2_validfraction);
   signal_tree->Branch("k2_vx", &the_k2_vx);
   signal_tree->Branch("k2_vy", &the_k2_vy);
   signal_tree->Branch("k2_vz", &the_k2_vz);
   signal_tree->Branch("k2_highpurityflag", &the_k2_highpurityflag);
   signal_tree->Branch("k2_islost", &the_k2_islost);
   signal_tree->Branch("k2_ispacked", &the_k2_ispacked);
   signal_tree->Branch("k2_ndof", &the_k2_ndof);
   signal_tree->Branch("k2_numberlosthits", &the_k2_numberlosthits);
   signal_tree->Branch("k2_numberpixelhits", &the_k2_numberpixelhits);
   signal_tree->Branch("k2_numberpixellayers", &the_k2_numberpixellayers);
   signal_tree->Branch("k2_numbertrackerlayers", &the_k2_numbertrackerlayers);
   signal_tree->Branch("k2_numberofvalidhits", &the_k2_numberofvalidhits);
   signal_tree->Branch("k2_numberofvalidpixelhits", &the_k2_numberofvalidpixelhits);
   signal_tree->Branch("k2_qualityindex", &the_k2_qualityindex);

   signal_tree->Branch("k3_dcasig", &the_k3_dcasig);
   signal_tree->Branch("k3_charge", &the_k3_charge);
   signal_tree->Branch("k3_covlamlam", &the_k3_covlamlam);
   signal_tree->Branch("k3_covlamphi", &the_k3_covlamphi);
   signal_tree->Branch("k3_covphiphi", &the_k3_covphiphi);
   signal_tree->Branch("k3_covqoplam", &the_k3_covqoplam);
   signal_tree->Branch("k3_covqopphi", &the_k3_covqopphi);
   signal_tree->Branch("k3_covqopqop", &the_k3_covqopqop);
   signal_tree->Branch("k3_dxy", &the_k3_dxy);
   signal_tree->Branch("k3_dxysig", &the_k3_dxysig);
   signal_tree->Branch("k3_dz", &the_k3_dz);
   signal_tree->Branch("k3_dzsig", &the_k3_dzsig);
   signal_tree->Branch("k3_normalisedchi2", &the_k3_normalisedchi2);
   signal_tree->Branch("k3_pt_err", &the_k3_pt_err);
   signal_tree->Branch("k3_validfraction", &the_k3_validfraction);
   signal_tree->Branch("k3_vx", &the_k3_vx);
   signal_tree->Branch("k3_vy", &the_k3_vy);
   signal_tree->Branch("k3_vz", &the_k3_vz);
   signal_tree->Branch("k3_highpurityflag", &the_k3_highpurityflag);
   signal_tree->Branch("k3_islost", &the_k3_islost);
   signal_tree->Branch("k3_ispacked", &the_k3_ispacked);
   signal_tree->Branch("k3_ndof", &the_k3_ndof);
   signal_tree->Branch("k3_numberlosthits", &the_k3_numberlosthits);
   signal_tree->Branch("k3_numberpixelhits", &the_k3_numberpixelhits);
   signal_tree->Branch("k3_numberpixellayers", &the_k3_numberpixellayers);
   signal_tree->Branch("k3_numbertrackerlayers", &the_k3_numbertrackerlayers);
   signal_tree->Branch("k3_numberofvalidhits", &the_k3_numberofvalidhits);
   signal_tree->Branch("k3_numberofvalidpixelhits", &the_k3_numberofvalidpixelhits);
   signal_tree->Branch("k3_qualityindex", &the_k3_qualityindex);

   signal_tree->Branch("k4_dcasig", &the_k4_dcasig);
   signal_tree->Branch("k4_charge", &the_k4_charge);
   signal_tree->Branch("k4_covlamlam", &the_k4_covlamlam);
   signal_tree->Branch("k4_covlamphi", &the_k4_covlamphi);
   signal_tree->Branch("k4_covphiphi", &the_k4_covphiphi);
   signal_tree->Branch("k4_covqoplam", &the_k4_covqoplam);
   signal_tree->Branch("k4_covqopphi", &the_k4_covqopphi);
   signal_tree->Branch("k4_covqopqop", &the_k4_covqopqop);
   signal_tree->Branch("k4_dxy", &the_k4_dxy);
   signal_tree->Branch("k4_dxysig", &the_k4_dxysig);
   signal_tree->Branch("k4_dz", &the_k4_dz);
   signal_tree->Branch("k4_dzsig", &the_k4_dzsig);
   signal_tree->Branch("k4_normalisedchi2", &the_k4_normalisedchi2);
   signal_tree->Branch("k4_pt_err", &the_k4_pt_err);
   signal_tree->Branch("k4_validfraction", &the_k4_validfraction);
   signal_tree->Branch("k4_vx", &the_k4_vx);
   signal_tree->Branch("k4_vy", &the_k4_vy);
   signal_tree->Branch("k4_vz", &the_k4_vz);
   signal_tree->Branch("k4_highpurityflag", &the_k4_highpurityflag);
   signal_tree->Branch("k4_islost", &the_k4_islost);
   signal_tree->Branch("k4_ispacked", &the_k4_ispacked);
   signal_tree->Branch("k4_ndof", &the_k4_ndof);
   signal_tree->Branch("k4_numberlosthits", &the_k4_numberlosthits);
   signal_tree->Branch("k4_numberpixelhits", &the_k4_numberpixelhits);
   signal_tree->Branch("k4_numberpixellayers", &the_k4_numberpixellayers);
   signal_tree->Branch("k4_numbertrackerlayers", &the_k4_numbertrackerlayers);
   signal_tree->Branch("k4_numberofvalidhits", &the_k4_numberofvalidhits);
   signal_tree->Branch("k4_numberofvalidpixelhits", &the_k4_numberofvalidpixelhits);
   signal_tree->Branch("k4_qualityindex", &the_k4_qualityindex);
}

Bool_t BsToPhiPhiTo4KDumper::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   // for data, we skip the event in case it doesn't pass the lumi mask
   if(!isMC && lumiMask(*run, *luminosityBlock) == false) return false;

   //std::cout << "in process" << std::endl;

   //the_event = *event; 
   //the_run = *run;
   //the_lumi = *luminosityBlock; 
   //the_pu_ntrueint = *Pileup_nTrueInt;

   //the_pv_npvs = *PV_npvs;

   // number of candidates in the event
   UInt_t nCand_sig = *nBsToPhiPhiTo4K; 
   if(nCand_sig > 0){ // at least one candidate per event

      // get the index of the candidate to save 
      //vector<pair<int,float>> pair_candidx_phipt = createPairWithDesc(nCand_sig, BsToPhiPhiTo4K_phi1_pt*BsToPhiPhiTo4K_phi2_2t);
      vector<pair<int,float>> pair_candidx_phipt = createPairWithMultDesc(nCand_sig, BsToPhiPhiTo4K_phi1_pt, BsToPhiPhiTo4K_phi2_pt);
      stable_sort(pair_candidx_phipt.begin(), pair_candidx_phipt.end(), sortcansbydesc);
      UInt_t cand_idx = pair_candidx_phipt[0].first;

      // for signal MC, only keep the matched events
      //if(isSignalMC && BsToPhiPhiTo4K_isMatched[selectedCandIdx_sig] != 1) return false;

      // fill the signal_tree
      the_event = *event; 
      the_run = *run;
      the_lumi = *luminosityBlock; 
      the_pu_ntrueint = *Pileup_nTrueInt;
      the_pv_npvs = *PV_npvs;
      //the_ntriggermuon = *nTriggerMuon;

      the_bs_beta = BsToPhiPhiTo4K_Bs_beta[cand_idx];
      the_bs_eta = fabs(BsToPhiPhiTo4K_Bs_eta[cand_idx]);
      the_bs_charge = BsToPhiPhiTo4K_Bs_charge[cand_idx];
      the_bs_cos2d = BsToPhiPhiTo4K_Bs_cos2D[cand_idx];
      the_bs_cxx = BsToPhiPhiTo4K_Bs_cxx[cand_idx];
      the_bs_cyx = BsToPhiPhiTo4K_Bs_cyx[cand_idx];
      the_bs_cyy = BsToPhiPhiTo4K_Bs_cyy[cand_idx];
      the_bs_czx = BsToPhiPhiTo4K_Bs_czx[cand_idx];
      the_bs_czy = BsToPhiPhiTo4K_Bs_czy[cand_idx];
      the_bs_czz = BsToPhiPhiTo4K_Bs_czz[cand_idx];
      the_bs_lxysig = BsToPhiPhiTo4K_Bs_lxy_sig[cand_idx];
      the_bs_mass = BsToPhiPhiTo4K_Bs_mass[cand_idx];
      the_bs_mass_err = BsToPhiPhiTo4K_Bs_massErr[cand_idx];
      the_bs_mass_corr = BsToPhiPhiTo4K_Bs_mass_corr[cand_idx];
      the_bs_phi = BsToPhiPhiTo4K_Bs_phi[cand_idx];
      the_bs_pt = BsToPhiPhiTo4K_Bs_pt[cand_idx];
      the_bs_sv_prob = BsToPhiPhiTo4K_Bs_sv_prob[cand_idx];
      the_bs_sv_chi2 = BsToPhiPhiTo4K_Bs_sv_chi2[cand_idx];
      the_bs_sv_ndof = BsToPhiPhiTo4K_Bs_sv_ndof[cand_idx];
      the_bs_vx = BsToPhiPhiTo4K_Bs_vx[cand_idx];
      the_bs_vy = BsToPhiPhiTo4K_Bs_vy[cand_idx];
      the_bs_vz = BsToPhiPhiTo4K_Bs_vz[cand_idx];

      the_beamspot_x = BsToPhiPhiTo4K_beamspot_x[cand_idx];
      the_beamspot_y = BsToPhiPhiTo4K_beamspot_y[cand_idx];
      the_beamspot_z = BsToPhiPhiTo4K_beamspot_z[cand_idx];

      the_cos_theta_star_phi1 = BsToPhiPhiTo4K_cos_theta_star_phi1[cand_idx];
      the_cos_theta_star_phi2 = BsToPhiPhiTo4K_cos_theta_star_phi2[cand_idx];

      the_deltar_k1k2 = BsToPhiPhiTo4K_deltaR_k1k2[cand_idx];
      the_deltar_k1k3 = BsToPhiPhiTo4K_deltaR_k1k3[cand_idx];
      the_deltar_k1k4 = BsToPhiPhiTo4K_deltaR_k1k4[cand_idx];
      the_deltar_k2k3 = BsToPhiPhiTo4K_deltaR_k2k3[cand_idx];
      the_deltar_k2k4 = BsToPhiPhiTo4K_deltaR_k2k4[cand_idx];
      the_deltar_k3k4 = BsToPhiPhiTo4K_deltaR_k3k4[cand_idx];
      the_deltar_min = BsToPhiPhiTo4K_deltaR_min[cand_idx];
      the_deltar_max = BsToPhiPhiTo4K_deltaR_max[cand_idx];
      the_deltar_phi1phi2 = BsToPhiPhiTo4K_deltaR_phi1phi2[cand_idx];

      the_k1_eta = fabs(BsToPhiPhiTo4K_k1_eta[cand_idx]); 
      the_k1_phi = BsToPhiPhiTo4K_k1_phi[cand_idx];
      the_k1_pt = BsToPhiPhiTo4K_k1_pt[cand_idx];
      the_k2_eta = fabs(BsToPhiPhiTo4K_k2_eta[cand_idx]); 
      the_k2_phi = BsToPhiPhiTo4K_k2_phi[cand_idx];
      the_k2_pt = BsToPhiPhiTo4K_k2_pt[cand_idx];
      the_k3_eta = fabs(BsToPhiPhiTo4K_k3_eta[cand_idx]); 
      the_k3_phi = BsToPhiPhiTo4K_k3_phi[cand_idx];
      the_k3_pt = BsToPhiPhiTo4K_k3_pt[cand_idx];
      the_k4_eta = fabs(BsToPhiPhiTo4K_k4_eta[cand_idx]); 
      the_k4_phi = BsToPhiPhiTo4K_k4_phi[cand_idx];
      the_k4_pt = BsToPhiPhiTo4K_k4_pt[cand_idx];

      the_k1k3_mass = BsToPhiPhiTo4K_k1k3_mass[cand_idx];
      the_k1k3_pt = BsToPhiPhiTo4K_k1k3_pt[cand_idx];
      the_k1k4_mass = BsToPhiPhiTo4K_k1k4_mass[cand_idx];
      the_k1k4_pt = BsToPhiPhiTo4K_k1k4_pt[cand_idx];
      the_k2k3_mass = BsToPhiPhiTo4K_k2k3_mass[cand_idx];
      the_k2k3_pt = BsToPhiPhiTo4K_k2k3_pt[cand_idx];
      the_k2k4_mass = BsToPhiPhiTo4K_k2k4_mass[cand_idx];
      the_k2k4_pt = BsToPhiPhiTo4K_k2k4_pt[cand_idx];

      phi1_idx = BsToPhiPhiTo4K_phi1_idx[cand_idx];
      phi2_idx = BsToPhiPhiTo4K_phi2_idx[cand_idx];
      k1_idx = BsToPhiPhiTo4K_k1_idx[cand_idx];
      k2_idx = BsToPhiPhiTo4K_k2_idx[cand_idx];
      k3_idx = BsToPhiPhiTo4K_k3_idx[cand_idx];
      k4_idx = BsToPhiPhiTo4K_k4_idx[cand_idx];

      ismatched = BsToPhiPhiTo4K_isMatched[cand_idx];

      the_phi1_eta = fabs(BsToPhiPhiTo4K_phi1_eta[cand_idx]);
      the_phi1_phi = BsToPhiPhiTo4K_phi1_phi[cand_idx];
      the_phi1_pt = BsToPhiPhiTo4K_phi1_pt[cand_idx];
      the_phi1_mass = BsToPhiPhiTo4K_phi1_mass[cand_idx];
      the_phi1_cos_theta_star_k1 = PhiToKK_cos_theta_star_k1[phi1_idx];
      the_phi1_cos_theta_star_k2 = PhiToKK_cos_theta_star_k2[phi1_idx];
      the_phi1_deltar = PhiToKK_deltaR_postfit[phi1_idx];
      the_phi1_cos2d = PhiToKK_phi_cos2D[phi1_idx];
      the_phi1_cxx = PhiToKK_phi_cxx[phi1_idx];
      the_phi1_cyx = PhiToKK_phi_cyx[phi1_idx];
      the_phi1_cyy = PhiToKK_phi_cyy[phi1_idx];
      the_phi1_czx = PhiToKK_phi_czx[phi1_idx];
      the_phi1_czy = PhiToKK_phi_czy[phi1_idx];
      the_phi1_czz = PhiToKK_phi_czz[phi1_idx];
      the_phi1_sv_chi2 = PhiToKK_phi_sv_chi2[phi1_idx];
      the_phi1_sv_ndof = PhiToKK_phi_sv_ndof[phi1_idx];
      the_phi1_sv_prob = PhiToKK_phi_sv_prob[phi1_idx];
      the_phi1_lxy = PhiToKK_phi_lxy[phi1_idx];
      the_phi1_lxysig = PhiToKK_phi_lxy_sig[phi1_idx];
      the_phi1_mass_err = PhiToKK_phi_masserr[phi1_idx];
      the_phi1_vx = PhiToKK_phi_vx[phi1_idx];
      the_phi1_vy = PhiToKK_phi_vy[phi1_idx];
      the_phi1_vz = PhiToKK_phi_vz[phi1_idx];
      the_phi1_is_matched = PhiToKK_isMatched[phi1_idx];

      the_phi2_eta = fabs(BsToPhiPhiTo4K_phi2_eta[cand_idx]);
      the_phi2_phi = BsToPhiPhiTo4K_phi2_phi[cand_idx];
      the_phi2_pt = BsToPhiPhiTo4K_phi2_pt[cand_idx];
      the_phi2_mass = BsToPhiPhiTo4K_phi2_mass[cand_idx];
      the_phi2_cos_theta_star_k1 = PhiToKK_cos_theta_star_k1[phi2_idx];
      the_phi2_cos_theta_star_k2 = PhiToKK_cos_theta_star_k2[phi2_idx];
      the_phi2_deltar = PhiToKK_deltaR_postfit[phi2_idx];
      the_phi2_cos2d = PhiToKK_phi_cos2D[phi2_idx];
      the_phi2_cxx = PhiToKK_phi_cxx[phi2_idx];
      the_phi2_cyx = PhiToKK_phi_cyx[phi2_idx];
      the_phi2_cyy = PhiToKK_phi_cyy[phi2_idx];
      the_phi2_czx = PhiToKK_phi_czx[phi2_idx];
      the_phi2_czy = PhiToKK_phi_czy[phi2_idx];
      the_phi2_czz = PhiToKK_phi_czz[phi2_idx];
      the_phi2_sv_chi2 = PhiToKK_phi_sv_chi2[phi2_idx];
      the_phi2_sv_ndof = PhiToKK_phi_sv_ndof[phi2_idx];
      the_phi2_sv_prob = PhiToKK_phi_sv_prob[phi2_idx];
      the_phi2_lxy = PhiToKK_phi_lxy[phi2_idx];
      the_phi2_lxysig = PhiToKK_phi_lxy_sig[phi2_idx];
      the_phi2_mass_err = PhiToKK_phi_masserr[phi2_idx];
      the_phi2_vx = PhiToKK_phi_vx[phi2_idx];
      the_phi2_vy = PhiToKK_phi_vy[phi2_idx];
      the_phi2_vz = PhiToKK_phi_vz[phi2_idx];
      the_phi2_is_matched = PhiToKK_isMatched[phi2_idx];

      the_phi1_pt_times_phi2_pt = BsToPhiPhiTo4K_phi1_pt[cand_idx] * BsToPhiPhiTo4K_phi2_pt[cand_idx];

      the_k1_dcasig = PhiToKK_phi_k1_DCASig[phi1_idx];
      the_k1_charge = PhiToKK_phi_k1_charge[phi1_idx];
      the_k1_covlamlam = PhiToKK_phi_k1_covLamLam[phi1_idx];
      the_k1_covlamphi = PhiToKK_phi_k1_covLamPhi[phi1_idx];
      the_k1_covphiphi = PhiToKK_phi_k1_covPhiPhi[phi1_idx];
      the_k1_covqoplam = PhiToKK_phi_k1_covQopLam[phi1_idx];
      the_k1_covqopphi = PhiToKK_phi_k1_covQopPhi[phi1_idx];
      the_k1_covqopqop = PhiToKK_phi_k1_covQopQop[phi1_idx];
      the_k1_dxy = PhiToKK_phi_k1_dxy[phi1_idx];
      the_k1_dxysig = PhiToKK_phi_k1_dxyS[phi1_idx];
      the_k1_dz = PhiToKK_phi_k1_dz[phi1_idx];
      the_k1_dzsig = PhiToKK_phi_k1_dzS[phi1_idx];
      the_k1_normalisedchi2 = PhiToKK_phi_k1_normalisedChi2[phi1_idx];
      the_k1_pt_err = PhiToKK_phi_k1_ptErr[phi1_idx];
      the_k1_validfraction = PhiToKK_phi_k1_validFraction[phi1_idx];
      the_k1_vx = PhiToKK_phi_k1_vx[phi1_idx];
      the_k1_vy = PhiToKK_phi_k1_vy[phi1_idx];
      the_k1_vz = PhiToKK_phi_k1_vz[phi1_idx];
      the_k1_highpurityflag = PhiToKK_phi_k1_highPurityFlag[phi1_idx];
      the_k1_islost = PhiToKK_phi_k1_islost[phi1_idx];
      the_k1_ispacked = PhiToKK_phi_k1_ispacked[phi1_idx];
      the_k1_ndof = PhiToKK_phi_k1_ndof[phi1_idx];
      the_k1_numberlosthits = PhiToKK_phi_k1_numberOfLostHits[phi1_idx];
      the_k1_numberpixelhits = PhiToKK_phi_k1_numberOfPixelHits[phi1_idx];
      the_k1_numberpixellayers = PhiToKK_phi_k1_numberOfPixelLayers[phi1_idx];
      the_k1_numbertrackerlayers = PhiToKK_phi_k1_numberOfTrackerLayers[phi1_idx];
      the_k1_numberofvalidhits = PhiToKK_phi_k1_numberOfValidHits[phi1_idx];
      the_k1_numberofvalidpixelhits = PhiToKK_phi_k1_numberOfValidPixelHits[phi1_idx];
      the_k1_qualityindex = PhiToKK_phi_k1_qualityIndex[phi1_idx];

      the_k2_dcasig = PhiToKK_phi_k2_DCASig[phi1_idx];
      the_k2_charge = PhiToKK_phi_k2_charge[phi1_idx];
      the_k2_covlamlam = PhiToKK_phi_k2_covLamLam[phi1_idx];
      the_k2_covlamphi = PhiToKK_phi_k2_covLamPhi[phi1_idx];
      the_k2_covphiphi = PhiToKK_phi_k2_covPhiPhi[phi1_idx];
      the_k2_covqoplam = PhiToKK_phi_k2_covQopLam[phi1_idx];
      the_k2_covqopphi = PhiToKK_phi_k2_covQopPhi[phi1_idx];
      the_k2_covqopqop = PhiToKK_phi_k2_covQopQop[phi1_idx];
      the_k2_dxy = PhiToKK_phi_k2_dxy[phi1_idx];
      the_k2_dxysig = PhiToKK_phi_k2_dxyS[phi1_idx];
      the_k2_dz = PhiToKK_phi_k2_dz[phi1_idx];
      the_k2_dzsig = PhiToKK_phi_k2_dzS[phi1_idx];
      the_k2_normalisedchi2 = PhiToKK_phi_k2_normalisedChi2[phi1_idx];
      the_k2_pt_err = PhiToKK_phi_k2_ptErr[phi1_idx];
      the_k2_validfraction = PhiToKK_phi_k2_validFraction[phi1_idx];
      the_k2_vx = PhiToKK_phi_k2_vx[phi1_idx];
      the_k2_vy = PhiToKK_phi_k2_vy[phi1_idx];
      the_k2_vz = PhiToKK_phi_k2_vz[phi1_idx];
      the_k2_highpurityflag = PhiToKK_phi_k2_highPurityFlag[phi1_idx];
      the_k2_islost = PhiToKK_phi_k2_islost[phi1_idx];
      the_k2_ispacked = PhiToKK_phi_k2_ispacked[phi1_idx];
      the_k2_ndof = PhiToKK_phi_k2_ndof[phi1_idx];
      the_k2_numberlosthits = PhiToKK_phi_k2_numberOfLostHits[phi1_idx];
      the_k2_numberpixelhits = PhiToKK_phi_k2_numberOfPixelHits[phi1_idx];
      the_k2_numberpixellayers = PhiToKK_phi_k2_numberOfPixelLayers[phi1_idx];
      the_k2_numbertrackerlayers = PhiToKK_phi_k2_numberOfTrackerLayers[phi1_idx];
      the_k2_numberofvalidhits = PhiToKK_phi_k2_numberOfValidHits[phi1_idx];
      the_k2_numberofvalidpixelhits = PhiToKK_phi_k2_numberOfValidPixelHits[phi1_idx];
      the_k2_qualityindex = PhiToKK_phi_k2_qualityIndex[phi1_idx];

      the_k3_dcasig = PhiToKK_phi_k1_DCASig[phi2_idx];
      the_k3_charge = PhiToKK_phi_k1_charge[phi2_idx];
      the_k3_covlamlam = PhiToKK_phi_k1_covLamLam[phi2_idx];
      the_k3_covlamphi = PhiToKK_phi_k1_covLamPhi[phi2_idx];
      the_k3_covphiphi = PhiToKK_phi_k1_covPhiPhi[phi2_idx];
      the_k3_covqoplam = PhiToKK_phi_k1_covQopLam[phi2_idx];
      the_k3_covqopphi = PhiToKK_phi_k1_covQopPhi[phi2_idx];
      the_k3_covqopqop = PhiToKK_phi_k1_covQopQop[phi2_idx];
      the_k3_dxy = PhiToKK_phi_k1_dxy[phi2_idx];
      the_k3_dxysig = PhiToKK_phi_k1_dxyS[phi2_idx];
      the_k3_dz = PhiToKK_phi_k1_dz[phi2_idx];
      the_k3_dzsig = PhiToKK_phi_k1_dzS[phi2_idx];
      the_k3_normalisedchi2 = PhiToKK_phi_k1_normalisedChi2[phi2_idx];
      the_k3_pt_err = PhiToKK_phi_k1_ptErr[phi2_idx];
      the_k3_validfraction = PhiToKK_phi_k1_validFraction[phi2_idx];
      the_k3_vx = PhiToKK_phi_k1_vx[phi2_idx];
      the_k3_vy = PhiToKK_phi_k1_vy[phi2_idx];
      the_k3_vz = PhiToKK_phi_k1_vz[phi2_idx];
      the_k3_highpurityflag = PhiToKK_phi_k1_highPurityFlag[phi2_idx];
      the_k3_islost = PhiToKK_phi_k1_islost[phi2_idx];
      the_k3_ispacked = PhiToKK_phi_k1_ispacked[phi2_idx];
      the_k3_ndof = PhiToKK_phi_k1_ndof[phi2_idx];
      the_k3_numberlosthits = PhiToKK_phi_k1_numberOfLostHits[phi2_idx];
      the_k3_numberpixelhits = PhiToKK_phi_k1_numberOfPixelHits[phi2_idx];
      the_k3_numberpixellayers = PhiToKK_phi_k1_numberOfPixelLayers[phi2_idx];
      the_k3_numbertrackerlayers = PhiToKK_phi_k1_numberOfTrackerLayers[phi2_idx];
      the_k3_numberofvalidhits = PhiToKK_phi_k1_numberOfValidHits[phi2_idx];
      the_k3_numberofvalidpixelhits = PhiToKK_phi_k1_numberOfValidPixelHits[phi2_idx];
      the_k3_qualityindex = PhiToKK_phi_k1_qualityIndex[phi2_idx];

      the_k4_dcasig = PhiToKK_phi_k2_DCASig[phi2_idx];
      the_k4_charge = PhiToKK_phi_k2_charge[phi2_idx];
      the_k4_covlamlam = PhiToKK_phi_k2_covLamLam[phi2_idx];
      the_k4_covlamphi = PhiToKK_phi_k2_covLamPhi[phi2_idx];
      the_k4_covphiphi = PhiToKK_phi_k2_covPhiPhi[phi2_idx];
      the_k4_covqoplam = PhiToKK_phi_k2_covQopLam[phi2_idx];
      the_k4_covqopphi = PhiToKK_phi_k2_covQopPhi[phi2_idx];
      the_k4_covqopqop = PhiToKK_phi_k2_covQopQop[phi2_idx];
      the_k4_dxy = PhiToKK_phi_k2_dxy[phi2_idx];
      the_k4_dxysig = PhiToKK_phi_k2_dxyS[phi2_idx];
      the_k4_dz = PhiToKK_phi_k2_dz[phi2_idx];
      the_k4_dzsig = PhiToKK_phi_k2_dzS[phi2_idx];
      the_k4_normalisedchi2 = PhiToKK_phi_k2_normalisedChi2[phi2_idx];
      the_k4_pt_err = PhiToKK_phi_k2_ptErr[phi2_idx];
      the_k4_validfraction = PhiToKK_phi_k2_validFraction[phi2_idx];
      the_k4_vx = PhiToKK_phi_k2_vx[phi2_idx];
      the_k4_vy = PhiToKK_phi_k2_vy[phi2_idx];
      the_k4_vz = PhiToKK_phi_k2_vz[phi2_idx];
      the_k4_highpurityflag = PhiToKK_phi_k2_highPurityFlag[phi2_idx];
      the_k4_islost = PhiToKK_phi_k2_islost[phi2_idx];
      the_k4_ispacked = PhiToKK_phi_k2_ispacked[phi2_idx];
      the_k4_ndof = PhiToKK_phi_k2_ndof[phi2_idx];
      the_k4_numberlosthits = PhiToKK_phi_k2_numberOfLostHits[phi2_idx];
      the_k4_numberpixelhits = PhiToKK_phi_k2_numberOfPixelHits[phi2_idx];
      the_k4_numberpixellayers = PhiToKK_phi_k2_numberOfPixelLayers[phi2_idx];
      the_k4_numbertrackerlayers = PhiToKK_phi_k2_numberOfTrackerLayers[phi2_idx];
      the_k4_numberofvalidhits = PhiToKK_phi_k2_numberOfValidHits[phi2_idx];
      the_k4_numberofvalidpixelhits = PhiToKK_phi_k2_numberOfValidPixelHits[phi2_idx];
      the_k4_qualityindex = PhiToKK_phi_k2_qualityIndex[phi2_idx];

      signal_tree->Fill();


   }// end at least one candidate in the event
   return kTRUE;
}

void BsToPhiPhiTo4KDumper::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void BsToPhiPhiTo4KDumper::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   my_file->cd();

   signal_tree->Write("", TObject::kOverwrite);

   my_file->Close();

   TString option = GetOption();
   TString outFileName = option;

   std::cout << outFileName << " created" << std::endl;

   cout << "- End Bs->PhiPhi->4K Dumper -" << endl;

}
