#define TrackMatcher_cxx
// The class definition in TrackMatcher.h has been generated automatically
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
// root> T->Process("TrackMatcher.C")
// root> T->Process("TrackMatcher.C","some options")
// root> T->Process("TrackMatcher.C+")
//


#include "TrackMatcher.h"
#include <TMath.h>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "utils.C"


using namespace std;


void TrackMatcher::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  cout << " --------------------------" << endl;
  cout << "       Track Matcher       " << endl;
  cout << " --------------------------" << endl;
}


void TrackMatcher::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  TString outFileName = option;

  if(outFileName.Contains("isMC")){
    isMC = true;
    outFileName.Resize(outFileName.Length()-5);
  }
  else isMC = false;

  // check if outputfile exists
  if(gSystem->AccessPathName(outFileName)){
    my_file = new TFile(outFileName, "RECREATE");  
  }
  else{
    my_file = new TFile(outFileName, "UPDATE");  
  }

  hist_deltaR_allmuons = new TH1F("hist_deltaR_allmuons", "DeltaR distribution_allmuons", 30, 0, 10);
  hist_deltaPtRel_allmuons = new TH1F("hist_deltaPtRel_allmuons", "DeltaPtRel distribution_allmuons", 30, 0, 10);
  //hist_2d_deltaR_deltaPtRel_matchingeff = new TH2F("hist_2d_deltaR_deltaPtRel_matchingeff", "2d matching efficiency map", 6, 0, 6, 6, 0, 6);
  //hist_2d_deltaR_deltaPtRel_matchingeff = new TH2F("hist_2d_deltaR_deltaPtRel_matchingeff", "2d matching efficiency map", 43, 0, 43, 27, 0, 27);
  //hist_2d_deltaR_deltaPtRel_matchingeff = new TH2F("hist_2d_deltaR_deltaPtRel_matchingeff", "2d matching efficiency map", 12, 0, 12, 11, 0, 11);
  hist_2d_deltaR_deltaPtRel_matchingeff = new TH2F("hist_2d_deltaR_deltaPtRel_matchingeff", "2d matching efficiency map", 6, 0, 6, 11, 0, 11);
  //hist_2d_deltaR_deltaPtRel_count = new TH2F("hist_2d_deltaR_deltaPtRel_count", "2d count map", 6, 0, 6, 6, 0, 6);
  hist_ismatched_tomuon = new TH1F("hist_ismatched_tomuon", "matching efficiency", 20, 0, 20);
  hist_count = new TH1F("hist_count", "count number of events", 20, 0, 1e9);


}


Bool_t TrackMatcher::Process(Long64_t entry)
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
  //if(entry > 100000) return false;
  //cout << endl << "--- Entry " << entry << " ---" << endl;

  // for data, we skip the event in case it doesn't pass the lumi mask
  if(!isMC && lumiMask(*run, *luminosityBlock) == false) return false;

  if(*nBToMuMuPi > 0){ // at least one candidate per event

    //std::cout << std::endl << "pion pt " << the_sig_pi_pt << " eta " << the_sig_pi_eta << " phi " << the_sig_pi_phi << std::endl;
    //for(unsigned int i(0); i<*nMuon; ++i){
    //  std::cout << "muon " << i << " pt " << Muon_pt[i] << " eta " << Muon_eta[i] << " phi " << Muon_phi[i] << std::endl;
    //}

    //TODO
    // make a more inclusive signal sample
    // study inclusive and SS
    // try to make ratio between efficiencies 
    // adapt setrangeuser?

    // select the best mumupi candidate
    vector<pair<int,float>> pair_candIdx_desc_cos2d_sig = createPairWithDesc(*nBToMuMuPi, BToMuMuPi_hnl_cos2D);
    stable_sort(pair_candIdx_desc_cos2d_sig.begin(), pair_candIdx_desc_cos2d_sig.end(), sortcansbydesc);

    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sig, BToMuMuPi_hnl_charge);
    stable_sort(pair_candIdx_desc_cos2d_sign_sig.begin(), pair_candIdx_desc_cos2d_sign_sig.end(), sortcansbydesc_opp);

    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_sig, BToMuMuPi_isMatched);
    stable_sort(pair_candIdx_desc_cos2d_sign_matched_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_sig.end(), sortcansbydesc);

    vector<pair<int,float>> pair_candIdx_desc_cos2d_sign_matched_muon_sig = updatePairWithDesc(pair_candIdx_desc_cos2d_sign_sig, BToMuMuPi_sel_mu_idx, Muon_isDSAMuon);
    stable_sort(pair_candIdx_desc_cos2d_sign_matched_muon_sig.begin(), pair_candIdx_desc_cos2d_sign_matched_muon_sig.end(), sortcansbydesc_opp);

    UInt_t selectedCandIdx_sig = pair_candIdx_desc_cos2d_sign_matched_sig[0].first;

    // get the information on the track (non-fitted)
    float the_trk_pt = BToMuMuPi_pi_pt[selectedCandIdx_sig];
    float the_trk_eta = BToMuMuPi_pi_eta[selectedCandIdx_sig];
    float the_trk_phi = BToMuMuPi_pi_phi[selectedCandIdx_sig];

    if(isMC && BToMuMuPi_isMatched[selectedCandIdx_sig]!=1) return false;
    //if(Muon_charge[BToMuMuPi_trg_mu_idx[selectedCandIdx_sig]] != Muon_charge[BToMuMuPi_sel_mu_idx[selectedCandIdx_sig]]) return false;

    // lxy0to1_SS
    int trgmu_idx = BToMuMuPi_trg_mu_idx[selectedCandIdx_sig];
    int mu_idx = BToMuMuPi_sel_mu_idx[selectedCandIdx_sig];
    if(Muon_softId[trgmu_idx]==1 && 
       Muon_looseId[mu_idx] && 
       (( Muon_charge[trgmu_idx]!= Muon_charge[mu_idx] && fabs(BToMuMuPi_trgmu_mu_mass[selectedCandIdx_sig]-3.097)>0.15 && fabs(BToMuMuPi_trgmu_mu_mass[selectedCandIdx_sig]-3.686)>0.08 && fabs(BToMuMuPi_trgmu_mu_mass[selectedCandIdx_sig]-1.019)>0.01) || (Muon_charge[trgmu_idx]== Muon_charge[mu_idx]))
       && Muon_charge[trgmu_idx] == Muon_charge[mu_idx] &&
       ((BToMuMuPi_hnl_mass[selectedCandIdx_sig] > 1.6 && BToMuMuPi_hnl_mass[selectedCandIdx_sig] < 1.9) || (BToMuMuPi_hnl_mass[selectedCandIdx_sig] > 3 && BToMuMuPi_hnl_mass[selectedCandIdx_sig] < 3.2)) &&
       //BToMuMuPi_sv_lxy[selectedCandIdx_sig] < 1 &&
       BToMuMuPi_hnl_charge[selectedCandIdx_sig] == 0 //&&
       //BToMuMuPi_fit_pi_pt[selectedCandIdx_sig] > 1.2 &&
       //BToMuMuPi_sv_lxy_sig[selectedCandIdx_sig] > 100 &&
       //fabs(Muon_dxyS[mu_idx]) > 12 &&
       //fabs(BToMuMuPi_pi_dxyS[selectedCandIdx_sig]) > 25
       )
       {   
      
      std::vector<pair<int, std::array<float, 2>>> pairs_muonIdx_deltaInfo;

      float deltaR_max = 0.3;
      float deltaPtRel_max = 10.;

      //vector<float> deltaR_cond {0.03, 0.08, 0.1, 0.3, 0.8, 1};
      //vector<float> deltaR_cond {0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35, 0.37, 0.39, 0.41, 0.43, 0.45, 0.47, 0.49, 0.51, 0.53, 0.55, 0.57, 0.59, 0.61, 0.63, 0.65, 0.67, 0.69, 0.71, 0.75, 0.78, 0.8, 0.85, 0.9, 0.95, 1.0};
      //vector<float> deltaR_cond {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};
      vector<float> deltaR_cond {0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
      //vector<float> deltaPtRel_cond {0.03, 0.1, 0.3, 0.8, 1, 1.5};
      //vector<float> deltaPtRel_cond {0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.5, 2, 2.5, 3, 3.5, 4};
      vector<float> deltaPtRel_cond {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

      //vector<TString> deltaR_label {"0.03", "0.08", "0.1", "0.3", "0.8", "1.0"};
      //vector<TString> deltaR_label {"0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6"};
      vector<TString> deltaR_label {"0.05", "0.1", "0.15", "0.2", "0.25", "0.3"};
      //vector<TString> deltaPtRel_label {"0.03", "0.1", "0.3", "0.8", "1.0", "1.5"};
      vector<TString> deltaPtRel_label {"0.05", "0.1", "0.2", "0.4", "0.3", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"};

      gStyle->SetOptStat(0);

      for(unsigned int iR(0); iR<deltaR_cond.size(); ++iR){
        for(unsigned int iPt(0); iPt<deltaPtRel_cond.size(); ++iPt){
          bool ismatched_tomuon = 0;

          // loop on the slimmed muon
          for(unsigned int iMuon(0); iMuon<*nMuon; ++iMuon){
            // exclude signal muons
            if(iMuon == trgmu_idx || iMuon == mu_idx) continue;
            // compute the deltaR and deltaPtRel between the given track and the muon
            float deltaR_track_muon = reco::deltaR(Muon_eta[iMuon], Muon_phi[iMuon], the_trk_eta, the_trk_phi);
            float deltaPtRel = fabs(Muon_pt[iMuon] - the_trk_pt) / the_trk_pt;

            hist_deltaR_allmuons->Fill(deltaR_track_muon);
            hist_deltaPtRel_allmuons->Fill(deltaPtRel);

            // in any muon satisfies the two conditions, we consider the track to be matched
            //if(deltaR_track_muon < deltaR_max && deltaPtRel < deltaPtRel_max){
            if(deltaR_track_muon < deltaR_cond[iR] && deltaPtRel < deltaPtRel_cond[iPt]){
              ismatched_tomuon = 1; 
              //pair<int, std::array<float, 2>> pairs_muonIdx_deltaInfo_tmp;
              //pairs_muonIdx_deltaInfo_tmp.first = iMuon;
              //pairs_muonIdx_deltaInfo_tmp.second[0] = deltaPtRel;
              //pairs_muonIdx_deltaInfo_tmp.second[1] = deltaR_track_muon;
              //pairs_muonIdx_deltaInfo.push_back(pairs_muonIdx_deltaInfo_tmp);
            }
          }
          //std::cout << deltaR_cond[iR] << " " << deltaPtRel_cond[iPt] << " "  << ismatched_tomuon << std::endl;
          hist_2d_deltaR_deltaPtRel_matchingeff->Fill(deltaR_label[iR], deltaPtRel_label[iPt], ismatched_tomuon);

          //hist_2d_deltaR_deltaPtRel_matchingeff->Fill(deltaR_cond[iR], deltaPtRel_cond[iPt], ismatched_tomuon);
          //hist_2d_deltaR_deltaPtRel_matchingeff->Fill(iR, iPt, ismatched_tomuon);
          //hist_2d_deltaR_deltaPtRel_count->Fill(iR, iPt, 1.);
        }
      }

      hist_2d_deltaR_deltaPtRel_matchingeff->GetXaxis()->SetTitle("max #DeltaR");
      hist_2d_deltaR_deltaPtRel_matchingeff->GetYaxis()->SetTitle("max #DeltaPtRel");
      hist_2d_deltaR_deltaPtRel_matchingeff->GetZaxis()->SetTitle("Matching efficiency");
      hist_2d_deltaR_deltaPtRel_matchingeff->GetZaxis()->SetRangeUser(0, 1);

      hist_count->Fill(1.);

      //hist_ismatched_tomuon->Fill(ismatched_tomuon);

      // make 2d scan in deltaR, deltaPt with matching efficiency in z axis
      // do this for both data and mc
      // For the data, we want the matching efficiency to be the largest, in mc the smallest (as otherwise events would be rejected)
      // Compare histograms and find sweet spot
      // For the signal, do not forget to ask the candidate to be matched
      
    } // end category selection
  }// end at least one candidate in the event

  return kTRUE;
}


void TrackMatcher::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}


void TrackMatcher::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  my_file->cd();
  my_file->Write();
  my_file->Close();

  cout << "- End TrackMatcher -" << endl;
}

