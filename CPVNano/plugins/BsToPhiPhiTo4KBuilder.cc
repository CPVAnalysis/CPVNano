#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "../interface/helper.h"
#include <limits>
#include <algorithm>
#include "../interface/KinVtxFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"


class BsToPhiPhiTo4KBuilder : public edm::global::EDProducer<> {

public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BsToPhiPhiTo4KBuilder(const edm::ParameterSet &cfg):
    // phi candidates
    Bs_candidates_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("phis") )},
    kaon_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("phisTransientTracks") )},

    // kaon tracks
    //kaons_ {consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("kaons"))},
    //kaons_ttracks_ {consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("kaonsTransientTracks"))},

    // selection first phi candidate (phi_1)
    //k1_selection_ {cfg.getParameter<std::string>("k1_selection")},
    //k2_selection_ {cfg.getParameter<std::string>("k2_selection")},
    //pre_vtx_selection_phi1_ {cfg.getParameter<std::string>("pre_vtx_selection_phi1")},
    //post_vtx_selection_phi1_ {cfg.getParameter<std::string>("post_vtx_selection_phi1")},

    // selection second phi candidate (phi_2)
    //k3_selection_ {cfg.getParameter<std::string>("k3_selection")},
    //k4_selection_ {cfg.getParameter<std::string>("k4_selection")},
    //pre_vtx_selection_phi2_ {cfg.getParameter<std::string>("pre_vtx_selection_phi2")},
    //post_vtx_selection_phi2_ {cfg.getParameter<std::string>("post_vtx_selection_phi2")},

    // selection PV
    PV_selection_ {cfg.getParameter<std::string>("PV_selection")},

    // selection Bs candidate
    pre_vtx_selection_Bs_ {cfg.getParameter<std::string>("pre_vtx_selection_Bs")},
    post_vtx_selection_Bs_ {cfg.getParameter<std::string>("post_vtx_selection_Bs")},

    genParticles_ {consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))}, 
    isMC_ {cfg.getParameter<bool>("isMC")},

    // vertices
    vertexToken_(consumes<reco::VertexCollection> (cfg.getParameter<edm::InputTag>( "vertices" ))), 

    // beamspot
    beamspot_{consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))}{ 
       produces<pat::CompositeCandidateCollection>();
    }

  ~BsToPhiPhiTo4KBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  // kaon tracks 
  //const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  //const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> Bs_candidates_;
  const edm::EDGetTokenT<TransientTrackCollection> kaon_ttracks_;

  //// selection first phi candidate (phi_1)
  //const StringCutObjectSelector<pat::CompositeCandidate> k1_selection_; 
  //const StringCutObjectSelector<pat::CompositeCandidate> k2_selection_; 
  //const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_phi1_;
  //const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_phi1_; 

  //// selection second phi candidate (phi_2)
  //const StringCutObjectSelector<pat::CompositeCandidate> k3_selection_; 
  //const StringCutObjectSelector<pat::CompositeCandidate> k4_selection_; 
  //const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_phi2_;
  //const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_phi2_; 
    
  // selection Bs candidate
  const StringCutObjectSelector<reco::Vertex> PV_selection_;
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_Bs_;
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_Bs_; 

  const edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  const bool isMC_;

  // vertices
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  // beamspot
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BsToPhiPhiTo4KBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //TODO study strategy for tag muon
  
  // input
  //edm::Handle<pat::CompositeCandidateCollection> kaons;
  //evt.getByToken(kaons_, kaons);

  edm::Handle<pat::CompositeCandidateCollection> Bs_candidates;
  evt.getByToken(Bs_candidates_, Bs_candidates);

  //std::cout << "number oh phi candidates: " << Bs_candidates->size() << std::endl;
  
  edm::Handle<TransientTrackCollection> kaon_ttracks;
  evt.getByToken(kaon_ttracks_, kaon_ttracks);

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexToken_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();
  const reco::VertexCollection& vertices = *vertexHandle;

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //std::cout << std::endl << std::endl << "start builder" << std::endl;

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

  //if(Bs_candidates->size() < 2) continue; // there must be at least two phi candidates to proceed
  
  // loop on first phi candidate (phi1)
  for(size_t phi1_idx = 0; phi1_idx < Bs_candidates->size(); ++phi1_idx) {
    edm::Ptr<pat::CompositeCandidate> phi1_ptr(Bs_candidates, phi1_idx);
  
    //std::cout << "found first phi" << std::endl;

    // retrieve kaons
    edm::Ptr<reco::Candidate> k1_ptr = phi1_ptr->userCand("k1");
    edm::Ptr<reco::Candidate> k2_ptr = phi1_ptr->userCand("k2");
    int k1_idx = phi1_ptr->userInt("k1_idx");
    int k2_idx = phi1_ptr->userInt("k2_idx");

    //math::PtEtaPhiMLorentzVector k1_p4(
    //  k1_ptr->pt(), 
    //  k1_ptr->eta(),
    //  k1_ptr->phi(),
    //  K_MASS // impose kaon mass
    //  );

    //math::PtEtaPhiMLorentzVector k2_p4(
    //  k2_ptr->pt(), 
    //  k2_ptr->eta(),
    //  k2_ptr->phi(),
    //  K_MASS // impose kaon mass
    //  );

    // loop on second phi candidate (phi2)
    for(size_t phi2_idx = phi1_idx + 1; phi2_idx < Bs_candidates->size(); ++phi2_idx) {
      edm::Ptr<pat::CompositeCandidate> phi2_ptr(Bs_candidates, phi2_idx);

      //std::cout << "found second phi" << std::endl;
        
      // retrieve kaons
      edm::Ptr<reco::Candidate> k3_ptr = phi2_ptr->userCand("k1");
      edm::Ptr<reco::Candidate> k4_ptr = phi2_ptr->userCand("k2");
      int k3_idx = phi2_ptr->userInt("k1_idx");
      int k4_idx = phi2_ptr->userInt("k2_idx");

      if(phi2_ptr->pt() > phi1_ptr->pt()){
        std::cout << "not ordered by pt!" << std::endl;
      }

      //math::PtEtaPhiMLorentzVector k3_p4(
      //  k3_ptr->pt(), 
      //  k3_ptr->eta(),
      //  k3_ptr->phi(),
      //  K_MASS // impose kaon mass
      //  );

      //math::PtEtaPhiMLorentzVector k4_p4(
      //  k4_ptr->pt(), 
      //  k4_ptr->eta(),
      //  k4_ptr->phi(),
      //  K_MASS // impose kaon mass
      //  );

      //std::cout << k1_ptr->mass() << " " << k2_ptr->mass() << " " << k3_ptr->mass() << " " << k4_ptr->mass() << std::endl;

      if(k3_idx == k1_idx || k3_idx == k2_idx) {
        //std::cout << "I skip candidate" << std::endl;
        continue;
      }
      if(k4_idx == k1_idx || k4_idx == k2_idx){
        //std::cout << "I skip candidate" << std::endl;
        continue;
      }

      // build Bs candidate
      pat::CompositeCandidate Bs_cand;
      Bs_cand.setP4(phi1_ptr->p4() + phi2_ptr->p4());
      Bs_cand.setCharge(phi1_ptr->charge() + phi2_ptr->charge());
      
      // apply pre-fit selection on Bs candidate
      //std::cout << "before prefit selection " << std::endl;
      if( !pre_vtx_selection_Bs_(Bs_cand) ) continue;
      //std::cout << "passed prefit selection " << std::endl;

      // save indices
      Bs_cand.addUserInt("k1_idx", k1_idx);
      Bs_cand.addUserInt("k2_idx", k2_idx);
      Bs_cand.addUserInt("k3_idx", k3_idx);
      Bs_cand.addUserInt("k4_idx", k4_idx);
      Bs_cand.addUserInt("phi1_idx", phi1_idx);
      Bs_cand.addUserInt("phi2_idx", phi2_idx);

      // save candidates
      Bs_cand.addUserCand("k1", k1_ptr);
      Bs_cand.addUserCand("k2", k2_ptr);
      Bs_cand.addUserCand("k3", k3_ptr);
      Bs_cand.addUserCand("k4", k4_ptr);
      Bs_cand.addUserCand("phi1", phi1_ptr);
      Bs_cand.addUserCand("phi2", phi2_ptr);

      // fit the four kaon tracks to a common vertex
      KinVtxFitter fitter_Bs(
        {kaon_ttracks->at(k1_idx), kaon_ttracks->at(k2_idx), kaon_ttracks->at(k3_idx), kaon_ttracks->at(k4_idx)},
        {K_MASS, K_MASS, K_MASS, K_MASS}, // force kaon mass
        {K_SIGMA, K_SIGMA, K_SIGMA, K_SIGMA} 
        );

      if(!fitter_Bs.success()) continue;
      //std::cout << "passed fitter " << std::endl;

      Bs_cand.setVertex( 
        reco::Candidate::Point( 
          fitter_Bs.fitted_vtx().x(),
          fitter_Bs.fitted_vtx().y(),
          fitter_Bs.fitted_vtx().z()
          )  
        );

      float Bs_vx = fitter_Bs.fitted_vtx().x();
      float Bs_vy = fitter_Bs.fitted_vtx().y();
      float Bs_vz = fitter_Bs.fitted_vtx().z();
      Bs_cand.addUserFloat("Bs_vx", Bs_vx);
      Bs_cand.addUserFloat("Bs_vy", Bs_vy);
      Bs_cand.addUserFloat("Bs_vz", Bs_vz);

      const auto& covMatrix = fitter_Bs.fitted_vtx_uncertainty();
      Bs_cand.addUserFloat("Bs_vtx_cxx", covMatrix.cxx());
      Bs_cand.addUserFloat("Bs_vtx_cyy", covMatrix.cyy());
      Bs_cand.addUserFloat("Bs_vtx_czz", covMatrix.czz());
      Bs_cand.addUserFloat("Bs_vtx_cyx", covMatrix.cyx());
      Bs_cand.addUserFloat("Bs_vtx_czx", covMatrix.czx());
      Bs_cand.addUserFloat("Bs_vtx_czy", covMatrix.czy());

      Bs_cand.addUserFloat("Bs_sv_chi2", fitter_Bs.chi2());
      Bs_cand.addUserFloat("Bs_sv_ndof", fitter_Bs.dof()); 
      Bs_cand.addUserFloat("Bs_sv_prob", fitter_Bs.prob());

      //TODO add all quantities needed to compute boosts, angles, etc.
      //TODO add ct and time

      auto fit_p4 = fitter_Bs.fitted_p4();

      float Bs_fitted_pt = fit_p4.pt();
      Bs_cand.addUserFloat("Bs_fitted_pt", Bs_fitted_pt); 
      
      // compute ptErr from error propagation
      float Bs_fitted_px = fit_p4.px();
      float Bs_fitted_py = fit_p4.py();
      float Bs_fitted_pz = fit_p4.pz();
      float var_px = fitter_Bs.fitted_candidate().kinematicParametersError().matrix()(3,3);
      float var_py = fitter_Bs.fitted_candidate().kinematicParametersError().matrix()(4,4);
      float cov_pxpy = fitter_Bs.fitted_candidate().kinematicParametersError().matrix()(3,4);
      float Bs_fitted_ptErr = sqrt((Bs_fitted_px*Bs_fitted_px * var_px + Bs_fitted_py*Bs_fitted_py * var_py + 2*Bs_fitted_px*Bs_fitted_py*cov_pxpy) / (Bs_fitted_pt*Bs_fitted_pt)); 
      Bs_cand.addUserFloat("Bs_fitted_px", Bs_fitted_px); 
      Bs_cand.addUserFloat("Bs_fitted_py", Bs_fitted_py); 
      Bs_cand.addUserFloat("Bs_fitted_ptErr", Bs_fitted_ptErr); 

      Bs_cand.addUserFloat("Bs_fitted_eta" , fit_p4.eta());
      Bs_cand.addUserFloat("Bs_fitted_phi" , fit_p4.phi());
      Bs_cand.addUserFloat("Bs_fitted_mass", fit_p4.mass());      
      Bs_cand.addUserFloat("Bs_fitted_mass_corr", fit_p4.mass() - (fitter_Bs.daughter_p4(0) + fitter_Bs.daughter_p4(1)).mass() + 1.0195 - (fitter_Bs.daughter_p4(2) + fitter_Bs.daughter_p4(3)).mass() + 1.0195); 
      Bs_cand.addUserFloat("Bs_fitted_massErr", sqrt(fitter_Bs.fitted_candidate().kinematicParametersError().matrix()(6,6))); 
      Bs_cand.addUserFloat("Bs_charge", Bs_cand.charge());
      Bs_cand.addUserFloat("Bs_cos_theta_2D", cos_theta_2D(fitter_Bs, *beamspot, fit_p4));
      auto Bs_lxy = l_xy(fitter_Bs, *beamspot);
      Bs_cand.addUserFloat("Bs_lxy", Bs_lxy.value());
      Bs_cand.addUserFloat("Bs_lxyErr", Bs_lxy.error());
      Bs_cand.addUserFloat("Bs_lxy_sig", Bs_lxy.value() / Bs_lxy.error());

      auto Bs_lxy_corr = l_xy_corr(fitter_Bs, *beamspot, PV);
      Bs_cand.addUserFloat("Bs_lxy_corr", Bs_lxy_corr.value());
      Bs_cand.addUserFloat("Bs_lxyErr_corr", Bs_lxy_corr.error());
      Bs_cand.addUserFloat("Bs_lxy_sig_corr", Bs_lxy_corr.value() / Bs_lxy_corr.error());

      auto Bs_lxyz = l_xyz(fitter_Bs, *beamspot);
      Bs_cand.addUserFloat("Bs_lxyz", Bs_lxyz.value());
      Bs_cand.addUserFloat("Bs_lxyzErr", Bs_lxyz.error());

      auto Bs_lxyz_corr = l_xyz_corr(fitter_Bs, *beamspot, PV);
      Bs_cand.addUserFloat("Bs_lxyz_corr", Bs_lxyz_corr.value());
      Bs_cand.addUserFloat("Bs_lxyzErr_corr", Bs_lxyz_corr.error());
    
      auto k1_p4 = fitter_Bs.daughter_p4(0);
      auto k2_p4 = fitter_Bs.daughter_p4(1);
      auto k3_p4 = fitter_Bs.daughter_p4(2);
      auto k4_p4 = fitter_Bs.daughter_p4(3);
      auto phi1_p4 = k1_p4 + k2_p4;
      auto phi2_p4 = k3_p4 + k4_p4;

      //TODO necessary to store this info, cannot be retrieved from PhiToKK collection?
      //TODO make sure that the index gymnastic with PhiToKK is sound
      // phi candidates
      std::vector<std::string> phi_names{"phi1", "phi2"};
      for (size_t i = 0; i < phi_names.size(); i++){
        int idx1 = -99;
        int idx2 = -99;
        if(phi_names[i] == "phi1"){
          idx1 = 0;
          idx2 = 1;
        }
        else if(phi_names[i] == "phi2"){
          idx1 = 2;
          idx2 = 3;
        }
        Bs_cand.addUserFloat(phi_names[i] + "_fitted_mass", (fitter_Bs.daughter_p4(idx1) + fitter_Bs.daughter_p4(idx2)).mass());
        Bs_cand.addUserFloat(phi_names[i] + "_fitted_pt", (fitter_Bs.daughter_p4(idx1) + fitter_Bs.daughter_p4(idx2)).pt());
        Bs_cand.addUserFloat(phi_names[i] + "_fitted_eta", (fitter_Bs.daughter_p4(idx1) + fitter_Bs.daughter_p4(idx2)).eta());
        Bs_cand.addUserFloat(phi_names[i] + "_fitted_phi", (fitter_Bs.daughter_p4(idx1) + fitter_Bs.daughter_p4(idx2)).phi());
      }
        // TODO get other quantities (sv_prob, cos2d, lxy, etc)

      // compute deltaR
      Bs_cand.addUserFloat("deltaR_phi1phi2", reco::deltaR(phi1_p4, phi2_p4));
      Bs_cand.addUserFloat("deltaR_k1k3", reco::deltaR(k1_p4, k3_p4));
      Bs_cand.addUserFloat("deltaR_k1k4", reco::deltaR(k1_p4, k4_p4));
      Bs_cand.addUserFloat("deltaR_k2k3", reco::deltaR(k2_p4, k3_p4));
      Bs_cand.addUserFloat("deltaR_k2k4", reco::deltaR(k2_p4, k4_p4));

      Bs_cand.addUserFloat("deltaR_k1k2", reco::deltaR(k1_p4, k2_p4));
      Bs_cand.addUserFloat("deltaR_k3k4", reco::deltaR(k3_p4, k4_p4));

      auto dr_info = min_max_dr({k1_ptr, k2_ptr, k3_ptr, k4_ptr});
      Bs_cand.addUserFloat("deltaR_min", dr_info.first);
      Bs_cand.addUserFloat("deltaR_max", dr_info.second);


      // apply pots-fit selection
      if(!post_vtx_selection_Bs_(Bs_cand)) continue; //TODO move to after definition of user floats?

      // compute invariant masses
      Bs_cand.addUserFloat("k1k3_mass", (k1_p4 + k3_p4).mass());
      Bs_cand.addUserFloat("k1k4_mass", (k1_p4 + k4_p4).mass());
      Bs_cand.addUserFloat("k2k3_mass", (k2_p4 + k3_p4).mass());
      Bs_cand.addUserFloat("k2k4_mass", (k2_p4 + k4_p4).mass());

      Bs_cand.addUserFloat("k1k3_pt", (k1_p4 + k3_p4).pt());
      Bs_cand.addUserFloat("k1k4_pt", (k1_p4 + k4_p4).pt());
      Bs_cand.addUserFloat("k2k3_pt", (k2_p4 + k3_p4).pt());
      Bs_cand.addUserFloat("k2k4_pt", (k2_p4 + k4_p4).pt());


      //TODO necessary to store this info, cannot be retrieved from PhiToKK collection?
      // kaons
      std::vector<std::string> kaon_names{"k1", "k2", "k3", "k4"};
      
      for (size_t i = 0; i < kaon_names.size(); i++){
	      Bs_cand.addUserFloat(kaon_names[i] + "_fitted_mass" ,fitter_Bs.daughter_p4(i).mass());
	      Bs_cand.addUserFloat(kaon_names[i] + "_fitted_pt" ,fitter_Bs.daughter_p4(i).pt());
        Bs_cand.addUserFloat(kaon_names[i] + "_fitted_eta",fitter_Bs.daughter_p4(i).eta());
        Bs_cand.addUserFloat(kaon_names[i] + "_fitted_phi",fitter_Bs.daughter_p4(i).phi());
      }
      
      // computation of cos(theta*), 
      // (angle between the Bs momentum direction in the lab frame and the daughter's momentum direction in the center of mass frame)
      float mass_Bs_fitted = fitter_Bs.fitted_candidate().mass();
      float mass_phi1 = (fitter_Bs.daughter_p4(0) + fitter_Bs.daughter_p4(1)).mass(); 
      float mass_phi2 = (fitter_Bs.daughter_p4(2) + fitter_Bs.daughter_p4(3)).mass(); 

      float energy_Bs_fitted_lab = sqrt(pow(fitter_Bs.fitted_candidate().globalMomentum().x(), 2) + pow(fitter_Bs.fitted_candidate().globalMomentum().y(), 2) + pow(fitter_Bs.fitted_candidate().globalMomentum().z(), 2) + pow(fitter_Bs.fitted_candidate().mass(), 2));
      //float energy_Bs_fitted_cm = mass_Bs_fitted;
      float energy_phi1_fitted_lab = (fitter_Bs.daughter_p4(0) + fitter_Bs.daughter_p4(1)).energy();
      float energy_phi1_fitted_cm = (pow(mass_phi1, 2) - pow(mass_phi2, 2) + pow(mass_Bs_fitted, 2)) / (2. * mass_Bs_fitted); 
      float energy_phi2_fitted_lab = (fitter_Bs.daughter_p4(1) + fitter_Bs.daughter_p4(2)).energy();
      float energy_phi2_fitted_cm = (pow(mass_phi2, 2) - pow(mass_phi1, 2) + pow(mass_Bs_fitted, 2)) / (2. * mass_Bs_fitted); 

      float momentum_phi1_fitted_cm = sqrt(pow(energy_phi1_fitted_cm, 2) - pow(mass_phi1, 2));
      float momentum_phi2_fitted_cm = sqrt(pow(energy_phi2_fitted_cm, 2) - pow(mass_phi2, 2));

      float beta = sqrt(pow(energy_Bs_fitted_lab, 2) - pow(mass_Bs_fitted, 2)) / energy_Bs_fitted_lab;
      float gamma = 1. / sqrt(1 - pow(beta, 2));

      float cos_theta_star_phi1 = (1. / (beta * momentum_phi1_fitted_cm)) * (energy_phi1_fitted_lab / gamma - energy_phi1_fitted_cm);
      float cos_theta_star_phi2 = (1. / (beta * momentum_phi2_fitted_cm)) * (energy_phi2_fitted_lab / gamma - energy_phi2_fitted_cm);

      Bs_cand.addUserFloat("Bs_beta", beta);
      Bs_cand.addUserFloat("Bs_gamma", gamma);
      Bs_cand.addUserFloat("cos_theta_star_phi1", cos_theta_star_phi1);
      Bs_cand.addUserFloat("cos_theta_star_phi2", cos_theta_star_phi2);

      // beamspot position
      float beamspot_x = beamspot->x0();
      float beamspot_y = beamspot->y0();
      float beamspot_z = beamspot->z0();
      Bs_cand.addUserFloat("beamspot_x", beamspot_x);
      Bs_cand.addUserFloat("beamspot_y", beamspot_y);
      Bs_cand.addUserFloat("beamspot_z", beamspot_z);

      // --  compute ct (transverse direction) -- // 
      const float mass_Bs_PDG = 5.36691; // GeV 

      GlobalPoint Bs_vtx = fitter_Bs.fitted_vtx();
      auto bs_pos = (*beamspot).position(Bs_vtx.z());
      auto bs_posbspv = (*beamspot).position(PV.z());

      // find PV with the smallest distance to Bs
      // compute distance between Bs momentum and PV as dist = |d x p| / |p|
      const reco::Vertex& the_PV = PV; // initialise it to the first primary vertex
      float the_dist = -99;
      double p_mag = TMath::Sqrt(Bs_fitted_px * Bs_fitted_px + Bs_fitted_py * Bs_fitted_py + Bs_fitted_pz * Bs_fitted_pz);

      for(const reco::Vertex& vertex: vertices){
        // apply PV selection 
        if(!PV_selection_(vertex)) continue;

        // compute distance
        double dx = static_cast<double>(vertex.position().x() - Bs_vtx.x());
        double dy = static_cast<double>(vertex.position().y() - Bs_vtx.y());
        double dz = static_cast<double>(vertex.position().z() - Bs_vtx.z());

        double cross_x = dy * Bs_fitted_pz - dz * Bs_fitted_py;
        double cross_y = dz * Bs_fitted_px - dx * Bs_fitted_pz;
        double cross_z = dx * Bs_fitted_py - dy * Bs_fitted_px;

        double cross_mag = TMath::Sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);

        double dist = cross_mag / p_mag;

        // get the PV with the smallest distance
        if(the_dist == -99 || dist < the_dist){
          the_dist = dist;
          auto the_PV = vertex;
        }
      }

      Bs_cand.addUserFloat("the_PV_chi2", the_PV.chi2());
      Bs_cand.addUserFloat("the_PV_ndof", the_PV.ndof());
      Bs_cand.addUserFloat("the_PV_covXX", the_PV.covariance(0, 0));
      Bs_cand.addUserFloat("the_PV_covYY", the_PV.covariance(1, 1));
      Bs_cand.addUserFloat("the_PV_covZZ", the_PV.covariance(2, 2));
      Bs_cand.addUserFloat("the_PV_covXY", the_PV.covariance(0, 1));
      Bs_cand.addUserFloat("the_PV_covXZ", the_PV.covariance(0, 2));
      Bs_cand.addUserFloat("the_PV_covYZ", the_PV.covariance(1, 2));
      Bs_cand.addUserFloat("the_PV_x", the_PV.position().x());
      Bs_cand.addUserFloat("the_PV_y", the_PV.position().y());
      Bs_cand.addUserFloat("the_PV_z", the_PV.position().z());

      float Lx = Bs_vtx.x() - beamspot_x;
      float Ly = Bs_vtx.y() - beamspot_y;

      float Lx_posbsz = Bs_vtx.x() - bs_pos.x();
      float Ly_posbsz = Bs_vtx.y() - bs_pos.y();

      float Lx_posbspv = Bs_vtx.x() - bs_posbspv.x();
      float Ly_posbspv = Bs_vtx.y() - bs_posbspv.y();

      float Lx_posthepv = Bs_vtx.x() - the_PV.position().x();
      float Ly_posthepv = Bs_vtx.y() - the_PV.position().y();

      float ct_2D_cm = mass_Bs_PDG * (Lx * Bs_fitted_px + Ly * Bs_fitted_py) / (Bs_fitted_px * Bs_fitted_px + Bs_fitted_py * Bs_fitted_py);
      float ct_2D_cm_posbsz = mass_Bs_PDG * (Lx_posbsz * Bs_fitted_px + Ly_posbsz * Bs_fitted_py) / (Bs_fitted_px * Bs_fitted_px + Bs_fitted_py * Bs_fitted_py);
      float ct_2D_cm_posbspv = mass_Bs_PDG * (Lx_posbspv * Bs_fitted_px + Ly_posbspv * Bs_fitted_py) / (Bs_fitted_px * Bs_fitted_px + Bs_fitted_py * Bs_fitted_py);
      float ct_2D_cm_posthepv = mass_Bs_PDG * (Lx_posthepv * Bs_fitted_px + Ly_posthepv * Bs_fitted_py) / (Bs_fitted_px * Bs_fitted_px + Bs_fitted_py * Bs_fitted_py);
      float ct_3D_cm = Bs_lxy.value() / (beta * gamma);  // take lxy computed at zpos

      Bs_cand.addUserFloat("Bs_lx", Lx);
      Bs_cand.addUserFloat("Bs_ly", Ly);
      Bs_cand.addUserFloat("Bs_lx_posbsz", Lx_posbsz);
      Bs_cand.addUserFloat("Bs_ly_posbsz", Ly_posbsz);
      Bs_cand.addUserFloat("Bs_lx_posbspv", Lx_posbspv);
      Bs_cand.addUserFloat("Bs_ly_posbspv", Ly_posbspv);
      Bs_cand.addUserFloat("Bs_lx_posthepv", Lx_posthepv);
      Bs_cand.addUserFloat("Bs_ly_posthepv", Ly_posthepv);
      Bs_cand.addUserFloat("Bs_ct_2D_cm", ct_2D_cm);
      Bs_cand.addUserFloat("Bs_ct_2D_cm_posbsz", ct_2D_cm_posbsz);
      Bs_cand.addUserFloat("Bs_ct_2D_cm_posbspv", ct_2D_cm_posbspv);
      Bs_cand.addUserFloat("Bs_ct_2D_cm_posthepv", ct_2D_cm_posthepv);
      Bs_cand.addUserFloat("Bs_ct_3D_cm", ct_3D_cm);
       
      // gen-matching (for MC only)
      int isMatched = 0;
      //int k1_isMatched(-1), k2_isMatched(-1), k3_isMatched(-1), k4_isMatched(-1);
      int k1_genIdx(-1), k2_genIdx(-1), k3_genIdx(-1), k4_genIdx(-1);
      int genKaon1Mother_genPdgId(-1), genKaon2Mother_genPdgId(-1), genKaon3Mother_genPdgId(-1), genKaon4Mother_genPdgId(-1);
      int genKaon1GrandMother_genPdgId(-1), genKaon2GrandMother_genPdgId(-1), genKaon3GrandMother_genPdgId(-1), genKaon4GrandMother_genPdgId(-1);

      if(isMC_ == true){
        // pdgId of the gen particle to which the final-state particles are matched
        int k1_genPdgId = phi1_ptr->userInt("k1_mcMatch");
        int k2_genPdgId = phi1_ptr->userInt("k2_mcMatch");
        int k3_genPdgId = phi2_ptr->userInt("k1_mcMatch");
        int k4_genPdgId = phi2_ptr->userInt("k2_mcMatch");

        // index of the gen particle to which the final-state particles are matched
        k1_genIdx = phi1_ptr->userInt("k1_mcMatchIndex"); 
        k2_genIdx = phi1_ptr->userInt("k2_mcMatchIndex"); 
        k3_genIdx = phi2_ptr->userInt("k1_mcMatchIndex"); 
        k4_genIdx = phi2_ptr->userInt("k2_mcMatchIndex"); 

        //std::cout << k1_genIdx << " " << k2_genIdx << " " << k3_genIdx << " " << k4_genIdx << std::endl;
        
        //if(k1_genIdx != -1 && k1_genIdx != k2_genIdx && k1_genIdx != k3_genIdx && k1_genIdx != k4_genIdx){
        if(k1_genIdx != -1 && k2_genIdx != -1 && k3_genIdx != -1 && k4_genIdx != -1){
          //TODO check indices are all different

          // getting the associated gen particles
          edm::Ptr<reco::GenParticle> genKaon1_ptr(genParticles, k1_genIdx);
          edm::Ptr<reco::GenParticle> genKaon2_ptr(genParticles, k2_genIdx);
          edm::Ptr<reco::GenParticle> genKaon3_ptr(genParticles, k3_genIdx);
          edm::Ptr<reco::GenParticle> genKaon4_ptr(genParticles, k4_genIdx);

          // index of the associated mother particle
          int genKaon1Mother_genIdx = -1;
          int genKaon2Mother_genIdx = -1;
          int genKaon3Mother_genIdx = -1;
          int genKaon4Mother_genIdx = -1;
          if(genKaon1_ptr->numberOfMothers()>0) genKaon1Mother_genIdx = genKaon1_ptr->motherRef(0).key();
          if(genKaon2_ptr->numberOfMothers()>0) genKaon2Mother_genIdx = genKaon2_ptr->motherRef(0).key();
          if(genKaon3_ptr->numberOfMothers()>0) genKaon3Mother_genIdx = genKaon3_ptr->motherRef(0).key();
          if(genKaon4_ptr->numberOfMothers()>0) genKaon4Mother_genIdx = genKaon4_ptr->motherRef(0).key();

          // getting the mother particle
          edm::Ptr<reco::GenParticle> genKaon1Mother_ptr(genParticles, genKaon1Mother_genIdx);
          edm::Ptr<reco::GenParticle> genKaon2Mother_ptr(genParticles, genKaon2Mother_genIdx);
          edm::Ptr<reco::GenParticle> genKaon3Mother_ptr(genParticles, genKaon3Mother_genIdx);
          edm::Ptr<reco::GenParticle> genKaon4Mother_ptr(genParticles, genKaon4Mother_genIdx);

          // pdgId of the mother particles
          genKaon1Mother_genPdgId = genKaon1Mother_ptr->pdgId();
          genKaon2Mother_genPdgId = genKaon2Mother_ptr->pdgId();
          genKaon3Mother_genPdgId = genKaon3Mother_ptr->pdgId();
          genKaon4Mother_genPdgId = genKaon4Mother_ptr->pdgId();
          
          // index of the grandmother particle
          int genKaon1GrandMother_genIdx = -1;
          int genKaon2GrandMother_genIdx = -1;
          int genKaon3GrandMother_genIdx = -1;
          int genKaon4GrandMother_genIdx = -1;
          if(genKaon1Mother_ptr->numberOfMothers()>0) genKaon1GrandMother_genIdx = genKaon1Mother_ptr->motherRef(0).key();
          if(genKaon2Mother_ptr->numberOfMothers()>0) genKaon2GrandMother_genIdx = genKaon2Mother_ptr->motherRef(0).key();
          if(genKaon3Mother_ptr->numberOfMothers()>0) genKaon3GrandMother_genIdx = genKaon3Mother_ptr->motherRef(0).key();
          if(genKaon4Mother_ptr->numberOfMothers()>0) genKaon4GrandMother_genIdx = genKaon4Mother_ptr->motherRef(0).key();
          
          // getting the grandmother particle
          edm::Ptr<reco::GenParticle> genKaon1GrandMother_ptr(genParticles, genKaon1GrandMother_genIdx);
          edm::Ptr<reco::GenParticle> genKaon2GrandMother_ptr(genParticles, genKaon2GrandMother_genIdx);
          edm::Ptr<reco::GenParticle> genKaon3GrandMother_ptr(genParticles, genKaon3GrandMother_genIdx);
          edm::Ptr<reco::GenParticle> genKaon4GrandMother_ptr(genParticles, genKaon4GrandMother_genIdx);

          // pdgId of the grandmother particles
          genKaon1GrandMother_genPdgId = genKaon1GrandMother_ptr->pdgId();
          genKaon2GrandMother_genPdgId = genKaon2GrandMother_ptr->pdgId();
          genKaon3GrandMother_genPdgId = genKaon3GrandMother_ptr->pdgId();
          genKaon4GrandMother_genPdgId = genKaon4GrandMother_ptr->pdgId();

          //std::cout << "k1: " << k1_genPdgId << " mother: " << genKaon1Mother_genPdgId << " grandmother: " << genKaon1GrandMother_genPdgId << std::endl;

          // candidate matching
          // match to kaon
          if(fabs(k1_genPdgId) == 321 && fabs(k2_genPdgId) == 321 && fabs(k3_genPdgId) == 321  && fabs(k4_genPdgId) == 321){
            // match to phi meson
            if(fabs(genKaon1Mother_genPdgId) == 333 && fabs(genKaon2Mother_genPdgId) == 333 && fabs(genKaon3Mother_genPdgId) == 333 && fabs(genKaon4Mother_genPdgId) == 333){
              // phi mesons must be different
              if(genKaon1Mother_genIdx == genKaon2Mother_genIdx && genKaon3Mother_genIdx == genKaon4Mother_genIdx && genKaon1Mother_genIdx != genKaon3Mother_genIdx){
                // match to Bs meson
                if(fabs(genKaon1GrandMother_genPdgId) == 531 && fabs(genKaon2GrandMother_genPdgId) == 531 && fabs(genKaon3GrandMother_genPdgId) == 531 && fabs(genKaon4GrandMother_genPdgId) == 531){
                  // common ancestor Bs meson
                  if(genKaon1GrandMother_genIdx == genKaon2GrandMother_genIdx && genKaon1GrandMother_genIdx == genKaon3GrandMother_genIdx && genKaon1GrandMother_genIdx == genKaon4GrandMother_genIdx){
                    isMatched = 1;
                  }
                }
              }
            }
          }
        }
      }

      Bs_cand.addUserInt("isMatched", isMatched);

      // save candidate
      ret_value->push_back(Bs_cand);
    }
  }
  
  evt.put(std::move(ret_value));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BsToPhiPhiTo4KBuilder);
