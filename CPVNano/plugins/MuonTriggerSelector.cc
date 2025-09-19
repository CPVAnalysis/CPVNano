// class to produce 2 pat::MuonCollections
// one matched to the Park triggers
// another fitered wrt the Park triggers


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TLorentzVector.h>
#include "../interface/helper.h"

using namespace std;

constexpr bool debug = false;

class MuonTriggerSelector : public edm::stream::EDProducer<> {

  public:

    explicit MuonTriggerSelector(const edm::ParameterSet &iConfig);

    ~MuonTriggerSelector() override {};


  private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<std::vector<reco::Track>> displacedStandaloneMuonSrc_;
    const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

    // added
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

    // trigger muon matching
    const double max_deltaR_trigger_matching_;
    const double max_deltaPtRel_trigger_matching_;
    
    // for the sel muon
    const double selmu_ptMin_;          // min pT in all muons for B candidates
    const double selmu_absEtaMax_;      // max eta ""
    std::vector<std::string> HLTPaths_;
    //std::vector<std::string> L1Seeds_;
};


MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet &iConfig):
  bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
  muonSrc_(consumes<std::vector<pat::Muon>> (iConfig.getParameter<edm::InputTag>("muonCollection"))),
  displacedStandaloneMuonSrc_(consumes<std::vector<reco::Track>> (iConfig.getParameter<edm::InputTag>("displacedStandaloneMuonCollection"))),
  vertexSrc_( consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  // added
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  //
  max_deltaR_trigger_matching_(iConfig.getParameter<double>("max_deltaR_trigger_matching")),
  max_deltaPtRel_trigger_matching_(iConfig.getParameter<double>("max_deltaPtRel_trigger_matching")),
  selmu_ptMin_(iConfig.getParameter<double>("selmu_ptMin")),
  selmu_absEtaMax_(iConfig.getParameter<double>("selmu_absEtaMax")),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths"))
  //L1Seeds_(iConfig.getParameter<std::vector<std::string>>("L1seeds"))
{
  // produce 2 collections: trgMuons (tags) and SelectedMuons (probes & tags if survive preselection cuts)
  // trigger muons are slimmedMuons only
  // selected muons are slimmedMuons and displacedStandaloneMuons
  // make sure that the indices between the selected muons and transient tracks are consistent
  produces<pat::MuonCollection>("trgMuons"); 
  produces<pat::MuonCollection>("SelectedMuons");
  produces<TransientTrackCollection>("SelectedTransientMuons");  
}


void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    const auto &bField = iSetup.getData(bFieldToken_);

    std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
    std::unique_ptr<pat::MuonCollection> selmuons_out  ( new pat::MuonCollection );
    std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );
    
    edm::Handle<std::vector<pat::Muon>> slimmed_muons;
    iEvent.getByToken(muonSrc_, slimmed_muons);

    edm::Handle<std::vector<reco::Track>> displaced_standalone_muons;
    iEvent.getByToken(displacedStandaloneMuonSrc_, displaced_standalone_muons);

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(vertexSrc_, vertexHandle);
    const reco::VertexCollection& vertices = *vertexHandle;

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    if ( ! beamSpotHandle.isValid() ) {
      edm::LogError("MuonTriggerSelector") << "No beam spot available from Event" ;
    }  
    const reco::BeamSpot& beamSpot = *beamSpotHandle;

    std::vector<int> muonIsTrigger(slimmed_muons->size(), 0);
    std::vector<float> muonDR(slimmed_muons->size(),10000.);
    std::vector<float> muonDPT(slimmed_muons->size(),10000.);
    std::vector<int> loose_id(slimmed_muons->size(),0);

    std::vector<int> matched_reco_flag(slimmed_muons->size(),-1);
    std::vector<int> matched_trg_index(slimmed_muons->size(),-1);
    std::vector<float> matched_dr(slimmed_muons->size(),10000.);
    std::vector<float> matched_dpt(slimmed_muons->size(),-10000.);
    std::vector<std::vector<int>> fires;
    std::vector<std::vector<int>> prescales;
    std::vector<std::vector<float>> matcher; 
    std::vector<std::vector<float>> DR;
    std::vector<std::vector<float>> DPT;    

    // fetching the primary vertex
    const reco::Vertex & PV = vertexHandle->front();

    //std::cout << std::endl;
    //std::cout << std::endl << "number of muons:" << slimmed_muons->size() << std::endl;

    //////////////////////////////
    /* method 1 */
    // proceeding to the trigger muon matching
    // for that, only use slimmedMuons
    for(const pat::Muon &muon : *slimmed_muons){
        if(debug) std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi() << endl;
        //std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi() << " nTriggerObjectMatch " << muon.triggerObjectMatches().size() << endl;

        std::vector<int> frs(HLTPaths_.size(),0); //path fires for each reco muon
        std::vector<int> prescale(HLTPaths_.size(),-1);
        std::vector<float> temp_matched_to(HLTPaths_.size(),1000.);
        std::vector<float> temp_DR(HLTPaths_.size(),1000.);
        std::vector<float> temp_DPT(HLTPaths_.size(),1000.);
        int ipath=-1;

        // for each muon, we study whether it fires a HLT line or not
        for (const std::string& path: HLTPaths_){
            if(debug) std::cout << path << " " << muon.triggerObjectMatches().size() << std::endl;
            ipath++;

            // the following vectors are used in order to find the minimum DR between a reco muon and all the HLT objects that is matched to it
            // as a reco muon will be matched with only one HLT object every time, so there is a one-to-one correspondence between the two collections
            // DPt_rel is not used to create this one-to-one correspondence but only to create a few plots, debugging and be sure that everything is working fine. 
            std::vector<float> temp_dr(muon.triggerObjectMatches().size(),1000.);
            std::vector<float> temp_dpt(muon.triggerObjectMatches().size(),1000.);
            std::vector<float> temp_pt(muon.triggerObjectMatches().size(),1000.);

            // put the HLT path name into a wildcard
            char cstr[ (path+"*").size() + 1];
            strcpy( cstr, (path+"*").c_str() );       

            // Here we find all the HLT objects from each HLT path each time that are matched with the reco muon.
            if(muon.triggerObjectMatches().size()!=0){
                for(size_t i=0; i<muon.triggerObjectMatches().size();i++){
                  if(debug) std::cout << "triggerObjMatch " << i << " " << cstr << " " <<  muon.triggerObjectMatch(i)->hasPathName(cstr,true,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,true,false) << std::endl;
                    //std::cout << "trgObjMatch " << i << " " << cstr << " " <<  muon.triggerObjectMatch(i)->hasPathName(cstr,true,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,false,true) << " " << muon.triggerObjectMatch(i)->hasPathName(cstr,true,false) << std::endl;
                    // first bool is pathLastFilterAccepted, second is pathL3FilterAccepted
                    if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasPathName(cstr,true,true)){
                    //if(muon.triggerObjectMatch(i)!=0 && muon.triggerObjectMatch(i)->hasPathName(cstr,false,false) && DR[iMuo][path]<max_deltaR_ && fabs(DPT[iMuo][path])<max_deltaPtRel_ && DR[iMuo][path]!=10000){
                        //if(abs(muon.triggerObjectMatch(i)->eta())>1.5) std::cout << "eta=" <<muon.triggerObjectMatch(i)->eta();
                        frs[ipath]=1;
                        float dr = reco::deltaR(muon.triggerObjectMatch(i)->p4(), muon.p4()); 
                        float dpt=(muon.triggerObjectMatch(i)->pt()-muon.pt())/muon.triggerObjectMatch(i)->pt();
                        temp_dr[i]=dr;
                        temp_dpt[i]=dpt;
                        temp_pt[i]=muon.triggerObjectMatch(i)->pt();                   
                        if(debug)std::cout <<"Path=" <<cstr << endl;
                        if(debug)std::cout <<"HLT  Pt="<<muon.triggerObjectMatch(i)->pt() <<" Eta="<<muon.triggerObjectMatch(i)->eta() <<" Phi="<<muon.triggerObjectMatch(i)->phi() << endl;
                        //std::cout <<"HLT " << cstr << " Pt="<<muon.triggerObjectMatch(i)->pt() <<" Eta="<<muon.triggerObjectMatch(i)->eta() <<" Phi="<<muon.triggerObjectMatch(i)->phi() << endl;
                        if(debug)std::cout <<"Muon Pt="<< muon.pt() << " Eta=" << muon.eta() << " Phi=" << muon.phi()  <<endl;
                        if(debug)std::cout <<"DR = " << temp_dr[i] <<endl;
                        //std::cout <<"DR = " << temp_dr[i] << " DPT = " << temp_dpt[i] <<endl;
                    }
                }
                // and now we find the real minimum between the reco muon and all its matched HLT objects. 
                temp_DR[ipath]=*min_element(temp_dr.begin(),temp_dr.end());
                int position=std::min_element(temp_dr.begin(),temp_dr.end()) - temp_dr.begin();
                temp_DPT[ipath]=temp_dpt[position];
                temp_matched_to[ipath]=temp_pt[position];
                }
            }
        //and now since we have found the minimum DR we save a few variables for plots       
        fires.push_back(frs); //This is used in order to see if a reco muon fired a Trigger (1) or not (0).
        //for(unsigned int i(0); i<frs.size(); ++i){
        //  std::cout << "fires " << frs[i] << " deltaR " << temp_DR[i] << " deltaPt " << temp_DPT[i] << std::endl;
        //}
        prescales.push_back(prescale); //This is used for the initialisation of the prescales.
        matcher.push_back(temp_matched_to); //This is used in order to see if a reco muon is matched with a HLT object. PT of the reco muon is saved in this vector. 
        DR.push_back(temp_DR);
        DPT.push_back(temp_DPT);

    }

    //now, check for different reco muons that are matched to the same HLTObject.
    for(unsigned int path=0; path<HLTPaths_.size(); path++){
        for(unsigned int iMuo=0; iMuo<slimmed_muons->size(); iMuo++){
            for(unsigned int im=(iMuo+1); im<slimmed_muons->size(); im++){
                if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path]){
                //if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path] && DR[iMuo][path]<max_deltaR_ && fabs(DPT[iMuo][path])<max_deltaPtRel_ && DR[im][path]<max_deltaR_ && fabs(DPT[im][path])<max_deltaPtRel_){
                    if(DR[iMuo][path]<DR[im][path]){ //Keep the one that has the minimum DR with the HLT object
                        fires[im][path]=0;
                        matcher[im][path]=1000.;
                        DR[im][path]=1000.;                       
                        DPT[im][path]=1000.;
                    }
                    else{
                        fires[iMuo][path]=0;
                        matcher[iMuo][path]=1000.;
                        DR[iMuo][path]=1000.;                       
                        DPT[iMuo][path]=1000.;
                    }
                }              
            }
            if(matcher[iMuo][path]!=1000. && DR[iMuo][path]<max_deltaR_trigger_matching_ && fabs(DPT[iMuo][path])<max_deltaPtRel_trigger_matching_ && DR[iMuo][path]!=10000){
              muonIsTrigger[iMuo]=1;
              muonDR[iMuo]=DR[iMuo][path];
              muonDPT[iMuo]=DPT[iMuo][path];                
              
              for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
                //if(names.triggerName(i).find(HLTPaths_[path]) != std::string::npos || names.triggerName(i)==HLTPaths_[path]){
                if(names.triggerName(i).find(HLTPaths_[path]) != std::string::npos){
                  prescales[iMuo][path] = triggerPrescales->getPrescaleForIndex<double>(i);
                }
              }
            }
            else{
              fires[iMuo][path]=0;
              //prescales[iMuo][path] = -1;
            }
        }
    }

    if(debug)std::cout << "number of Muons=" <<slimmed_muons->size() << endl;
    // And now create a collection with trg muons from bParking line
    for(const pat::Muon & muon : *slimmed_muons){
        unsigned int iMuo(&muon -&(slimmed_muons->at(0)));

        if(muon.pt()<selmu_ptMin_) continue;
        if(fabs(muon.eta())>selmu_absEtaMax_) continue;

        if(muonIsTrigger[iMuo]==1){ 
            pat::Muon recoTriggerMuonCand(muon);
            trgmuons_out->emplace_back(recoTriggerMuonCand);
        }
    }

    //bool at_least_one_triggering_muon = 0;

    // add the slimmed muons to the collection 
    for(const pat::Muon & slimmed_muon : *slimmed_muons){
      unsigned int iMuo(&slimmed_muon - &(slimmed_muons->at(0)));
      if(slimmed_muon.pt()<selmu_ptMin_) continue;
      if(fabs(slimmed_muon.eta())>selmu_absEtaMax_) continue;
      //std::cout << "after the kin selection" << std::endl;

      const reco::TransientTrack muonTT((*(slimmed_muon.bestTrack())), &bField);
      if(!muonTT.isValid()) continue; // GM: and why do we skip this muon if muonTT is invalid? This seems to have no effect so I kept it.
      //std::cout << "after the muonTT selection" << std::endl;

      // transient track collection
      trans_muons_out->emplace_back(muonTT);

      //std::cout << "slimmed muon " << iMuo << " pt " << slimmed_muon.pt()  << " eta " << slimmed_muon.eta() << std::endl;
      
      // muon collection
      selmuons_out->emplace_back(slimmed_muon);

      for(unsigned int i=0; i<HLTPaths_.size(); i++){
        selmuons_out->back().addUserInt(HLTPaths_[i], fires[iMuo][i]);
        selmuons_out->back().addUserInt(HLTPaths_[i] + "_prescale", prescales[iMuo][i]);
      }
      selmuons_out->back().addUserInt("isTriggering", muonIsTrigger[iMuo]);
      selmuons_out->back().addUserFloat("DR", muonDR[iMuo]);
      selmuons_out->back().addUserFloat("DPT" ,muonDPT[iMuo]);

      //TODO check that the variables are called at the correct point
      //TODO add further variables
      selmuons_out->back().addUserFloat("dz", slimmed_muon.dB(slimmed_muon.PVDZ));
      selmuons_out->back().addUserFloat("dzS", slimmed_muon.dB(slimmed_muon.PVDZ)/slimmed_muon.edB(slimmed_muon.PVDZ));
      selmuons_out->back().addUserFloat("dxy", slimmed_muon.dB(slimmed_muon.PV2D));
      selmuons_out->back().addUserFloat("dxyS", slimmed_muon.dB(slimmed_muon.PV2D)/slimmed_muon.edB(slimmed_muon.PV2D));
      selmuons_out->back().addUserInt("softID", slimmed_muon.reco::Muon::passed(reco::Muon::SoftCutBasedId));

      // computing the IP with respect to the beamspot
      // first, a la R(D*)
      float dz_alaRdst = -99.;
      float dxy_BS_alaRdst = -99.;
      float dxyS_BS_alaRdst = -99.;
      if(!slimmed_muon.innerTrack().isNull()){
        auto tk = slimmed_muon.innerTrack();
        dz_alaRdst = tk->dz(PV.position());
        dxy_BS_alaRdst = tk->dxy(beamSpot);
        dxyS_BS_alaRdst = tk->dxy(beamSpot) / dxyError(*tk, beamSpot.position(), beamSpot.rotatedCovariance3D());
      }

      selmuons_out->back().addUserFloat("dz_alaRdst", dz_alaRdst);
      selmuons_out->back().addUserFloat("dxy_BS_alaRdst", dxy_BS_alaRdst);
      selmuons_out->back().addUserFloat("dxyS_BS_alaRdst", dxyS_BS_alaRdst);

      // then, by accounting for that the axis of the beamspot can be shifted
      // get the vertex the closest in dz to the muon
      const reco::Vertex& the_PV = PV; // initialise it to the first primary vertex
      float dist = -99;
      for(const reco::Vertex& vertex: vertices){
        // compute muon dz at this vertex
        float slimmed_muon_dz = abs(slimmed_muon.bestTrack()->dz(vertex.position()));
        if(dist == -99 || slimmed_muon_dz < dist){
          dist = slimmed_muon_dz;
          auto the_PV = vertex;
        }
      }
      auto bs_point = reco::Vertex::Point(beamSpot.x(the_PV.z()), beamSpot.y(the_PV.z()), beamSpot.z0());
      auto bs_error = beamSpot.covariance3D();
      float chi2 = 0.;
      float ndof = 0.;
      auto the_BS = reco::Vertex(bs_point, bs_error, chi2, ndof, 3);
        
      float mu_dxy_BS = -99.;
      float mu_dxyS_BS = -99.;
      if(slimmed_muon.bestTrack()){
        mu_dxy_BS = slimmed_muon.bestTrack()->dxy(the_BS.position());
        mu_dxyS_BS = slimmed_muon.bestTrack()->dxy(the_BS.position()) / dxyError(*slimmed_muon.bestTrack(), the_BS.position(), the_BS.error());
      }

      selmuons_out->back().addUserFloat("dxy_BS", mu_dxy_BS);
      selmuons_out->back().addUserFloat("dxyS_BS", mu_dxyS_BS);

      selmuons_out->back().addUserFloat("ip3d", fabs(slimmed_muon.dB(slimmed_muon.PV3D)));
      selmuons_out->back().addUserFloat("sip3d", fabs(slimmed_muon.dB(slimmed_muon.PV3D)/slimmed_muon.edB(slimmed_muon.PV3D)));

      selmuons_out->back().addUserFloat("segmentCompatibility", slimmed_muon.segmentCompatibility());
      selmuons_out->back().addUserFloat("validHitFraction", slimmed_muon.isGlobalMuon() || slimmed_muon.isTrackerMuon() ? slimmed_muon.innerTrack()->validFraction(): -1.);
      selmuons_out->back().addUserFloat("kinkFinderChi2", slimmed_muon.combinedQuality().trkKink);
      selmuons_out->back().addUserFloat("globalNormalisedChi2", slimmed_muon.isGlobalMuon() ? slimmed_muon.globalTrack()->normalizedChi2(): -1.);
      selmuons_out->back().addUserFloat("localPositionChi2", slimmed_muon.combinedQuality().chi2LocalPosition);
      selmuons_out->back().addUserFloat("caloCompatibility", slimmed_muon.caloCompatibility());
      selmuons_out->back().addUserInt("numberOfValidMuonHits", slimmed_muon.isGlobalMuon() ? slimmed_muon.globalTrack()->hitPattern().numberOfValidMuonHits(): -1);
      selmuons_out->back().addUserInt("numberOfValidPixelHits", slimmed_muon.isGlobalMuon() || slimmed_muon.isTrackerMuon() ? slimmed_muon.innerTrack()->hitPattern().numberOfValidPixelHits(): -1);
      selmuons_out->back().addUserInt("numberOfTrackerLayers", slimmed_muon.isGlobalMuon() || slimmed_muon.isTrackerMuon() ? slimmed_muon.innerTrack()->hitPattern().trackerLayersWithMeasurement(): -1);
      selmuons_out->back().addUserInt("numberOfPixelLayers", slimmed_muon.isGlobalMuon() || slimmed_muon.isTrackerMuon() ? slimmed_muon.innerTrack()->hitPattern().pixelLayersWithMeasurement(): -1);
      selmuons_out->back().addUserInt("trackerHighPurityFlag", slimmed_muon.isGlobalMuon() || slimmed_muon.isTrackerMuon() ? slimmed_muon.innerTrack()->quality(reco::TrackBase::highPurity): -1);
      selmuons_out->back().addUserInt("nStations", slimmed_muon.numberOfMatchedStations());
    }

    //std::cout << "at least one triggering muon: " << at_least_one_triggering_muon << std::endl;

    //std::cout << std::endl;

    iEvent.put(std::move(trgmuons_out), "trgMuons"); 
    iEvent.put(std::move(selmuons_out), "SelectedMuons");
    iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
}



DEFINE_FWK_MODULE(MuonTriggerSelector);
