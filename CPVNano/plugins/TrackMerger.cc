// Merges the PFPackedCandidates and Lost tracks
// beam spot readout in case dcasig to be calculated wrt beam spot
// currently computed wrt triggeringMuon vertex

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include "../interface/helper.h"

class TrackMerger : public edm::global::EDProducer<> {


public:

  explicit TrackMerger(const edm::ParameterSet &cfg):
    bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    beamSpotSrc_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
    tracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    lostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    trgMuonToken_(consumes<std::vector<pat::Muon>>(cfg.getParameter<edm::InputTag>("trgMuon"))),
    muonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
    eleToken_(consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("pfElectrons"))),
    vertexToken_(consumes<reco::VertexCollection> (cfg.getParameter<edm::InputTag>( "vertices" ))), 
    trkPtCut_(cfg.getParameter<double>("trkPtCut")),
    trkEtaCut_(cfg.getParameter<double>("trkEtaCut")),
    dzTrg_cleaning_(cfg.getParameter<double>("dzTrg_cleaning")),
    do_trgmu_cleaning_(cfg.getParameter<bool>("do_trgmu_cleaning")),
    do_mu_cleaning_(cfg.getParameter<bool>("do_mu_cleaning")),
    do_el_cleaning_(cfg.getParameter<bool>("do_el_cleaning")),
    do_trk_highpurity_(cfg.getParameter<bool>("do_trk_highpurity")),
    dcaSig_(cfg.getParameter<double>("dcaSig")),
    trkNormChiMin_(cfg.getParameter<int>("trkNormChiMin")),
    trkNormChiMax_(cfg.getParameter<int>("trkNormChiMax")) 
{
    produces<pat::CompositeCandidateCollection>("SelectedTracks");  
    produces<TransientTrackCollection>("SelectedTransientTracks");  
}

  ~TrackMerger() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  

private:
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
  const edm::EDGetTokenT<std::vector<pat::Muon>> trgMuonToken_;
  const edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> eleToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  //selections                                                                 
  const double trkPtCut_;
  const double trkEtaCut_;
  const double dzTrg_cleaning_;
  const bool do_trgmu_cleaning_;
  const bool do_mu_cleaning_;
  const bool do_el_cleaning_;
  const bool do_trk_highpurity_;
  const double dcaSig_;
  const int trkNormChiMin_;
  const int trkNormChiMax_;
};






void TrackMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &stp) const {

  //input
  const auto &bField = stp.getData(bFieldToken_);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  evt.getByToken(beamSpotSrc_, beamSpotHandle);
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKstllProducer") << "No beam spot available from Event" ;
  }  
  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  edm::Handle<pat::PackedCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);
  edm::Handle<pat::PackedCandidateCollection> lostTracks;
  evt.getByToken(lostTracksToken_, lostTracks);
  edm::Handle<std::vector<pat::Muon>> trgMuons;
  evt.getByToken(trgMuonToken_, trgMuons);

  edm::Handle<std::vector<pat::Muon>> muons;
  evt.getByToken(muonToken_, muons);
  edm::Handle<pat::ElectronCollection> pfele;
  evt.getByToken(eleToken_, pfele);
  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexToken_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();
  const reco::VertexCollection& vertices = *vertexHandle;

  //for lost tracks / pf discrimination
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();

  //ok this was CompositeCandidateCollection 
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_out      (new pat::CompositeCandidateCollection);
  std::unique_ptr<TransientTrackCollection>          trans_tracks_out(new TransientTrackCollection);

   std::vector< std::pair<pat::CompositeCandidate,reco::TransientTrack> > vectrk_ttrk; 
 
  // for loop is better to be range based - especially for large ensembles  
  for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
    const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk-nTracks];

    //arranging cuts for speed
    if(!trk.hasTrackDetails()) continue;
    if(abs(trk.pdgId()) != 211) continue; //do we want also to keep muons?
    if(trk.pt() < trkPtCut_ ) continue;
    if(fabs(trk.eta()) > trkEtaCut_) continue;

    if( (trk.bestTrack()->normalizedChi2() < trkNormChiMin_ &&
         trkNormChiMin_>=0 ) ||
        (trk.bestTrack()->normalizedChi2() > trkNormChiMax_ &&
         trkNormChiMax_>0)  )    continue; 

    // keep only tracks close to trigger muon
    bool skipTrack = true;
    float dzTrg = -99.;
    float drTrg = -99.;
    if(do_trgmu_cleaning_){
      for (const pat::Muon & mu: *trgMuons){
        drTrg = reco::deltaR(trk, mu);
        dzTrg = trk.vz() - mu.vz();

        if (fabs(dzTrg) < dzTrg_cleaning_){
          skipTrack = false;
          break;
        }
      }
      if (skipTrack) continue;
    }


    // high purity requirement applied only in packedCands
    if( do_trk_highpurity_ && iTrk < nTracks && !trk.trackHighPurity()) continue;
    const reco::TransientTrack trackTT((*trk.bestTrack()), &bField);

    // distance closest approach in x,y wrt beam spot
    std::pair<double,double> DCA = computeDCA(trackTT, beamSpot);
    float DCABS = DCA.first;
    float DCABSErr = DCA.second;
    float DCASig = (DCABSErr != 0 && float(DCABSErr) == DCABSErr) ? fabs(DCABS/DCABSErr) : -1;
    if (DCASig >  dcaSig_  && dcaSig_ >0) continue;
    
    // Compute corrected DCA
    // By accounting for that the axis of the beamspot can be shifted get the vertex the closest in dz to the muon
    const reco::Vertex& the_PV = PV; // initialise it to the first primary vertex
    float dist = -99;
    for(const reco::Vertex& vertex: vertices){
      // compute muon dz at this vertex
      float track_dz = abs(trk.dz(vertex.position()));
      if(dist == -99 || track_dz < dist){
        dist = track_dz;
        auto the_PV = vertex;
      }
    }

    std::pair<double,double> DCA_corr = computeDCA(trackTT, beamSpot, the_PV);
    float DCABS_corr = DCA_corr.first;
    float DCABSErr_corr = DCA_corr.second;
    float DCASig_corr = (DCABSErr_corr != 0 && float(DCABSErr_corr) == DCABSErr_corr) ? fabs(DCABS_corr/DCABSErr_corr) : -1;

    // clean tracks wrt to all muons
    int matchedToMuon       = 0;
    if(do_mu_cleaning_){
      for (const pat::Muon &imutmp : *muons) {
          for (unsigned int i = 0; i < imutmp.numberOfSourceCandidatePtrs(); ++i) {
              if (! ((imutmp.sourceCandidatePtr(i)).isNonnull() && 
                     (imutmp.sourceCandidatePtr(i)).isAvailable())
                 )   continue;
              
              const edm::Ptr<reco::Candidate> & source = imutmp.sourceCandidatePtr(i);
              if (source.id() == tracks.id() && source.key() == iTrk){
                  matchedToMuon =1;
                  break;
              }
          }
      }
    }

    // clean tracks wrt to all pf electrons
    int matchedToEle = 0;
    if(do_el_cleaning_){
      for (const pat::Electron &ietmp : *pfele) {
          for (unsigned int i = 0; i < ietmp.numberOfSourceCandidatePtrs(); ++i) {
              if (! ((ietmp.sourceCandidatePtr(i)).isNonnull() && 
                     (ietmp.sourceCandidatePtr(i)).isAvailable())
                 )   continue;
              const edm::Ptr<reco::Candidate> & source = ietmp.sourceCandidatePtr(i);
              if (source.id() == tracks.id() && source.key() == iTrk){
                  matchedToEle =1;
                  break;
              }        
          }
       }
    }

    // clean tracks wrt to all low pT electrons
    //edm::Ptr<pat::PackedCandidate> track;
    //if ( iTrk < nTracks ) { track = edm::Ptr<pat::PackedCandidate>(tracks,iTrk); }
    //else { track = edm::Ptr<pat::PackedCandidate>(lostTracks,iTrk-nTracks); }
    //int matchedToLowPtEle = 0;
    //for ( auto const& ele : *lowptele ) {
    //  reco::GsfTrackRef gsf = ele.gsfTrack();
    //  if ( iTrk < nTracks ) {
	//edm::Ptr<pat::PackedCandidate> packed = edm::refToPtr( (*gsf2packed)[gsf] );
	//if ( track == packed ) { matchedToLowPtEle = 1; break; }
  //    } else {
	//edm::Ptr<pat::PackedCandidate> lost = edm::refToPtr( (*gsf2lost)[gsf] );
	//if ( track == lost ) { matchedToLowPtEle = 1; break; }
  //    }
  //  }

    pat::CompositeCandidate pcand;
    pcand.setP4(trk.p4());
    pcand.setCharge(trk.charge());
    pcand.setVertex(trk.vertex());
    pcand.setPdgId(trk.pdgId());
    pcand.addUserInt("isPacked", (iTrk < nTracks));
    pcand.addUserInt("isLostTrk", (iTrk < nTracks) ? 0 : 1);      
    // IP computed wrt to PV ref, as in https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/PackedCandidate.h?v=CMSSW_12_5_X_2022-08-01-2300&%21v=CMSSW_10_6_20
    // checked that it is compatible (avg e-3 accuracy) to the IP evaluated at the first PV

    pcand.addUserFloat("dxy", trk.bestTrack()->dxy(the_PV.position()));
    pcand.addUserFloat("dxyS", trk.bestTrack()->dxy(the_PV.position())/trk.bestTrack()->dxyError(the_PV.position(), the_PV.covariance()));
    pcand.addUserFloat("dz", trk.bestTrack()->dz(the_PV.position())); 
    pcand.addUserFloat("dzS", trk.bestTrack()->dz(the_PV.position())/trk.bestTrack()->dzError());
    pcand.addUserFloat("DCASig", DCASig);
    pcand.addUserFloat("DCASig_corr", DCASig_corr);

    pcand.addUserFloat("ptErr", trk.bestTrack()->ptError());
    pcand.addUserFloat("drTrg", drTrg);
    pcand.addUserFloat("dzTrg", dzTrg);
    pcand.addUserFloat("chi2", trk.bestTrack()->chi2()); 
    pcand.addUserInt("ndof", trk.bestTrack()->ndof()); 
    pcand.addUserFloat("normalisedChi2", trk.bestTrack()->normalizedChi2()); 
    pcand.addUserInt("numberOfValidHits", trk.bestTrack()->numberOfValidHits()); 
    pcand.addUserInt("numberOfLostHits", trk.bestTrack()->numberOfLostHits()); 
    pcand.addUserInt("numberOfPixelHits", trk.numberOfPixelHits()); 
    pcand.addUserInt("numberOfValidPixelHits", trk.bestTrack()->hitPattern().numberOfValidPixelHits()); 
    pcand.addUserInt("numberOfTrackerLayers", trk.bestTrack()->hitPattern().trackerLayersWithMeasurement()); 
    pcand.addUserInt("numberOfPixelLayers", trk.bestTrack()->hitPattern().pixelLayersWithMeasurement()); 
    pcand.addUserInt("qualityIndex", trk.bestTrack()->qualityMask()); 
    pcand.addUserInt("highPurityFlag", trk.bestTrack()->quality(reco::TrackBase::highPurity)); 
    pcand.addUserFloat("validFraction", trk.bestTrack()->validFraction()); 
    pcand.addUserInt("isMatchedToMuon", matchedToMuon);
    pcand.addUserInt("isMatchedToEle", matchedToEle);
    //pcand.addUserInt("isMatchedToLowPtEle", matchedToLowPtEle);
    pcand.addUserInt("nValidHits", trk.bestTrack()->found());
    pcand.addUserInt("keyPacked", iTrk);

    // Covariance matrix elements for helix parameters for decay time uncertainty
    pcand.addUserFloat("covQopQop", trk.bestTrack()->covariance(0, 0));
    pcand.addUserFloat("covLamLam", trk.bestTrack()->covariance(1, 1));
    pcand.addUserFloat("covPhiPhi", trk.bestTrack()->covariance(2, 2));
    pcand.addUserFloat("covQopLam", trk.bestTrack()->covariance(0, 1));
    pcand.addUserFloat("covQopPhi", trk.bestTrack()->covariance(0, 2));
    pcand.addUserFloat("covLamPhi", trk.bestTrack()->covariance(1, 2));

    
    if ( iTrk < nTracks )
      pcand.addUserCand( "cand", edm::Ptr<pat::PackedCandidate> ( tracks, iTrk ));
    else 
      pcand.addUserCand( "cand", edm::Ptr<pat::PackedCandidate> ( lostTracks, iTrk-nTracks ));   
 
  //in order to avoid revoking the sxpensive ttrack builder many times and still have everything sorted, we add them to vector of pairs
   vectrk_ttrk.emplace_back( std::make_pair(pcand,trackTT ) );   
  }

  // sort to be uniform with leptons
  std::sort( vectrk_ttrk.begin(), vectrk_ttrk.end(), 
             [] ( auto & trk1, auto & trk2) -> 
                  bool {return (trk1.first).pt() > (trk2.first).pt();} 
           );

  // finnaly save ttrks and trks to the correct _out vectors
  for ( auto & trk: vectrk_ttrk){
    tracks_out -> emplace_back( trk.first);
    trans_tracks_out -> emplace_back(trk.second);
  }

  evt.put(std::move(tracks_out),       "SelectedTracks");
  evt.put(std::move(trans_tracks_out), "SelectedTransientTracks");
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackMerger);
