// -*- C++ -*-
//
// Package:    L1PixelTrigger
// Class:      L1PixelTrigger
// 
/**\class L1PixelTrigger L1PixelTrigger.h 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jose Cupertino Ruiz Vargas,32 4-C14,+41227674949,
//         Created:  Wed Aug  7 11:57:33 CEST 2013
// $Id$
//
//

#include "L1PixelTrigger.h"

L1PixelTrigger::L1PixelTrigger(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  eGamma_tag = iConfig.getParameter<edm::InputTag>("eGamma");
  eGammaCrystal_tag = iConfig.getParameter<edm::InputTag>("eGammaCrystal");
  edm::Service<TFileService> fs;
  t = fs->make<TTree>("t","t");
  t->Branch("run",               &run);
  t->Branch("event",             &event);
  t->Branch("bunchN",            &bunch_n);
  t->Branch("pileup",            &pileup);
  t->Branch("beamSpotX0",        &beamspot_x0);
  t->Branch("beamSpotY0",        &beamspot_y0);
  t->Branch("beamSpotZ0",        &beamspot_z0);
  t->Branch("beamSpotX0Error",   &beamspot_x0Error);
  t->Branch("beamSpotY0Error",   &beamspot_y0Error);
  t->Branch("beamSpotZ0Error",   &beamspot_z0Error);
  t->Branch("beamWidthX",        &beam_widthX);
  t->Branch("beamWidthY",        &beam_widthY);
  t->Branch("beamSigmaZ",        &beam_sigmaZ);
  t->Branch("beamWidthXError",   &beam_widthXError);
  t->Branch("beamWidthYError",   &beam_widthYError);
  t->Branch("beamSigmaZError",   &beam_sigmaZError);      
  t->Branch("genPartN",          &genpart_n);
  t->Branch("egCrysN",           &egammaC_n);
  t->Branch("egN",               &egamma_n);
  t->Branch("genPartEt",         &genpart_et);
  t->Branch("genPartPt",         &genpart_pt);
  t->Branch("genPartEta",        &genpart_eta);
  t->Branch("genPartPhi",        &genpart_phi);
  t->Branch("genPartCharge",     &genpart_charge);
  t->Branch("genPartId",         &genpart_id);
  t->Branch("egCrysE",           &egammaC_e);
  t->Branch("egCrysEt",          &egammaC_et);
  t->Branch("egCrysEta",         &egammaC_eta);
  t->Branch("egCrysPhi",         &egammaC_phi);
  t->Branch("egCrysGx",          &egammaC_gx);
  t->Branch("egCrysGy",          &egammaC_gy);
  t->Branch("egCrysGz",          &egammaC_gz);
  t->Branch("egCrysCharge",      &egammaC_charge);
  t->Branch("egE",               &egamma_e);
  t->Branch("egEt",              &egamma_et);
  t->Branch("egEta",             &egamma_eta);
  t->Branch("egPhi",             &egamma_phi);
  t->Branch("egGx",              &egamma_gx);
  t->Branch("egGy",              &egamma_gy);
  t->Branch("egGz",              &egamma_gz);
  t->Branch("egCharge",          &egamma_charge);
  t->Branch("fHitN",             &fHit_n);
  t->Branch("bHitN",             &bHit_n);
  t->Branch("fHitSubid",         &fHit_subid);
  t->Branch("fHitDisk",          &fHit_disk);
  t->Branch("fHitBlade",         &fHit_blade);
  t->Branch("fHitModule",        &fHit_module);
  t->Branch("fHitSide",          &fHit_side);
  t->Branch("bHitSubid",         &bHit_subid);
  t->Branch("bHitLayer",         &bHit_layer);
  t->Branch("bHitLadder",        &bHit_ladder);
  t->Branch("bHitModule",        &bHit_module);
  t->Branch("fClSize",           &fCl_size);
  t->Branch("bClSize",           &bCl_size);
  t->Branch("fHitGx",            &fHit_gx);
  t->Branch("fHitGy",            &fHit_gy);
  t->Branch("fHitGz",            &fHit_gz);
  t->Branch("bHitGx",            &bHit_gx);
  t->Branch("bHitGy",            &bHit_gy);
  t->Branch("bHitGz",            &bHit_gz);
}

L1PixelTrigger::~L1PixelTrigger()
{
}

// ------------ method called for each event  ------------
void L1PixelTrigger::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  InitializeVectors();
  run = e.id().run();
  event = e.id().event();
  ///////////////////////////////////////////////////////////
  // Pileup  
  //////////////////////////////////////////////////////////
  edm::Handle< std::vector<PileupSummaryInfo> > puinfo;
  e.getByLabel( "addPileupInfo", puinfo );
  bunch_n = puinfo->size();
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  for(PVI = puinfo->begin(); PVI != puinfo->end(); ++PVI) {
     pileup.push_back(PVI->getPU_NumInteractions());
  }
  ///////////////////////////////////////////////////////////
  // Beam Spot  
  //////////////////////////////////////////////////////////
  edm::Handle<reco::BeamSpot> thebeamSpot;
  e.getByLabel( "BeamSpotFromSim", "BeamSpot", thebeamSpot );
  beamspot_x0      = thebeamSpot->x0();
  beamspot_y0      = thebeamSpot->y0();
  beamspot_z0      = thebeamSpot->z0();
  beamspot_x0Error = thebeamSpot->x0Error();
  beamspot_y0Error = thebeamSpot->y0Error();
  beamspot_z0Error = thebeamSpot->z0Error();
  beam_widthX      = thebeamSpot->BeamWidthX();
  beam_widthY      = thebeamSpot->BeamWidthY();
  beam_sigmaZ      = thebeamSpot->sigmaZ();
  beam_widthXError = thebeamSpot->BeamWidthXError();
  beam_widthYError = thebeamSpot->BeamWidthYError();
  beam_sigmaZError = thebeamSpot->sigmaZ0Error();
  ///////////////////////////////////////////////////////////
  // Gen Particles  
  //////////////////////////////////////////////////////////
  edm::Handle< reco::GenParticleCollection > genParticles;
  e.getByLabel( "genParticles", genParticles );
  genpart_n = genParticles->size() ;
  for ( int i = 0; i < genpart_n; ++i ){
      const reco::GenParticle & genParticle = genParticles->at(i);
      genpart_et.push_back(     genParticle.et()     );
      genpart_pt.push_back(     genParticle.pt()     );
      genpart_eta.push_back(    genParticle.eta()    );
      genpart_phi.push_back(    genParticle.phi()    );
      genpart_charge.push_back( genParticle.charge() );
      genpart_id.push_back(     genParticle.pdgId()  );
  }
  ///////////////////////////////////////////////////////////
  // Ecal Single Crystal 
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > EgammaC;
  e.getByLabel( eGammaCrystal_tag, EgammaC );
  egammaC_n = EgammaC->size();
  for(l1extra::L1EmParticleCollection::const_iterator it = EgammaC->begin(); it!=EgammaC->end(); ++it){
    egammaC_e.push_back(it->energy());
    egammaC_et.push_back(it->et());
    egammaC_eta.push_back(it->eta());
    egammaC_phi.push_back(it->phi());
    egammaC_charge.push_back(it->charge());
    GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
    egammaC_gx.push_back(pos.x());
    egammaC_gy.push_back(pos.y());
    egammaC_gz.push_back(pos.z());
  }
  ///////////////////////////////////////////////////////////
  // L1 Ecal Trigger Primitives
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > Egamma;
  e.getByLabel( eGamma_tag, Egamma );
  egamma_n = Egamma->size();
  for(l1extra::L1EmParticleCollection::const_iterator it = Egamma->begin(); it!=Egamma->end(); ++it){
    egamma_e.push_back(it->energy());
    egamma_et.push_back(it->et());
    egamma_eta.push_back(it->eta());
    egamma_phi.push_back(it->phi());
    egamma_charge.push_back(it->charge());
    GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
    egamma_gx.push_back(pos.x());
    egamma_gy.push_back(pos.y());
    egamma_gz.push_back(pos.z());
  }
  //////////////////////////////////////////////////////////
  // Geometry 
  //////////////////////////////////////////////////////////
  edm::ESHandle<TrackerGeometry> geom;
  es.get<TrackerDigiGeometryRecord>().get(geom);
  //////////////////////////////////////////////////////////
  // RecHits 
  //////////////////////////////////////////////////////////
  edm::Handle<SiPixelRecHitCollection> recHits;
  e.getByLabel( "siPixelRecHits", recHits );
  bHit_n = 0;
  fHit_n = 0;
  SiPixelRecHitCollection::const_iterator detUnitIt    = recHits->begin();
  SiPixelRecHitCollection::const_iterator detUnitItEnd = recHits->end();
  for ( ; detUnitIt != detUnitItEnd; detUnitIt++ ) {
      DetId detId = DetId(detUnitIt->detId()); 
      int subid = detId.subdetId();
      SiPixelRecHitCollection::DetSet::const_iterator recHitIt    = detUnitIt->begin();
      SiPixelRecHitCollection::DetSet::const_iterator recHitItEnd = detUnitIt->end();
      for ( ; recHitIt != recHitItEnd; ++recHitIt) {
          LocalPoint  lp = recHitIt->localPosition();
          GlobalPoint gp = ( (geom.product())->idToDet(detId) )->surface().toGlobal(lp);
          SiPixelRecHit::ClusterRef const& Cluster = recHitIt->cluster();
          if ( gp.perp() < 20 ){ // reject outer tracker
             if ( subid==PixelSubdetector::PixelBarrel ){
                bHit_n ++;
                bHit_gx.push_back( gp.x() );
                bHit_gy.push_back( gp.y() );
                bHit_gz.push_back( gp.z() );
                bHit_subid.push_back(subid);
                bHit_layer.push_back(  PXBDetId(detId).layer() ); 
                bHit_ladder.push_back( PXBDetId(detId).ladder() ); 
                bHit_module.push_back( PXBDetId(detId).module() ); 
                bCl_size.push_back(Cluster->size());
             }
             if ( subid==PixelSubdetector::PixelEndcap ){
                fHit_n ++;
                fHit_gx.push_back( gp.x() );
                fHit_gy.push_back( gp.y() );
                fHit_gz.push_back( gp.z() );
                fHit_subid.push_back(subid);
                fHit_disk.push_back(   PXFDetId(detId).disk()   );
                fHit_blade.push_back(  PXFDetId(detId).blade()  );
                fHit_module.push_back( PXFDetId(detId).module() );
                fHit_side.push_back(   PXFDetId(detId).side()   );
                fCl_size.push_back(Cluster->size());
             }
          }
      } // close recHits loop
  } // close detUnits loop
  
  t->Fill();

}// close L1PixelTrigger::analyze

GlobalPoint L1PixelTrigger::getCalorimeterPosition(double phi, double eta, double e) {
  double x = 0;
  double y = 0;
  double z = 0;
  double depth = 0.89*(7.7+ log(e) );
  double theta = 2*atan(exp(-1*eta));
  double r = 0;
  if( fabs(eta) > 1.479 )
    {
      double ecalZ = 314*fabs(eta)/eta;

      r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  else
    {
      double rperp = 129.0;
      double zface =  sqrt( cos( theta ) * cos( theta ) /
                            ( 1 - cos( theta ) * cos( theta ) ) *
                            rperp * rperp ) * abs( eta ) / eta;
      r = sqrt( rperp * rperp + zface * zface ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  GlobalPoint pos(x,y,z);
  return pos;
}
void L1PixelTrigger::InitializeVectors()
{
             pileup.clear();
         genpart_et.clear();
         genpart_pt.clear();
        genpart_eta.clear();
        genpart_phi.clear();
     genpart_charge.clear();
         genpart_id.clear();
          egammaC_e.clear();
         egammaC_et.clear();
        egammaC_eta.clear();
        egammaC_phi.clear();
         egammaC_gx.clear();
         egammaC_gy.clear();
         egammaC_gz.clear();
     egammaC_charge.clear();
           egamma_e.clear();
          egamma_et.clear();
         egamma_eta.clear();
         egamma_phi.clear();
          egamma_gx.clear();
          egamma_gy.clear();
          egamma_gz.clear();
      egamma_charge.clear();
         bHit_subid.clear();
         bHit_layer.clear();
        bHit_ladder.clear();
        bHit_module.clear();
           bCl_size.clear();
            bHit_gx.clear();
            bHit_gy.clear();
            bHit_gz.clear();
         fHit_subid.clear();
          fHit_disk.clear();
         fHit_blade.clear();
          fHit_side.clear();
           fCl_size.clear();
            fHit_gx.clear();
            fHit_gy.clear();
            fHit_gz.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PixelTrigger);
