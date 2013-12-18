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

L1PixelTrigger::L1PixelTrigger(const edm::ParameterSet& iConfig) : conf_(iConfig) 

{
   //now do what ever initialization is needed
  eGamma_tag = iConfig.getParameter<edm::InputTag>("eGamma");
  eGammaCrystal_tag = iConfig.getParameter<edm::InputTag>("eGammaCrystal");
  edm::Service<TFileService> fs;
  t = fs->make<TTree>("t","t");
  t->Branch("run",                 &run);
  t->Branch("event",               &event);
  t->Branch("bunchN",              &bunch_n);
  t->Branch("bunchId",             &bunch_id);
  t->Branch("pileup",              &pileup);
  t->Branch("beamSpotX0",          &beamspot_x0);
  t->Branch("beamSpotY0",          &beamspot_y0);
  t->Branch("beamSpotZ0",          &beamspot_z0);
  t->Branch("beamSpotX0Error",     &beamspot_x0Error);
  t->Branch("beamSpotY0Error",     &beamspot_y0Error);
  t->Branch("beamSpotZ0Error",     &beamspot_z0Error);
  t->Branch("beamWidthX",          &beam_widthX);
  t->Branch("beamWidthY",          &beam_widthY);
  t->Branch("beamSigmaZ",          &beam_sigmaZ);
  t->Branch("beamWidthXError",     &beam_widthXError);
  t->Branch("beamWidthYError",     &beam_widthYError);
  t->Branch("beamSigmaZError",     &beam_sigmaZError);      
  t->Branch("genPartN",            &genpart_n);
  t->Branch("genPartEt",           &genpart_et);
  t->Branch("genPartPt",           &genpart_pt);
  t->Branch("genPartEta",          &genpart_eta);
  t->Branch("genPartPhi",          &genpart_phi);
  t->Branch("genPartCharge",       &genpart_charge);
  t->Branch("genPartId",           &genpart_id);
  t->Branch("egCrysN",             &egammaC_n);
  t->Branch("egCrysE",             &egammaC_e);
  t->Branch("egCrysEt",            &egammaC_et);
  t->Branch("egCrysEta",           &egammaC_eta);
  t->Branch("egCrysPhi",           &egammaC_phi);
  t->Branch("egCrysGx",            &egammaC_gx);
  t->Branch("egCrysGy",            &egammaC_gy);
  t->Branch("egCrysGz",            &egammaC_gz);
  t->Branch("egCrysCharge",        &egammaC_charge);
  t->Branch("egN",                 &egamma_n);
  t->Branch("egE",                 &egamma_e);
  t->Branch("egEt",                &egamma_et);
  t->Branch("egEta",               &egamma_eta);
  t->Branch("egPhi",               &egamma_phi);
  t->Branch("egGx",                &egamma_gx);
  t->Branch("egGy",                &egamma_gy);
  t->Branch("egGz",                &egamma_gz);
  t->Branch("egCharge",            &egamma_charge);
  t->Branch("bRecHitN",            &bRecHit_n);
  t->Branch("bRecHitSubid",        &bRecHit_subid);
  t->Branch("bRecHitLayer",        &bRecHit_layer);
  t->Branch("bRecHitLadder",       &bRecHit_ladder);
  t->Branch("bRecHitModule",       &bRecHit_module);
  t->Branch("bRecHitLx",           &bRecHit_lx);
  t->Branch("bRecHitLy",           &bRecHit_ly);
  t->Branch("bRecHitLxx",          &bRecHit_lxx);
  t->Branch("bRecHitLxy",          &bRecHit_lxy);
  t->Branch("bRecHitLyy",          &bRecHit_lyy);
  t->Branch("bRecHitGx",           &bRecHit_gx);
  t->Branch("bRecHitGy",           &bRecHit_gy);
  t->Branch("bRecHitGz",           &bRecHit_gz);
  t->Branch("bClSize",             &bCl_size);
  t->Branch("bClSizeX",            &bCl_sizeX);
  t->Branch("bClSizeY",            &bCl_sizeY);
  t->Branch("bClCharge",           &bCl_charge);
  t->Branch("bClX",                &bCl_x);
  t->Branch("bClY",                &bCl_y);
  t->Branch("bSimHitN",            &bSimHit_n);
  t->Branch("bSimHitSubid",        &bSimHit_subid);
  t->Branch("bSimHitLayer",        &bSimHit_layer);
  t->Branch("bSimHitLadder",       &bSimHit_ladder);
  t->Branch("bSimHitModule",       &bSimHit_module);
  t->Branch("bSimHitPiD",          &bSimHit_pid);
  t->Branch("bSimHitEloss",        &bSimHit_eloss);
  t->Branch("bSimHitToF",          &bSimHit_tof);
  t->Branch("bSimHitLx",           &bSimHit_lx);
  t->Branch("bSimHitLy",           &bSimHit_ly);
  t->Branch("bSimHitGx",           &bSimHit_gx);
  t->Branch("bSimHitGy",           &bSimHit_gy);
  t->Branch("bSimHitGz",           &bSimHit_gz);
  t->Branch("bClosestSimHitPiD",   &bClosestSimHit_pid);
  t->Branch("bClosestSimHitEloss", &bClosestSimHit_eloss);
  t->Branch("bClosestSimHitToF",   &bClosestSimHit_tof);
  t->Branch("bClosestSimHitLx",    &bClosestSimHit_lx);
  t->Branch("bClosestSimHitLy",    &bClosestSimHit_ly);
  t->Branch("bClosestSimHitGx",    &bClosestSimHit_gx);
  t->Branch("bClosestSimHitGy",    &bClosestSimHit_gy);
  t->Branch("bClosestSimHitGz",    &bClosestSimHit_gz);
  t->Branch("fRecHitN",            &fRecHit_n);
  t->Branch("fRecHitSubid",        &fRecHit_subid);
  t->Branch("fRecHitDisk",         &fRecHit_disk);
  t->Branch("fRecHitBlade",        &fRecHit_blade);
  t->Branch("fRecHitModule",       &fRecHit_module);
  t->Branch("fRecHitSide",         &fRecHit_side);
  t->Branch("fRecHitLx",           &fRecHit_lx);
  t->Branch("fRecHitLy",           &fRecHit_ly);
  t->Branch("fRecHitLxx",          &fRecHit_lxx);
  t->Branch("fRecHitLxy",          &fRecHit_lxy);
  t->Branch("fRecHitLyy",          &fRecHit_lyy);
  t->Branch("fRecHitGx",           &fRecHit_gx);
  t->Branch("fRecHitGy",           &fRecHit_gy);
  t->Branch("fRecHitGz",           &fRecHit_gz);
  t->Branch("fClSize",             &fCl_size);
  t->Branch("fClSizeX",            &fCl_sizeX);
  t->Branch("fClSizeY",            &fCl_sizeY);
  t->Branch("fClCharge",           &fCl_charge);
  t->Branch("fClX",                &fCl_x);
  t->Branch("fClY",                &fCl_y);
  t->Branch("fSimHitN",            &fSimHit_n);
  t->Branch("fSimHitSubid",        &fSimHit_subid);
  t->Branch("fSimHitDisk",         &fSimHit_disk);
  t->Branch("fSimHitBlade",        &fSimHit_blade);
  t->Branch("fSimHitModule",       &fSimHit_module);
  t->Branch("fSimHitSide",         &fSimHit_side);
  t->Branch("fSimHitPiD",          &fSimHit_pid);
  t->Branch("fSimHitEloss",        &fSimHit_eloss);
  t->Branch("fSimHitToF",          &fSimHit_tof);
  t->Branch("fSimHitLx",           &fSimHit_lx);
  t->Branch("fSimHitLy",           &fSimHit_ly);
  t->Branch("fSimHitGx",           &fSimHit_gx);
  t->Branch("fSimHitGy",           &fSimHit_gy);
  t->Branch("fSimHitGz",           &fSimHit_gz);
  t->Branch("fClosestSimHitPiD",   &fClosestSimHit_pid);
  t->Branch("fClosestSimHitEloss", &fClosestSimHit_eloss);
  t->Branch("fClosestSimHitToF",   &fClosestSimHit_tof);
  t->Branch("fClosestSimHitLx",    &fClosestSimHit_lx);
  t->Branch("fClosestSimHitLy",    &fClosestSimHit_ly);
  t->Branch("fClosestSimHitGx",    &fClosestSimHit_gx);
  t->Branch("fClosestSimHitGy",    &fClosestSimHit_gy);
  t->Branch("fClosestSimHitGz",    &fClosestSimHit_gz);
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
     bunch_id.push_back(PVI->getBunchCrossing());
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
  // for finding matched simhit
  TrackerHitAssociator associate(e,conf_);
  bRecHit_n = 0;
  fRecHit_n = 0;
  bSimHit_n = 0; 
  fSimHit_n = 0; 
  SiPixelRecHitCollection::const_iterator detUnitIt    = recHits->begin();
  SiPixelRecHitCollection::const_iterator detUnitItEnd = recHits->end();
  for ( ; detUnitIt != detUnitItEnd; detUnitIt++ ) {
      DetId detId = DetId(detUnitIt->detId()); 
      int subid = detId.subdetId();
      SiPixelRecHitCollection::DetSet::const_iterator recHitIt    = detUnitIt->begin();
      SiPixelRecHitCollection::DetSet::const_iterator recHitItEnd = detUnitIt->end();
      std::vector<PSimHit> matched;
      std::vector<PSimHit>::const_iterator closest_simhit;
      for ( ; recHitIt != recHitItEnd; ++recHitIt) {
          LocalPoint lp = recHitIt->localPosition();
          LocalError le = recHitIt->localPositionError();
          GlobalPoint gp = ( (geom.product())->idToDet(detId) )->surface().toGlobal(lp);
          if ( gp.perp() < 20 ){ // reject outer tracker
             matched.clear();
             matched = associate.associateHit(*recHitIt);
             if ( !matched.empty() ) {
                float closest = 9999.9;
                std::vector<PSimHit>::const_iterator closestit = matched.begin();
                //loop over simhits and find closest
                for ( std::vector<PSimHit>::const_iterator m = matched.begin(); m<matched.end(); m++) {
                    LocalPoint lpos = m->localPosition();  
                    GlobalPoint gpos = ( (geom.product())->idToDet(detId) )->surface().toGlobal(lpos);
                    if ( subid==PixelSubdetector::PixelBarrel ){
                         bSimHit_n ++;
                         bSimHit_subid.push_back(subid);
                         bSimHit_layer.push_back(  PXBDetId(detId).layer() ); 
                         bSimHit_ladder.push_back( PXBDetId(detId).ladder() ); 
                         bSimHit_module.push_back( PXBDetId(detId).module() ); 
                         bSimHit_pid.push_back(m->particleType());
                         bSimHit_eloss.push_back(m->energyLoss());
                         bSimHit_tof.push_back(m->timeOfFlight());
                         bSimHit_lx.push_back( lpos.x() );
                         bSimHit_ly.push_back( lpos.y() );
                         bSimHit_gx.push_back( gpos.x() );
                         bSimHit_gy.push_back( gpos.y() );
                         bSimHit_gz.push_back( gpos.z() );
                    }
                    if ( subid==PixelSubdetector::PixelEndcap ){
                         fSimHit_n ++;
                         fSimHit_subid.push_back(subid);
                         fSimHit_disk.push_back(   PXFDetId(detId).disk()   );
                         fSimHit_blade.push_back(  PXFDetId(detId).blade()  );
                         fSimHit_module.push_back( PXFDetId(detId).module() );
                         fSimHit_side.push_back(   PXFDetId(detId).side()   );
                         fSimHit_pid.push_back(m->particleType());
                         fSimHit_eloss.push_back(m->energyLoss());
                         fSimHit_tof.push_back(m->timeOfFlight());
                         fSimHit_lx.push_back( lpos.x() );
                         fSimHit_ly.push_back( lpos.y() );
                         fSimHit_gx.push_back( gpos.x() );
                         fSimHit_gy.push_back( gpos.y() );
                         fSimHit_gz.push_back( gpos.z() );
                    }
                    float x_res = fabs(lpos.x() - lp.x());
                    float y_res = fabs(lpos.y() - lp.y());
                    float dist = sqrt(x_res*x_res + y_res*y_res);
                    if ( dist < closest ) {
                          closest = dist;
                          closestit = m;
                    }
                } // end of simhit loop
                closest_simhit = closestit;
             } // end matched emtpy
             LocalPoint lpos = closest_simhit->localPosition();  
             GlobalPoint gpos = ( (geom.product())->idToDet(detId) )->surface().toGlobal(lpos);
             SiPixelRecHit::ClusterRef const& Cluster = recHitIt->cluster();
             if ( subid==PixelSubdetector::PixelBarrel ){
                bRecHit_n ++;
                bRecHit_lx.push_back( lp.x() );
                bRecHit_ly.push_back( lp.y() );
                bRecHit_lxx.push_back( le.xx() );
                bRecHit_lxy.push_back( le.xy() );
                bRecHit_lyy.push_back( le.yy() );
                bRecHit_gx.push_back( gp.x() );
                bRecHit_gy.push_back( gp.y() );
                bRecHit_gz.push_back( gp.z() );
                bRecHit_subid.push_back(subid);
                bRecHit_layer.push_back(  PXBDetId(detId).layer() ); 
                bRecHit_ladder.push_back( PXBDetId(detId).ladder() ); 
                bRecHit_module.push_back( PXBDetId(detId).module() ); 
                bCl_size.push_back(Cluster->size());
                bCl_sizeX.push_back(Cluster->sizeX());
                bCl_sizeY.push_back(Cluster->sizeY());
                bCl_charge.push_back(Cluster->charge());
                bCl_x.push_back(Cluster->x());
                bCl_y.push_back(Cluster->y());
                bClosestSimHit_pid.push_back( closest_simhit->particleType() );
                bClosestSimHit_eloss.push_back( closest_simhit->energyLoss() );
                bClosestSimHit_tof.push_back( closest_simhit->timeOfFlight() );
                bClosestSimHit_lx.push_back( lpos.x() );
                bClosestSimHit_ly.push_back( lpos.y() );
                bClosestSimHit_gx.push_back( gpos.x() );
                bClosestSimHit_gy.push_back( gpos.y() );
                bClosestSimHit_gz.push_back( gpos.z() );
             }
             if ( subid==PixelSubdetector::PixelEndcap ){
                fRecHit_n ++;
                fRecHit_lx.push_back( lp.x() );
                fRecHit_ly.push_back( lp.y() );
                fRecHit_lxx.push_back( le.xx() );
                fRecHit_lxy.push_back( le.xy() );
                fRecHit_lyy.push_back( le.yy() );
                fRecHit_gx.push_back( gp.x() );
                fRecHit_gy.push_back( gp.y() );
                fRecHit_gz.push_back( gp.z() );
                fRecHit_subid.push_back(subid);
                fRecHit_disk.push_back(   PXFDetId(detId).disk()   );
                fRecHit_blade.push_back(  PXFDetId(detId).blade()  );
                fRecHit_module.push_back( PXFDetId(detId).module() );
                fRecHit_side.push_back(   PXFDetId(detId).side()   );
                fCl_size.push_back(Cluster->size());
                fCl_sizeX.push_back(Cluster->sizeX());
                fCl_sizeY.push_back(Cluster->sizeY());
                fCl_charge.push_back(Cluster->charge());
                fCl_x.push_back(Cluster->x());
                fCl_y.push_back(Cluster->y());
                fClosestSimHit_pid.push_back( closest_simhit->particleType() );
                fClosestSimHit_eloss.push_back( closest_simhit->energyLoss() );
                fClosestSimHit_tof.push_back( closest_simhit->timeOfFlight() );
                fClosestSimHit_lx.push_back( lpos.x() );
                fClosestSimHit_ly.push_back( lpos.y() );
                fClosestSimHit_gx.push_back( gpos.x() );
                fClosestSimHit_gy.push_back( gpos.y() );
                fClosestSimHit_gz.push_back( gpos.z() );
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
            bunch_id.clear();
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
       bRecHit_subid.clear();
       bRecHit_layer.clear();
      bRecHit_ladder.clear();
      bRecHit_module.clear();
          bRecHit_lx.clear();
          bRecHit_ly.clear();
         bRecHit_lxx.clear();
         bRecHit_lxy.clear();
         bRecHit_lyy.clear();
          bRecHit_gx.clear();
          bRecHit_gy.clear();
          bRecHit_gz.clear();
            bCl_size.clear();
           bCl_sizeX.clear();
           bCl_sizeY.clear();
          bCl_charge.clear();
              bCl_x.clear();
              bCl_y.clear();
         bSimHit_pid.clear();
       bSimHit_subid.clear();
       bSimHit_layer.clear();
      bSimHit_ladder.clear();
      bSimHit_module.clear();
       bSimHit_eloss.clear();
         bSimHit_tof.clear();
          bSimHit_lx.clear();
          bSimHit_ly.clear();
          bSimHit_gx.clear();
          bSimHit_gy.clear();
          bSimHit_gz.clear();
  bClosestSimHit_pid.clear();
bClosestSimHit_eloss.clear();
  bClosestSimHit_tof.clear();
   bClosestSimHit_lx.clear();
   bClosestSimHit_ly.clear();
   bClosestSimHit_gx.clear();
   bClosestSimHit_gy.clear();
   bClosestSimHit_gz.clear();
       fRecHit_subid.clear();
        fRecHit_disk.clear();
       fRecHit_blade.clear();
        fRecHit_side.clear();
          fRecHit_lx.clear();
          fRecHit_ly.clear();
         fRecHit_lxx.clear();
         fRecHit_lxy.clear();
         fRecHit_lyy.clear();
          fRecHit_gx.clear();
          fRecHit_gy.clear();
          fRecHit_gz.clear();
            fCl_size.clear();
           fCl_sizeX.clear();
           fCl_sizeY.clear();
          fCl_charge.clear();
              fCl_x.clear();
              fCl_y.clear();
         fSimHit_pid.clear();
       fSimHit_subid.clear();
        fSimHit_disk.clear();
       fSimHit_blade.clear();
      fSimHit_module.clear();
        fSimHit_side.clear();
       fSimHit_eloss.clear();
         fSimHit_tof.clear();
          fSimHit_lx.clear();
          fSimHit_ly.clear();
          fSimHit_gx.clear();
          fSimHit_gy.clear();
          fSimHit_gz.clear();
  fClosestSimHit_pid.clear();
fClosestSimHit_eloss.clear();
  fClosestSimHit_tof.clear();
   fClosestSimHit_lx.clear();
   fClosestSimHit_ly.clear();
   fClosestSimHit_gx.clear();
   fClosestSimHit_gy.clear();
   fClosestSimHit_gz.clear();
 }

//define this as a plug-in
DEFINE_FWK_MODULE(L1PixelTrigger);
