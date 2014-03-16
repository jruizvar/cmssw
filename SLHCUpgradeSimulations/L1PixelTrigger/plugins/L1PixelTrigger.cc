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
  for ( PVI = puinfo->begin(); PVI != puinfo->end(); ++PVI ) {
      pileup.push_back(PVI->getPU_NumInteractions());
  }
  ///////////////////////////////////////////////////////////
  // Gen Particles  
  //////////////////////////////////////////////////////////
  edm::Handle< reco::GenParticleCollection > genParticles;
  e.getByLabel( "genParticles", genParticles );
  genpart_n = genParticles->size() ;
  for ( int i = 0; i < genpart_n; ++i ){
      const reco::GenParticle & genParticle = genParticles->at(i);
      genpart_e.push_back(      genParticle.energy() );
      genpart_pt.push_back(     genParticle.pt()     );
      genpart_eta.push_back(    genParticle.eta()    );
      genpart_phi.push_back(    genParticle.phi()    );
      genpart_charge.push_back( genParticle.charge() );
      genpart_id.push_back(     genParticle.pdgId()  );
  }
  ///////////////////////////////////////////////////////////
  // Hcal Trigger Primitives 
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1JetParticleCollection > tau;
  e.getByLabel( tau_tag, tau );
  tau_n = tau->size();
  l1extra::L1JetParticleCollection::const_iterator tauIt = tau->begin();
  for ( ; tauIt!=tau->end(); ++tauIt ){
      tau_e.push_back(   tauIt->energy() );
      tau_et.push_back(  tauIt->et()     );
      tau_eta.push_back( tauIt->eta()    );
      tau_phi.push_back( tauIt->phi()    );
  }
  ///////////////////////////////////////////////////////////
  // Ecal Trigger Primitives
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > egamma;
  e.getByLabel( egamma_tag, egamma );
  egamma_n = egamma->size();
  l1extra::L1EmParticleCollection::const_iterator egIt = egamma->begin();
  for ( ; egIt!=egamma->end(); ++egIt ){
      egamma_e.push_back(   egIt->energy() );
      egamma_et.push_back(  egIt->et()     );
      egamma_eta.push_back( egIt->eta()    );
      egamma_phi.push_back( egIt->phi()    );
  }
  //////////////////////////////////////////////////////////
  // Geometry record 
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
          if ( gp.perp() < 20 ) { // drop outer tracker 
                if ( subid==PixelSubdetector::PixelEndcap ){
                     fHit_n ++;
                     fHit_disk.push_back(  PXFDetId(detId).disk()  );
                     fHit_blade.push_back( PXFDetId(detId).blade() );
                     fHit_side.push_back(  PXFDetId(detId).side()  );
                     fHit_gx.push_back(    gp.x() );
                     fHit_gy.push_back(    gp.y() );
                     fHit_gz.push_back(    gp.z() );
                     fHit_eta.push_back(   gp.eta() );
                     fHit_phi.push_back(   gp.phi() );
                     fCl_size.push_back(   Cluster->size()  );
                     fCl_sizex.push_back(  Cluster->sizeX() );
                     fCl_sizey.push_back(  Cluster->sizeY() );
                }
                if ( subid==PixelSubdetector::PixelBarrel ){
                     bHit_n ++;
                     bHit_layer.push_back(  PXBDetId(detId).layer() ); 
                     bHit_ladder.push_back( PXBDetId(detId).ladder() ); 
                     bHit_gx.push_back(     gp.x() );
                     bHit_gy.push_back(     gp.y() );
                     bHit_gz.push_back(     gp.z() );
                     bHit_eta.push_back(    gp.eta() );
                     bHit_phi.push_back(    gp.phi() );
                     bCl_size.push_back(    Cluster->size()  );
                     bCl_sizex.push_back(   Cluster->sizeX() );
                     bCl_sizey.push_back(   Cluster->sizeY() );
                }
          } 
      } // close recHits loop
  } // close detUnits loop

  //////////////////////////////////////////////////////////
  // Fill Tree 
  //////////////////////////////////////////////////////////
  t->Fill();

}// close L1PixelTrigger::analyze

void L1PixelTrigger::InitializeVectors()
{
             pileup.clear();
          genpart_e.clear();
         genpart_pt.clear();
        genpart_eta.clear();
        genpart_phi.clear();
     genpart_charge.clear();
         genpart_id.clear();
              tau_e.clear();
             tau_et.clear();
            tau_eta.clear();
            tau_phi.clear();
           egamma_e.clear();
          egamma_et.clear();
         egamma_eta.clear();
         egamma_phi.clear();
          fHit_disk.clear();
         fHit_blade.clear();
          fHit_side.clear();
            fHit_gx.clear();
            fHit_gy.clear();
            fHit_gz.clear();
           fHit_eta.clear();
           fHit_phi.clear();
           fCl_size.clear();
          fCl_sizex.clear();
          fCl_sizey.clear();
         bHit_layer.clear();
        bHit_ladder.clear();
            bHit_gx.clear();
            bHit_gy.clear();
            bHit_gz.clear();
           bCl_size.clear();
          bCl_sizex.clear();
          bCl_sizey.clear();
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1PixelTrigger);
