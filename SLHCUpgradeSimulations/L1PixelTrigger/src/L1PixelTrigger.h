#ifndef L1PixelTrigger_h
#define L1PixelTrigger_h

// system include files
#include <cmath>


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

// Tracker data formats
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// --- L1 Egamma dataFormats
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

// SiPixelRecHit dataFormat
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// --- GenParticles and SimHits
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

// --- Beam Spot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// ROOT includes
#include "TTree.h"

class L1PixelTrigger  : public edm::EDAnalyzer {
 public:

      explicit L1PixelTrigger(const edm::ParameterSet&);
      ~L1PixelTrigger();

 private:

      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      int                                      run;
      int                                    event;
      int                                  bunch_n;
      std::vector<int>                    bunch_id;
      std::vector<int>                      pileup;
      float                            beamspot_x0;
      float                            beamspot_y0;
      float                            beamspot_z0;
      float                       beamspot_x0Error;
      float                       beamspot_y0Error;
      float                       beamspot_z0Error;
      float                            beam_widthX;
      float                            beam_widthY;
      float                            beam_sigmaZ;
      float                       beam_widthXError;
      float                       beam_widthYError;
      float                       beam_sigmaZError;
      int                                genpart_n;
      std::vector<float>                genpart_et;
      std::vector<float>                genpart_pt;
      std::vector<float>               genpart_eta;
      std::vector<float>               genpart_phi;
      std::vector<int>              genpart_charge;
      std::vector<int>                  genpart_id;
      int                                egammaC_n;
      std::vector<float>                 egammaC_e;
      std::vector<float>                egammaC_et;
      std::vector<float>               egammaC_eta;
      std::vector<float>               egammaC_phi;
      std::vector<float>                egammaC_gx;
      std::vector<float>                egammaC_gy;
      std::vector<float>                egammaC_gz;
      std::vector<int>              egammaC_charge;
      int                                 egamma_n;
      std::vector<float>                  egamma_e;
      std::vector<float>                 egamma_et;
      std::vector<float>                egamma_eta;
      std::vector<float>                egamma_phi;
      std::vector<float>                 egamma_gx;
      std::vector<float>                 egamma_gy;
      std::vector<float>                 egamma_gz;
      std::vector<int>               egamma_charge;
      int                                bRecHit_n;
      std::vector<int>               bRecHit_subid;
      std::vector<int>               bRecHit_layer;
      std::vector<int>              bRecHit_ladder;
      std::vector<int>              bRecHit_module;
      std::vector<float>                bRecHit_lx;
      std::vector<float>                bRecHit_ly;
      std::vector<float>               bRecHit_lxx;
      std::vector<float>               bRecHit_lxy;
      std::vector<float>               bRecHit_lyy;
      std::vector<float>                bRecHit_gx;
      std::vector<float>                bRecHit_gy;
      std::vector<float>                bRecHit_gz;
      std::vector<int>                    bCl_size;
      std::vector<int>                   bCl_sizeX;
      std::vector<int>                   bCl_sizeY;
      std::vector<float>                bCl_charge;  
      std::vector<float>                     bCl_x;  
      std::vector<float>                     bCl_y;  
      int                                bSimHit_n;
      std::vector<int>               bSimHit_subid;
      std::vector<int>               bSimHit_layer;
      std::vector<int>              bSimHit_ladder;
      std::vector<int>              bSimHit_module;
      std::vector<int>                 bSimHit_pid;
      std::vector<float>             bSimHit_eloss;
      std::vector<float>               bSimHit_tof;
      std::vector<float>                bSimHit_lx;
      std::vector<float>                bSimHit_ly;
      std::vector<float>                bSimHit_gx;
      std::vector<float>                bSimHit_gy;
      std::vector<float>                bSimHit_gz;
      std::vector<int>          bClosestSimHit_pid;
      std::vector<float>      bClosestSimHit_eloss;
      std::vector<float>        bClosestSimHit_tof;
      std::vector<float>         bClosestSimHit_lx;
      std::vector<float>         bClosestSimHit_ly;
      std::vector<float>         bClosestSimHit_gx;
      std::vector<float>         bClosestSimHit_gy;
      std::vector<float>         bClosestSimHit_gz;
      int                                fRecHit_n;
      std::vector<int>               fRecHit_subid;
      std::vector<int>                fRecHit_disk;
      std::vector<int>               fRecHit_blade;
      std::vector<int>              fRecHit_module;
      std::vector<int>                fRecHit_side;
      std::vector<float>                fRecHit_lx;
      std::vector<float>                fRecHit_ly;
      std::vector<float>               fRecHit_lxx;
      std::vector<float>               fRecHit_lxy;
      std::vector<float>               fRecHit_lyy;
      std::vector<float>                fRecHit_gx;
      std::vector<float>                fRecHit_gy;
      std::vector<float>                fRecHit_gz;
      std::vector<int>                    fCl_size;
      std::vector<int>                   fCl_sizeX;
      std::vector<int>                   fCl_sizeY;
      std::vector<float>                fCl_charge;  
      std::vector<float>                     fCl_x;  
      std::vector<float>                     fCl_y;  
      int                                fSimHit_n;
      std::vector<int>               fSimHit_subid;
      std::vector<int>                fSimHit_disk;
      std::vector<int>               fSimHit_blade;
      std::vector<int>              fSimHit_module;
      std::vector<int>                fSimHit_side;
      std::vector<int>                 fSimHit_pid;
      std::vector<float>             fSimHit_eloss;
      std::vector<float>               fSimHit_tof;
      std::vector<float>                fSimHit_lx;
      std::vector<float>                fSimHit_ly;
      std::vector<float>                fSimHit_gx;
      std::vector<float>                fSimHit_gy;
      std::vector<float>                fSimHit_gz;
      std::vector<int>          fClosestSimHit_pid;
      std::vector<float>      fClosestSimHit_eloss;
      std::vector<float>        fClosestSimHit_tof;
      std::vector<float>         fClosestSimHit_lx;
      std::vector<float>         fClosestSimHit_ly;
      std::vector<float>         fClosestSimHit_gx;
      std::vector<float>         fClosestSimHit_gy;
      std::vector<float>         fClosestSimHit_gz;
      void InitializeVectors();
      TTree* t;
      edm::InputTag eGamma_tag;
      edm::InputTag eGammaCrystal_tag;
      edm::ParameterSet conf_;
      GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
};

#endif
