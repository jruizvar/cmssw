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

// --- GenParticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


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
      int                                   bHit_n;
      std::vector<int>                  bHit_subid;
      std::vector<int>                  bHit_layer;
      std::vector<int>                 bHit_ladder;
      std::vector<int>                 bHit_module;
      std::vector<int>                    bCl_size;
      std::vector<float>                   bHit_gx;
      std::vector<float>                   bHit_gy;
      std::vector<float>                   bHit_gz;
      int                                   fHit_n;
      std::vector<int>                  fHit_subid;
      std::vector<int>                   fHit_disk;
      std::vector<int>                  fHit_blade;
      std::vector<int>                 fHit_module;
      std::vector<int>                   fHit_side;
      std::vector<int>                    fCl_size;
      std::vector<float>                   fHit_gx;
      std::vector<float>                   fHit_gy;
      std::vector<float>                   fHit_gz;
      void InitializeVectors();
      TTree* t;
      edm::InputTag eGamma_tag;
      edm::InputTag eGammaCrystal_tag;
      GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
};

#endif
