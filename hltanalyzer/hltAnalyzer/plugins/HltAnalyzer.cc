// -*- C++ -*-
//
// Package:    HLTAnalyzer/HltAnalyzer
// Class:      HltAnalyzer
// 
/**\class HltAnalyzer HltAnalyzer.cc HLTAnalyzer/HltAnalyzer/plugins/HltAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jose Cupertino Ruiz Vargas
//         Created:  Sun, 09 Mar 2014 18:38:18 GMT
//
//


// system include files
#include <memory>
#include <map>
#include <string>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "TH1.h"
#include "TEfficiency.h"
#include "TObjString.h"
#include "TString.h"
#include "TObject.h"
#include "TPRegexp.h"

//
// class declaration
//

class HltAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HltAnalyzer(const edm::ParameterSet&);
      ~HltAnalyzer();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
      virtual void endJob() ;

      void beginRun(const edm::Run& run, const edm::EventSetup& c);
      void endRun(const edm::Run& run, const edm::EventSetup& c);


      // ----------member data ---------------------------
      HLTConfigProvider hltConfig_;
      edm::EDGetTokenT<edm::TriggerResults> hlt_token_;
      std::string processName_;
      std::vector<std::string> hltPaths_;
      std::vector<std::string> hltNames;
      edm::InputTag eles_tag_;
      edm::InputTag e1_tag_;
      edm::InputTag e2_tag_;
      std::map<std::string,TEfficiency*> effContainer;
      TH1I *elesMult;
      TH1D *elesPt;
      TH1D *elesEta;
      TH1D *elesPhi;
      TH1I *elesCharge;
      TEfficiency *eff_e1Pt;
      TEfficiency *eff_e1Eta;
      TEfficiency *eff_e2Pt;
      TEfficiency *eff_e2Eta;
      TEfficiency *eff_eedR;
      TEfficiency *eff_eedPhi;
      TEfficiency *eff_eedEta;
      TEfficiency *eff_eePt;
      TEfficiency *eff_eeEta;
      TEfficiency *eff_eeMass;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HltAnalyzer::HltAnalyzer(const edm::ParameterSet& iConfig):
   hlt_token_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltLabel"))),
   processName_(iConfig.getParameter<std::string>("procName")),
   hltPaths_(iConfig.getParameter<std::vector<std::string>>("hltPaths")),
   eles_tag_(iConfig.getParameter<edm::InputTag>("elesTag")),
   e1_tag_(iConfig.getParameter<edm::InputTag>("e1Tag")),
   e2_tag_(iConfig.getParameter<edm::InputTag>("e2Tag"))
{
   //now do what ever initialization is needed
}

HltAnalyzer::~HltAnalyzer()
{
}

// ------------ method called for each event  ------------
void
HltAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   // ------ analize all electrons----------//
   Handle<reco::CandidateView> allEles;
   iEvent.getByLabel(eles_tag_, allEles);
   elesMult->Fill(allEles->size());
   reco::CandidateView::const_iterator it=allEles->begin();
   for (; it!=allEles->end(); ++it) {
       elesPt    ->Fill(it->pt()    );
       elesEta   ->Fill(it->eta()   );
       elesPhi   ->Fill(it->phi()   );
       elesCharge->Fill(it->charge());
   }
   // ------ handle positive and negative electrons----------//
   Handle<reco::CandidateView> posEle;
   Handle<reco::CandidateView> negEle;
   iEvent.getByLabel(e1_tag_, posEle);
   iEvent.getByLabel(e2_tag_, negEle);
   reco::CandidateView::const_iterator e1=posEle->begin();
   reco::CandidateView::const_iterator e2=negEle->begin();
   // ------ deltaPhi, deltaR, and Invariant mass -------//
   const reco::Candidate::LorentzVector & e1p = e1->p4();
   const reco::Candidate::LorentzVector & e2p = e2->p4();
   double eedR   = reco::deltaR(e1p,e2p);
   double eedPhi = reco::deltaPhi(e1->phi(),e2->phi());
   double eedEta = std::sqrt(eedR*eedR-eedPhi*eedPhi);
   double eePt   = (e1p + e2p).pt();
   double eeEta  = (e1p + e2p).eta();
   double eeMass = (e1p + e2p).mass();
   // ------ analize trigger results ----------//
   Handle<TriggerResults> trigRes;
   iEvent.getByToken(hlt_token_, trigRes);
   for (size_t i=0; i< hltNames.size(); i++){
       const unsigned int path_index = hltConfig_.triggerIndex(hltNames[i]);
       bool path_bit = trigRes->accept(path_index);
       effContainer[ "eff_e1Pt_"   +hltNames[i] ] -> Fill( path_bit, e1->pt()  );
       effContainer[ "eff_e1Eta_"  +hltNames[i] ] -> Fill( path_bit, e1->eta() );
       effContainer[ "eff_e2Pt_"   +hltNames[i] ] -> Fill( path_bit, e2->pt()  );
       effContainer[ "eff_e2Eta_"  +hltNames[i] ] -> Fill( path_bit, e2->eta() );
       effContainer[ "eff_eedR_"   +hltNames[i] ] -> Fill( path_bit, eedR      );
       effContainer[ "eff_eedPhi_" +hltNames[i] ] -> Fill( path_bit, eedPhi    );
       effContainer[ "eff_eedEta_" +hltNames[i] ] -> Fill( path_bit, eedEta    );
       effContainer[ "eff_eePt_"   +hltNames[i] ] -> Fill( path_bit, eePt      );
       effContainer[ "eff_eeEta_"  +hltNames[i] ] -> Fill( path_bit, eeEta     );
       effContainer[ "eff_eeMass_" +hltNames[i] ] -> Fill( path_bit, eeMass    );
   }// close trigger results
}// close analyze method

// ------------ method called once each job just before starting event loop  ------------
void HltAnalyzer::beginJob() { }

void HltAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& c)
{
   bool changed;
   if (!hltConfig_.init(run,c,processName_,changed)) {
     edm::LogError("HltAnalysis") << "Initialization of HLTConfigProvider failed!!";
     return;
   }
   for (size_t i = 0; i < hltPaths_.size(); i++) {
       TPRegexp pattern(hltPaths_[i]);
       for ( size_t j = 0; j < hltConfig_.triggerNames().size(); j++ )
         if ( TString(hltConfig_.triggerNames()[j]).Contains(pattern) )
           hltNames.push_back( hltConfig_.triggerNames()[j] );
   }
   edm::Service<TFileService> fs;
   elesMult  = fs->make<TH1I>("elesMult"  ,"Electrons multiplicity;number of electrons", 6, 0, 6);
   elesPt    = fs->make<TH1D>("elesPt"    ,"Electrons transverse momentum;p_{T} [GeV]", 60, 0, 3000);
   elesEta   = fs->make<TH1D>("elesEta"   ,"Electrons pseudo-rapidity;#eta", 80, -4., 4.);
   elesPhi   = fs->make<TH1D>("elesPhi"   ,"Electrons azimuth;#phi", 80, -4., 4.);
   elesCharge= fs->make<TH1I>("elesCharge","Electrons charge;charge",5, -2, 3);
   for (size_t i = 0; i < hltNames.size(); i++) {
       eff_e1Pt   = fs->make<TEfficiency>(Form("eff_e1Pt_%s"  ,hltNames[i].c_str()), Form("MC Efficiency of %s;p_{T}(e plus) [GeV];Efficiency",hltNames[i].c_str()),  60, 0., 3000.);
       eff_e1Eta  = fs->make<TEfficiency>(Form("eff_e1Eta_%s" ,hltNames[i].c_str()), Form("MC Efficiency of %s;#eta (e plus);Efficiency"      ,hltNames[i].c_str()),  40, -4.,   4.);
       effContainer[ "eff_e1Pt_"  + hltNames[i] ] = eff_e1Pt;
       effContainer[ "eff_e1Eta_" + hltNames[i] ] = eff_e1Eta;
       eff_e2Pt   = fs->make<TEfficiency>(Form("eff_e2Pt_%s"  ,hltNames[i].c_str()), Form("MC Efficiency of %s;p_{T}(e minus) [GeV];Efficiency",hltNames[i].c_str()),  60, 0., 3000.);
       eff_e2Eta  = fs->make<TEfficiency>(Form("eff_e2Eta_%s" ,hltNames[i].c_str()), Form("MC Efficiency of %s;#eta (e minus);Efficiency"      ,hltNames[i].c_str()),  40,-4.,    4.);
       effContainer[ "eff_e2Pt_"  + hltNames[i] ] = eff_e2Pt;
       effContainer[ "eff_e2Eta_" + hltNames[i] ] = eff_e2Eta;
       eff_eedR   = fs->make<TEfficiency>(Form("eff_eedR_%s"  ,hltNames[i].c_str()), Form("MC Efficiency of %s;dR(ee);Efficiency"         ,hltNames[i].c_str()),  60,  0.,   1.5);
       eff_eedPhi = fs->make<TEfficiency>(Form("eff_eedPhi_%s",hltNames[i].c_str()), Form("MC Efficiency of %s;dPhi(ee);Efficiency"       ,hltNames[i].c_str()),  60,  0.,   1.5);
       eff_eedEta = fs->make<TEfficiency>(Form("eff_eedEta_%s",hltNames[i].c_str()), Form("MC Efficiency of %s;dEta(ee);Efficiency"       ,hltNames[i].c_str()),  60,  0.,   1.5);
       eff_eePt   = fs->make<TEfficiency>(Form("eff_eePt_%s"  ,hltNames[i].c_str()), Form("MC Efficiency of %s;p_{T}(ee) [GeV];Efficiency",hltNames[i].c_str()),  60,  0., 3000.);
       eff_eeEta  = fs->make<TEfficiency>(Form("eff_eeEta_%s" ,hltNames[i].c_str()), Form("MC Efficiency of %s;#eta (ee) ;Efficiency"     ,hltNames[i].c_str()),  40, -4.,    4.);
       eff_eeMass = fs->make<TEfficiency>(Form("eff_eeMass_%s",hltNames[i].c_str()), Form("MC Efficiency of %s;m(ee) [GeV];Efficiency"    ,hltNames[i].c_str()),  30, 70.,  130.);
       effContainer[ "eff_eedR_"   + hltNames[i] ] = eff_eedR;
       effContainer[ "eff_eedPhi_" + hltNames[i] ] = eff_eedPhi;
       effContainer[ "eff_eedEta_" + hltNames[i] ] = eff_eedEta;
       effContainer[ "eff_eePt_"   + hltNames[i] ] = eff_eePt;
       effContainer[ "eff_eeEta_"  + hltNames[i] ] = eff_eeEta;
       effContainer[ "eff_eeMass_" + hltNames[i] ] = eff_eeMass;
   }

}// close beginRun

// ------------ method called once each job just after ending the event loop  ------------
void HltAnalyzer::endRun(const edm::Run& run, const edm::EventSetup& c) { }

void HltAnalyzer::endJob() { }

//define this as a plug-in
DEFINE_FWK_MODULE(HltAnalyzer);
