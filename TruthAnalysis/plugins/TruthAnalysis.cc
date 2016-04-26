// -*- C++ -*-
//
// Package:    TruthAnalysis
// Class:      TruthAnalysis
// 
/**\class TruthAnalysis TruthAnalysis.cc test/TruthAnalysis/src/TruthAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Hamer,,,
//         Created:  Wed Nov  5 15:23:02 CET 2014
// $Id$
//
//

#include <iostream>
#include <iomanip>
#include <fstream>

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/EventID.h"

//#include "CommonTools/UtilsAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"


class TruthAnalysis : public edm::EDAnalyzer {
   public:
      explicit TruthAnalysis(const edm::ParameterSet&);
      ~TruthAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      TFile *fOut;
      TTree *tOut;
      
      std::vector<double> *mct_px = new std::vector<double>;
      std::vector<double> *mct_py = new std::vector<double>;
      std::vector<double> *mct_pz = new std::vector<double>;
      std::vector<double> *mct_e = new std::vector<double>;
      std::vector<double> *mct_vx = new std::vector<double>;
      std::vector<double> *mct_vy = new std::vector<double>;
      std::vector<double> *mct_vz = new std::vector<double>;
      std::vector<int> *mct_id = new std::vector<int>;
      std::vector<int> *mct_status = new std::vector<int>;
      std::vector<int> *mct_parent = new std::vector<int>;
      std::vector<std::vector<int> > *mct_daughters = new std::vector<std::vector<int> >;

      int eventNumber = 0;
      int lumiBlock = 0;
      int runNumber = 0;
      int EventCounter = 0;
      std::vector<int> pRunNumbers;
      std::vector<int> pEventNumbers;
      std::vector<int> pLumiBlocks;
      
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken;

};

//
// constructors and destructor
//
TruthAnalysis::TruthAnalysis(const edm::ParameterSet& iConfig):
genParticleToken(consumes<reco::GenParticleCollection>( iConfig.getParameter<edm::InputTag>("GenParticles" )))
{
  //now do what ever initialization is needed
  fOut = new TFile("output.root", "RECREATE");
  tOut = new TTree("output", "output");
  tOut->Branch("EventNumber", &eventNumber );
  tOut->Branch("LumiBlock", &lumiBlock );
  tOut->Branch("RunNumber", &runNumber );
  tOut->Branch("Particle_px", &mct_px );
  tOut->Branch("Particle_py", &mct_py );
  tOut->Branch("Particle_pz", &mct_pz );
  tOut->Branch("Particle_vx", &mct_vx );
  tOut->Branch("Particle_vy", &mct_vy );
  tOut->Branch("Particle_vz", &mct_vz );
  tOut->Branch("Particle_E", &mct_e );
  tOut->Branch("Particle_ID", &mct_id );
  tOut->Branch("Particle_Status", &mct_status );
  tOut->Branch("Particle_Parent", &mct_parent );
  tOut->Branch("Particle_Daughter", &mct_daughters );

  /*
  std::string pfilename = iConfig.getParameter<std::string>("PassedFileName");
  std::ifstream pfile( pfilename, std::ios::in );
  while( pfile.good() ) {
    std::string prn, plb, pen;
    pfile >> prn >> std::ws >> plb >> std::ws >> pen;
    if( pfile.eof() ) break;
    pRunNumbers.push_back( atoi( prn.c_str() ) );
    pEventNumbers.push_back( atoi( pen.c_str() ) );
    pLumiBlocks.push_back( atoi( plb.c_str() ) );
  }
  */
}


TruthAnalysis::~TruthAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   gDirectory = fOut;
   tOut -> Write();
   fOut -> Close();


}


//
// member functions
//

// ------------ method called for each event  ------------
void
TruthAnalysis::analyze(const edm::Event& e,const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  mct_px -> clear();
  mct_py -> clear();
  mct_pz -> clear();
  mct_vx -> clear();
  mct_vy -> clear();
  mct_vz -> clear();
  mct_e -> clear();
  mct_id -> clear();
  mct_status -> clear();
  mct_parent -> clear();
  mct_daughters -> clear();


  EventCounter += 1;

  edm::EventAuxiliary aux = e.eventAuxiliary();
  edm::EventID id = aux.id();
  
  eventNumber = id.event();
  runNumber = id.run();
  lumiBlock = id.luminosityBlock();

  /*
  bool procEvent = false;
  for( unsigned int ic = 0; ic < pRunNumbers.size(); ++ic ) {
    if( eventNumber == pEventNumbers.at(ic) && runNumber == pRunNumbers.at(ic) && lumiBlock == pLumiBlocks.at(ic) ) {
      procEvent = true;
      break;
    }
  }
  
  if( !procEvent ) return;
  */

  Handle<reco::GenParticleCollection> genParticles;
  e.getByToken( genParticleToken, genParticles);

  for(size_t i = 0; i < genParticles->size(); ++i) {
       const GenParticle & p = (*genParticles)[i];
       mct_px->push_back( p.px() );
       mct_py->push_back( p.py() );
       mct_pz->push_back( p.pz() );
       mct_e->push_back( p.energy() );
       mct_vx->push_back( p.vx() );
       mct_vy->push_back( p.vy() );
       mct_vz->push_back( p.vz() );
       mct_id->push_back( p.pdgId() );
       mct_status->push_back( p.status() );
       

        const Candidate * mom = p.mother();
        if( !mom ) mct_parent->push_back( -2 );
        else {
          int imother = -1;
          for( size_t j = 0; j < genParticles->size(); ++j ) {
            const GenParticle &p2 = (*genParticles)[j];
            // now check if p2 is the mom
            // the 4-momentum should be the same as four the mom
            // one of the children should be the particle
            if( fabs( p2.px() - mom->px() ) < 1.e-10 && 
                fabs( p2.py() - mom->py() ) < 1.e-10 &&
                fabs( p2.pz() - mom->pz() ) < 1.e-10 &&
                p2.pdgId() == mom->pdgId() &&
                p2.status() == mom->status() ) {
                  bool isMother = false;
                  int targetId = p.pdgId();
                  // now i have a particle with the same four vector as the mother
                  // check if the target id is amongst the daughter's 
                  // if yes: it's the mother
                  size_t nD = p2.numberOfDaughters();
                  for( size_t j = 0; j < nD; ++j ) {
                    if( p2.daughter(j)->pdgId() == targetId && p2.daughter(j)->status() == p.status() ) {
                      isMother = true;
                      break;
                    }
                  }
                  if( isMother ) {
                    imother = j;
                    break;
                  }
                }
          }
          mct_parent->push_back( imother );
       }
      
       std::vector<int> daughters;
       size_t n = p.numberOfDaughters();
       for(size_t j = 0; j < n; ++ j) {
          const Candidate * d = p.daughter( j );
          int iDaughter = -1;
          for( size_t k = 0; k < genParticles->size(); ++k ) {
            const GenParticle &p2 = (*genParticles)[k];
            bool isDaughterCandidate = false;
            if( p2.mother() ) {
              if( p2.mother()->pdgId() == p.pdgId() && p2.mother()->status() == p.status() ) isDaughterCandidate = true;
            }
            if( !isDaughterCandidate ) continue;
            if( fabs( p2.px() - d->px() ) < 1.e-10 &&
                fabs( p2.py() - d->py() ) < 1.e-10 &&
                fabs( p2.pz() - d->pz() ) < 1.e-10 &&
                p2.pdgId() == d->pdgId() &&
                p2.status() == d->status() ) {
                iDaughter = k;
                break;
            }
          }
          daughters.push_back( iDaughter );
       }
       mct_daughters->push_back( daughters );
  }
  
  tOut->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TruthAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TruthAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TruthAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TruthAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TruthAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TruthAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TruthAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TruthAnalysis);
