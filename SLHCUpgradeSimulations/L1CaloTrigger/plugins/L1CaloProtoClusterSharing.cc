
#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/L1CaloAlgoBase.h"

#include "SimDataFormats/SLHC/interface/L1CaloClusterWithSeed.h"
#include "SimDataFormats/SLHC/interface/L1CaloClusterWithSeedFwd.h"

#include "SLHCUpgradeSimulations/L1CaloTrigger/interface/TriggerTowerGeometry.h"

#include "SimDataFormats/SLHC/interface/L1TowerNav.h"


class L1CaloProtoClusterSharing:public L1CaloAlgoBase < l1slhc::L1CaloClusterWithSeedCollection , l1slhc::L1CaloClusterWithSeedCollection  > 
{
    public:
        L1CaloProtoClusterSharing( const edm::ParameterSet & );
        ~L1CaloProtoClusterSharing(  );

        void initialize(  );

        void algorithm( const int &, const int & );

    private:
        int mHoECutEB, mHoECutEE;


};

L1CaloProtoClusterSharing::L1CaloProtoClusterSharing( const edm::ParameterSet & aConfig ):
L1CaloAlgoBase < l1slhc::L1CaloClusterWithSeedCollection , l1slhc::L1CaloClusterWithSeedCollection > ( aConfig )
{
}

L1CaloProtoClusterSharing::~L1CaloProtoClusterSharing(  )
{
}

/*
void L1CaloProtoClusterSharing::initialize(  )
{
}
*/

void L1CaloProtoClusterSharing::initialize(  )
{
    mHoECutEB = 40; // 0-1000 -> 0-1 : 40 := 0.04
    mHoECutEE = 15; // := 0.015
}


void L1CaloProtoClusterSharing::algorithm( const int &aEta, const int &aPhi )
{

    // Look if there is a cluster here
    l1slhc::L1CaloClusterWithSeedCollection::const_iterator lClusterItr = fetch( aEta, aPhi );
    if ( lClusterItr != mInputCollection->end(  ) )
    {

        l1slhc::L1CaloClusterWithSeed lSharedCluster( *lClusterItr );
        // loop over cluster constituents
        for ( int lTowerEta = aEta-1; lTowerEta <= aEta + 1; ++lTowerEta )
        {
            for ( int lTowerPhi = aPhi-1; lTowerPhi <= aPhi + 1; ++lTowerPhi )
            {
                if(lTowerEta==aEta && lTowerPhi==aPhi)
                {
                    continue;
                }
                // look at clusters around constituents, skipping the one at (aEta,aPhi)
                // And find max and 2nd max neighbor cluster
                int maxE = 0;
                int secondMaxE = 0;
                l1slhc::L1CaloClusterWithSeedCollection::const_iterator maxCluster = mInputCollection->end(  );
                l1slhc::L1CaloClusterWithSeedCollection::const_iterator secondMaxCluster = mInputCollection->end(  );
                for(int lClusterEta = lTowerEta-1; lClusterEta <= lTowerEta+1; ++lClusterEta)
                {
                    for(int lClusterPhi = lTowerPhi-1; lClusterPhi <= lTowerPhi+1; ++lClusterPhi)
                    {
                        if((lClusterEta==aEta && lClusterPhi==aPhi) || (lClusterEta==lTowerEta && lClusterPhi==lTowerPhi) )
                        {
                            continue;
                        }
                        l1slhc::L1CaloClusterWithSeedCollection::const_iterator lNeighborItr = fetch( lClusterEta, lClusterPhi );
                        if ( lNeighborItr != mInputCollection->end(  ) )
                        {
                            if(lNeighborItr->EmEt()>secondMaxE)
                            {
                                if(lNeighborItr->EmEt()>maxE)
                                {
                                    secondMaxE       = maxE;
                                    secondMaxCluster = maxCluster;
                                    maxE             = lNeighborItr->EmEt();
                                    maxCluster       = lNeighborItr;
                                }
                                else
                                {
                                    secondMaxE       = lNeighborItr->EmEt();
                                    secondMaxCluster = lNeighborItr;
                                }
                            }

                        }
                    }
                }
                // Share energy depending on the rank of the cluster
                if(secondMaxE>lClusterItr->EmEt()) // 3rd or more -> give all the tower energy
                {
                    lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 5);
                }
                else if(secondMaxE==lClusterItr->EmEt() && 
                        (secondMaxCluster->iEta() > lClusterItr->iEta() || (secondMaxCluster->iEta()==lClusterItr->iEta() && secondMaxCluster->iPhi()>lClusterItr->iPhi()) )
                       ) // 2nd ex-aequo, favor cluster with biggest ieta, and biggest iphi -> give all the tower energy
                {
                    lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 5);
                }
                else // 2nd or 1st
                {
                    if(lClusterItr->EmEt() >= 4*maxE) // -> keep all the tower energy
                    {
                        lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 0);
                    }
                    else if(lClusterItr->EmEt() >= 2*maxE) // -> keep 3/4 of the tower energy
                    {
                        lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 1);
                    }
                    else if(lClusterItr->EmEt() > maxE) // -> keep 1/2+ of the tower energy (+1 if tower energy is odd)
                    {
                        lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 2);
                    }
                    else if(lClusterItr->EmEt() == maxE && 
                            (maxCluster->iEta() < lClusterItr->iEta() || (maxCluster->iEta()==lClusterItr->iEta() && maxCluster->iPhi()<lClusterItr->iPhi()))
                           ) // 1st ex-aequo, favor cluster with biggest ieta, and biggest iphi
                    {
                        lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 2); // -> keep 1/2+ of the tower energy (+1 if tower energy is odd)
                    }
                    else if(2*lClusterItr->EmEt() >= maxE) // -> keep 1/2- of the tower energy (-1 if tower energy is odd)
                    {
                        lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 3);
                    }
                    else if(4*lClusterItr->EmEt() >= maxE) // -> keep 1/4 of the tower energy
                    {
                        lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 4);
                    }
                    else // -> give all the tower energy
                    {
                        lSharedCluster.shareConstituent(lTowerEta-aEta, lTowerPhi-aPhi, 5);
                    }
                }
            }
        }


        // Calculate Electron Cut and Save it in the Cluster
        int lElectronValue = ( int )( 1000. * ( ( double )lSharedCluster.HadEt() ) / ( ( double )lSharedCluster.EmEt() ) );
        lSharedCluster.setEGammaValue( lElectronValue );


        // Electron Bit Decision
        bool egammaBitEB =  (abs(lSharedCluster.iEta())<=16 && lElectronValue<=mHoECutEB); 
        bool egammaBitEE =  (abs(lSharedCluster.iEta())>16 && lElectronValue<=mHoECutEE);


        // FineGrain bit already set in the initialization of the cluster using the FG of the seed tower
        lSharedCluster.setEGamma( (egammaBitEB || egammaBitEE) );

        int lIndex = mCaloTriggerSetup->getBin( aEta, aPhi );
        std::pair < int, int >lEtaPhi = mCaloTriggerSetup->getTowerEtaPhi( lIndex );
        mOutputCollection->insert( lEtaPhi.first , lEtaPhi.second , lSharedCluster );
    }
}



DEFINE_EDM_PLUGIN (edm::MakerPluginFactory,edm::WorkerMaker<L1CaloProtoClusterSharing>,"L1CaloProtoClusterSharing");
DEFINE_FWK_PSET_DESC_FILLER(L1CaloProtoClusterSharing);

