// ===========================================================================
// 
//       Filename:  DelHATS.cc
// 
//    Description:  An example code for Delphes exercise of JetMET performance
//    1. Basic Jet distribution comparison
//    2. MET and MET resolution 
//    3. Jet response simply using TProfile
//    4. Jet resolution with Jet (Pt>30) for different eta range (0-2.5),(2.5-4)
// 
//        Version:  1.0
//        Created:  02/11/2014 01:44:39 PM
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu, John Stupak
//        Company:  HATS@LPC
// 
// ===========================================================================

// Classes from STL
#include <cstdlib>
#include <iostream>  
#include <string>
#include <vector>
#include <sstream>
#include <cassert>

// Classes from ROOT
#include "TH1.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TProfile.h"

// Classes from Delphes
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "classes/DelphesClasses.h"

#include "../interface/JetCorrectorParameters.h"
#include "../interface/FactorizedJetCorrector.h"

double correction(Jet &iJet, FactorizedJetCorrector *iJetCorr, const TClonesArray *branchRho);

// ===  FUNCTION  ============================================================
//         Name:  main
//  Description:  Wrap up of Example1.C
// ===========================================================================
int main ( int argc, char *argv[] )
{
  if(argc != 3)
  {
    std::cout << " Usage: DelHATS input_file output_file" << std::endl;
    std::cout << " input_file - input file in ROOT format ('Delphes' tree)," << std::endl;
    std::cout << " output_file - output file in ROOT format" << std::endl;
    return EXIT_FAILURE;
  }

  // Getting the input filename
  const std::string inputFile_name  = argv[1];
  const std::string outputFile_name = argv[2];

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile_name.c_str());

  // Create the output file
  TFile outputfile(outputFile_name.c_str(), "RECREATE");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchRawJet         = treeReader->UseBranch("RawJet");
  TClonesArray *branchRawJetNoPU     = treeReader->UseBranch("RawJetNoPU");
  TClonesArray *branchGenJet         = treeReader->UseBranch("GenJet");
  TClonesArray *branchJet            = treeReader->UseBranch("Jet");
  TClonesArray *branchPileUpJetIDJet = treeReader->UseBranch("PileUpJetIDJet");
  TClonesArray *branchPuppiJet       = treeReader->UseBranch("PuppiJet");

  TClonesArray *branchElectron   = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon       = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton     = treeReader->UseBranch("Photon");
  TClonesArray *branchMet        = treeReader->UseBranch("MissingET");
  TClonesArray *branchParticle   = treeReader->UseBranch("Particle");
  TClonesArray *branchRho        = treeReader->UseBranch("Rho");
  
  // Book histograms
//----------------------------------------------------------------------------
//  Example 
//---------------------------------------------------------------------------
  TH1 *histRawJetPT = new TH1F("rawjet_pt", "Rawjet P_{T}", 200, 0.0, 200);
  TH1 *histRawJetPTCorr = new TH1F("rawjet_ptCorr", "Rawjet corrected P_{T}", 200, 0.0, 200);
  TH1 *histRawJetNoPUPT = new TH1F("rawjetNoPU_pt", "Rawjet NoPU P_{T}", 200, 0.0, 200);
  TH1 *histGenJetPT = new TH1F("genjet_pt", "Genjet P_{T}", 200, 0.0, 200);
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 200, 0.0, 200);
  TH1 *histPUIDJetPT = new TH1F("PUIDjet_pt", "PUIDjet P_{T}", 200, 0.0, 200);
  TH1 *histPUPPIJetPT = new TH1F("PUPPIjet_pt", "PUPPIjet P_{T}", 200, 0.0, 200);
  
//----------------------------------------------------------------------------
//  JEC
//----------------------------------------------------------------------------

  JetCorrectorParameters *L1JetPar;
  JetCorrectorParameters *L2JetPar;
  JetCorrectorParameters *L3JetPar;
  std::vector<JetCorrectorParameters> vPar;

  std::string path = "../data/"; 
  //std::string era = "Delphes_V1_MC";
  //std::string era = "PhaseI_50PU_V1_MC";
  std::string era = "PhaseII_140PU_V1_MC";
  std::string alias = "RawJet";
  L1JetPar = new JetCorrectorParameters(path + era + "_L1FastJet_"    + alias + ".txt");
  L2JetPar = new JetCorrectorParameters(path + era + "_L2Relative_"   + alias + ".txt");
  L3JetPar = new JetCorrectorParameters(path + era + "_L3Absolute_"   + alias + ".txt");
  vPar.push_back(*L1JetPar);
  vPar.push_back(*L2JetPar);
  vPar.push_back(*L3JetPar);

  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);
//----------------------------------------------------------------------------
//   Loop over all events
//----------------------------------------------------------------------------
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (entry % 500 == 0)
      std::cout << "--------------------" << entry << std::endl;

    for (int i = 0; i < branchRawJet->GetEntries(); ++i)
    {
      // Take first jet
      Jet *jet = (Jet*) branchRawJet->At(i);
      
      // Plot jet transverse momentum
      histRawJetPT->Fill(jet->PT);
      histRawJetPTCorr->Fill(jet->PT * correction(*jet, JetCorrector, branchRho));
    }
    
    for (int i = 0; i < branchRawJetNoPU->GetEntries(); ++i)
    {
      // Take first jet
      Jet *jet = (Jet*) branchRawJetNoPU->At(i);
      
      // Plot jet transverse momentum
      histRawJetNoPUPT->Fill(jet->PT);
    }

    for (int i = 0; i < branchGenJet->GetEntries(); ++i)
    {
      // Take first jet
      Jet *jet = (Jet*) branchGenJet->At(i);
      
      // Plot jet transverse momentum
      histGenJetPT->Fill(jet->PT);
      
    }

    for (int i = 0; i < branchJet->GetEntries(); ++i)
    {
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(i);
      
      // Plot jet transverse momentum
      histJetPT->Fill(jet->PT);
      
    }
    for (int i = 0; i < branchPileUpJetIDJet->GetEntries(); ++i)
    {
      // Take first jet
      Jet *jet = (Jet*) branchPileUpJetIDJet->At(i);
      
      // Plot jet transverse momentum
      histPUIDJetPT->Fill(jet->PT);
      
    }

    for (int i = 0; i < branchPuppiJet->GetEntries(); ++i)
    {
      // Take first jet
      Jet *jet = (Jet*) branchPuppiJet->At(i);
      
      // Plot jet transverse momentum
      histPUPPIJetPT->Fill(jet->PT);
    }

  } // End of looping event

  // Saving resulting histograms
  histRawJetPT->Write();
  histRawJetPTCorr->Write();
  histRawJetNoPUPT->Write();
  histGenJetPT->Write();
  histJetPT->Write();
  histPUIDJetPT->Write();
  histPUPPIJetPT->Write();

  outputfile.Close();
  return EXIT_SUCCESS;
}				// ----------  end of function main  ----------



double correction(Jet &iJet, FactorizedJetCorrector *iJetCorr, const TClonesArray *branchRho)
{ 
  iJetCorr->setJetPt (iJet.PT);
  iJetCorr->setJetEta(iJet.Eta);
  iJetCorr->setJetPhi(iJet.Phi);
  iJetCorr->setJetE  (iJet.P4().E());
  iJetCorr->setJetA  (iJet.AreaP4().Pt());

  double rho = 0.0;
  int avgcount = 0;
  for (int i = 0; i < branchRho->GetEntries(); ++i)
  {
    Rho *iRho = (Rho*) branchRho->At(i);
    rho += iRho->Rho;
    avgcount++;
  }
  rho = rho / avgcount;

  assert(rho != -999.);
  iJetCorr->setRho(rho);

  return iJetCorr->getCorrection();
}

