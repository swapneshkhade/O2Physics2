// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFCorrelatorD0Hadrons.cxx
 //Service-Oriented Architecture (SOA)
/// \author Shyam Kumar <shyam.kumar@cern.ch>
/// \author Swapnesh Khade <swapnesh.santosh.khade@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
//#include "Common/Core/TrackSelection.h"
//#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_correlation_d0hadron;
using namespace o2::constants::math;

// Additions for D0
using namespace o2::aod::hf_cand_prong2;
using namespace o2::analysis::hf_cuts_d0_topik;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies

double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PI / 2.);
}

/// definition of variables for D0 hadron pairs vs eta acceptance studies (hD0HadronVsEtaCut, in data-like, MC-reco and MC-kine tasks)
const double maxEtaCut = 5.;
const double ptThresholdForMaxEtaCut = 10.;
const double incrementEtaCut = 0.1;
const double incrementPtThreshold = 0.5;
const double epsilon = 1E-5;
// npTBinsMassAndEfficiency = 12 for D+ and 25 for D0
const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_d0_topik::npTBins;
const double efficiencyDmesonDefault[npTBinsMassAndEfficiency] = {};
auto efficiencyDmeson_v = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsMassAndEfficiency};

// histogram binning definition 
const int massAxisBins = 350;
const double massAxisMin = 1.7;
const double massAxisMax = 2.05;
const int phiAxisBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = 2. * o2::constants::math::PI;
const int yAxisBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int ptDAxisBins = 100;  //changed by swapnesh : 180 to 100
const double ptDAxisMin = 0.;
const double ptDAxisMax = 10.; //changed by swapnesh : 36 to 10

using MCParticlesPlus2Prong = soa::Join<aod::McParticles, aod::HfCandProng2MCGen>;

/// D0-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorD0HadronsOld {
  Produces<aod::D0HadronPair> entryD0HadronPair;
  Produces<aod::D0HadronRecoInfo> entryD0HadronRecoInfo;

  HistogramRegistry registry{
    "registry",
    //NOTE: use hMassD0 for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {
     {"hPtCand", "D0,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "D0,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "D0,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "D0,Hadron candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "D0,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "D0,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "D0,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPtCandMCRec", "D0,Hadron candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0MCRec", "D0,Hadron candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1MCRec", "D0,Hadron candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMCRec", "D0,Hadron candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEtaMCRec", "D0,Hadron candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCRec", "D0,Hadron candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCRec", "D0,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "D0,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaMCGen", "D0,Hadron particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCGen", "D0,Hadron particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCGen", "D0,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hcountD0HadronPerEvent", "D0,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hD0DaughtersEtaCut", "D0 daughters vs #eta cut on D daughters;#eta_{max};entries", {HistType::kTH2F, {{(int)(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {(int)(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hD0HadronVsEtaCut", "D0,Hadron pairs vs #eta cut;#eta_{max};entries", {HistType::kTH2F, {{(int)(maxEtaCut / incrementEtaCut), 0., maxEtaCut}, {(int)(ptThresholdForMaxEtaCut / incrementPtThreshold), 0., ptThresholdForMaxEtaCut}}}}
   }
 };

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<double> cutPtCandMin{"cutPtCandMin", -1., "min. cand. pT"};
  Configurable<double> cutPtCandMax{"cutPtCandMax", -1., "max. cand. pT"};
  Configurable<std::vector<double>> bins{"ptBinsForMassAndEfficiency", std::vector<double>{o2::analysis::hf_cuts_d0_topik::pTBins_v}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for D0 meson"};
  Configurable<int> flagApplyEfficiency{"efficiencyFlagD", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hMassD0_2D", "D0 candidates;inv. mass (K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0Data", "D0 candidates;inv. mass (K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD0MCRec", "D0 candidates;inv. mass (K^{-}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD0MCRecSig", "D0 signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecBkg", "D0 background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcountD0triggersMCGen", "D0 trigger particles - MC gen;;N of trigger D0", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});  
  }

   // Select D0
  Filter filterSelectCandidates = (aod::hf_selcandidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= selectionFlagD0bar);     

  void processData(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>> const& candidates)
  {
    int nTracks = 0; 
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < -4.0 || track.eta() > 4.0) {
          continue;
        }
        if (abs(track.dcaXY()) > 0.0025 || abs(track.dcaZ()) > 0.0025) { // change for D0
          continue;
        }
        nTracks++;
      }
    }

    //printf("cuts on Y max %f \n",(double)cutYCandMax);
    // printf("cuts on pT max %f \n",(double)cutPtCandMax);
    //  printf("cuts on pT min %f \n",(double)cutPtCandMin);
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks); 
    if (nTracks < multMin || nTracks > multMax) {    // Not clear purpose of this condition
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    for (auto& candidate1 : candidates) {
      if (cutYCandMax >= 0. && std::abs(YD0(candidate1)) > cutYCandMax) {
        continue;
      }
      /*
      if (cutPtCandMin >= 0. && candidate1.pt() < cutPtCandMin) {
        continue;
      }
      if (candidate1.pt() > cutPtCandMax) {
        continue;
      }
      */
      //check decay channel flag for candidate1  // Possible error 2 : DecayType::D0ToPiK not found
      if (!(candidate1.hfflag() & 1 << DecayType::D0ToPiK)) { // check Left Shift and Right Shift Operators C++
        continue;
      }
      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
        //LOGF(info, "efficiencyDmeson: %f",efficiencyWeight);
      }
      //fill invariant mass plots and generic info from all D0 candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMassD0_2D"), InvMassD0(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassD0Data"), InvMassD0(candidate1), efficiencyWeight);
      }
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      //registry.fill(HIST("hPtProng2"), candidate1.ptProng2());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), YD0(candidate1));
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelD0());

      //D0-Hadron correlation dedicated section
      //if the candidate is a D0, search for Hadrons and evaluate correlations
      if (candidate1.isSelD0() < selectionFlagD0) {
        continue;
      }
       for (const auto& track : tracks) {
        //Removing D0 daughters by checking track indices
        if ((candidate1.index0Id() == track.mRowIndex) || (candidate1.index1Id() == track.mRowIndex)) {
          continue;
        }
        if (abs(track.dcaXY()) >= 1.0 || abs(track.dcaZ()) >= 1.0) continue; // Remove secondary tracks
        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                         track.eta() - candidate1.eta(),
                         candidate1.pt(),
                         track.pt());
        entryD0HadronRecoInfo(InvMassD0(candidate1),0);
        double etaCut = 0.;
        double ptCut = 0.;
        do { //fill pairs vs etaCut plot
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do { //fill pairs vs etaCut plot
            if (std::abs(candidate1.eta()) < etaCut && std::abs(track.eta()) < etaCut && candidate1.pt() > ptCut && track.pt() > ptCut) {
              registry.fill(HIST("hD0HadronVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            ptCut += incrementPtThreshold;
          } while (ptCut < ptThresholdForMaxEtaCut - epsilon);
        } while (etaCut < maxEtaCut - epsilon);
        //note: candidates selected as both D0 and D0bar are used, and considered in both situation (but not auto-correlated): reflections could play a relevant role.
        //another, more restrictive, option, could be to consider only candidates selected with a single option (D0 xor D0bar)

      } // Hadron Tracks loop

    } //end outer loop
  }

PROCESS_SWITCH(HfCorrelatorD0HadronsOld, processData, "Process data", false);


/// D0-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
void processMcRec(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::HfCandProng2MCRec>> const& candidates)
  {
    
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < -4.0 || track.eta() > 4.0) {
          continue;
        }
        if (abs(track.dcaXY()) > 0.0025 || abs(track.dcaZ()) > 0.0025) {
          continue;
        }
        nTracks++;
      }
    }

    //LOGP(info, "nTracks = {}", nTracks); // Added by swapnesh to scheck error
    //printf("cuts on Y max %f \n",(double)cutYCandMax);
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    //MC reco level
     bool flagD0Signal = false;
    for (auto& candidate1 : candidates) {
      //check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YD0(candidate1)) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && candidate1.pt() < cutPtCandMin) {
        continue;
      }
      if (candidate1.pt() >= cutPtCandMax) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }

      if (std::abs(candidate1.flagMCMatchRec()) == 1 << DecayType::D0ToPiK) {
        //fill per-candidate distributions from D0 true candidates
        registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
        //registry.fill(HIST("hPtProng2MCRec"), candidate1.ptProng2());
        registry.fill(HIST("hEtaMCRec"), candidate1.eta());
        registry.fill(HIST("hPhiMCRec"), candidate1.phi());
        registry.fill(HIST("hYMCRec"), YD0(candidate1));
        registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelD0());
      }
      //fill invariant mass plots from D0 signal and background candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {                  //only reco as D0
        registry.fill(HIST("hMassD0MCRec"), InvMassD0(candidate1), efficiencyWeight);
        if (candidate1.flagMCMatchRec() == 1 << DecayType::D0ToPiK) { //also matched as D0
          registry.fill(HIST("hMassD0MCRecSig"), InvMassD0(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0MCRecBkg"), InvMassD0(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }

      //D0-Hadron correlation dedicated section
      //if the candidate is selected as D0, search for Hadron and evaluate correlations
      if (candidate1.isSelD0() < selectionFlagD0) { //discard candidates not selected as D0 in outer loop
        continue;
      }
      flagD0Signal = candidate1.flagMCMatchRec() == 1 << DecayType::D0ToPiK;   //flagD0Signal 'true' if candidate1 matched to D0 (particle)
      for (const auto& track : tracks) {

        //Removing D0 daughters by checking track indices
        if ((candidate1.index0Id()== track.mRowIndex) || (candidate1.index1Id() == track.mRowIndex)) {
          continue;
        }
        if (abs(track.dcaXY()) >= 1.0 || abs(track.dcaZ()) >= 1.0) continue; // Remove secondary tracks
        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                         track.eta() - candidate1.eta(),
                         candidate1.pt(),
                         track.pt());
        entryD0HadronRecoInfo(InvMassD0(candidate1), (int)flagD0Signal);
        double etaCut = 0.;
        double ptCut = 0.;
        do { //fill pairs vs etaCut plot
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do { //fill pairs vs etaCut plot
            if (std::abs(candidate1.eta()) < etaCut && std::abs(track.eta()) < etaCut && candidate1.pt() > ptCut && track.pt() > ptCut) {
              registry.fill(HIST("hD0HadronVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            ptCut += incrementPtThreshold;
          } while (ptCut < ptThresholdForMaxEtaCut - epsilon);
        } while (etaCut < maxEtaCut - epsilon);
      } // end inner loop (Tracks)

    } //end outer loop 
  }
  

 PROCESS_SWITCH(HfCorrelatorD0HadronsOld, processMcRec, "Process MC Reco mode", true)


/// D0-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::McCollision const& mccollision, MCParticlesPlus2Prong const& particlesMC)  // swapnesh : Find the source of mcparticlesplus2prongs
  {
    int counterD0Hadron = 0;
    registry.fill(HIST("hMCEvtCount"), 0);
    //MC gen level
    for (auto& particle1 : particlesMC) {
      //check if the particle is D0  (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != pdg::Code::kD0) {
        continue;
      }
      double yD = RecoDecay::y(array{particle1.px(), particle1.py(), particle1.pz()}, RecoDecay::getMassPDG(particle1.pdgCode()));
      if (cutYCandMax >= 0. && std::abs(yD) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && particle1.pt() < cutPtCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMCGen"), particle1.eta());
      registry.fill(HIST("hPhiMCGen"), particle1.phi());
      registry.fill(HIST("hYMCGen"), yD);
      counterD0Hadron++;
      //D0 Hadron correlation dedicated section
      //if it's a D0 particle, search for Hadron and evaluate correlations
      if (particle1.pdgCode() != pdg::Code::kD0) { //just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hcountD0triggersMCGen"), 0, particle1.pt()); //to count trigger D0 for normalisation)
      for (auto& particle2 : particlesMC) {
        if ((std::abs(particle2.pdgCode()) != 11) && (std::abs(particle2.pdgCode()) != 13) && (std::abs(particle2.pdgCode()) != 211) && (std::abs(particle2.pdgCode()) != 321) && (std::abs(particle2.pdgCode()) != 2212)) { //check that inner particle is D0
          continue;
        }
        entryD0HadronPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                         particle2.eta() - particle1.eta(),
                         particle1.pt(),
                         particle2.pt());
        double etaCut = 0.;
        double ptCut = 0.;

        //fill pairs vs etaCut plot
        bool rightDecayChannels = false;
        if ((std::abs(particle1.flagMCMatchGen()) == 1 << DecayType::D0ToPiK)) {
          rightDecayChannels = true;
        }
        do {
          ptCut = 0.;
          etaCut += incrementEtaCut;
          do {                                                                                                                                  //fill pairs vs etaCut plot
            if (std::abs(particle1.eta()) < etaCut && std::abs(particle2.eta()) < etaCut && particle1.pt() > ptCut && particle2.pt() > ptCut) { //fill with D0 Hadron in acceptance checks
              registry.fill(HIST("hD0HadronVsEtaCut"), etaCut - epsilon, ptCut + epsilon);
            }
            if (rightDecayChannels) { //fill with D0 daughter particles acceptance checks
              bool candidate1DauInAcc = true;
              for (auto& dau : particle1.daughters_as<MCParticlesPlus2Prong>()) {
                if (std::abs(dau.eta()) > etaCut) {
                  candidate1DauInAcc = false;
                  break;
                }
              }
              if (candidate1DauInAcc && particle1.pt() > ptCut && particle1.pt() > ptCut) {
              registry.fill(HIST("hD0DaughtersEtaCut"), etaCut - epsilon, ptCut + epsilon);
              }
              }
            ptCut += incrementPtThreshold;
          } while (ptCut < ptThresholdForMaxEtaCut - epsilon);
        } while (etaCut < maxEtaCut - epsilon);
      } //end inner loop
    }   //end outer loop
    registry.fill(HIST("hcountD0HadronPerEvent"), counterD0Hadron);
  }
   PROCESS_SWITCH(HfCorrelatorD0HadronsOld, processMcGen, "Process MC Gen mode", false);
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorD0HadronsOld>(cfgc)};
}