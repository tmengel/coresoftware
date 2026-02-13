// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef DIJETQA_H
#define DIJETQA_H

#include "JetQADefs.h"

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <utility>

class JetContainer;
class Fun4AllHistoManager;
class TriggerAnalyzer;
class PHCompositeNode;
class TH1;
class TH2;

class DijetQA : public SubsysReco
{
 public:
  DijetQA(const std::string &name = "DijetQA", const std::string &recojetname = "AntiKt_Tower_r04");

  ~DijetQA() override = default;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;
  void FindPairs(JetContainer *);
  int grab_zvrtx(PHCompositeNode *topNode); // grab the zvrtx from the global vertex node and fill the member variable, return 0 if successful, 1 if not

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Set the eta range for reco jets
  void setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }

  /// Set the leading pt range for the reco jets
  void setPtLeadRange(double low, double high)
  {
    m_ptLeadRange.first = low;
    m_ptLeadRange.second = high;
  }

  /// Set the sub-leading pt range for the reco jets
  void setPtSubRange(double low, double high)
  {
    m_ptSubRange.first = low;
    m_ptSubRange.second = high;
  }

  /// Specifies a trigger to require
  void setTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSJet1)
  {
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }

  //////////////////////////////////////////////////////////////
  //							    //
  //        X_j = (p_(T, 1))/(p_(T,2))			    //
  //        A_j = (p_(T, 1) -p_(T,2))/P_T			    //
  //							    //
  //////////////////////////////////////////////////////////////

 private:
  Fun4AllHistoManager *m_manager{nullptr};
  TriggerAnalyzer *m_analyzer{nullptr};

  TH1 *h_Ajj{nullptr};
  TH1 *h_xj{nullptr};
  TH1 *h_pt{nullptr};
  TH1 *h_dphi{nullptr};
  TH2 *h_Ajj_pt{nullptr};
  TH2 *h_xj_pt{nullptr};
  TH2 *h_dphi_pt{nullptr};
  TH2 *h_dphi_Ajj{nullptr};
  TH1 *h_Ajj_l{nullptr};
  TH1 *h_xj_l{nullptr};
  TH1 *h_pt_l{nullptr};
  TH1 *h_dphi_l{nullptr};
  TH2 *h_Ajj_pt_l{nullptr};
  TH2 *h_xj_pt_l{nullptr};
  TH2 *h_dphi_pt_l{nullptr};
  TH2 *h_dphi_Ajj_l{nullptr};

  int m_nJet{-1};
  int m_nJetPair{-1};
  uint32_t m_trgToSelect{false};

  float DeltaPhi{M_PI * 0.25};  // cut on the opening angle of phi for the identified jets
                                // using same value as the dijet analysis
  float m_zvtx{-1};
  float m_Ajj{-1};
  float m_xj{-1};
  float m_ptl{-1};
  float m_ptsl{-1};
  float m_phil{-1};
  float m_phisl{-1};
  float m_dphil{-1};
  float m_etal{-1};
  float m_etasl{-1};
  float m_deltaeta{-1};

  bool m_doTrgSelect = (JetQADefs::GL1::MBDNSJet1);  // maps an int to a boolean, needs static cast if {} is used

  std::string m_moduleName;
  std::string m_recoJetName;

  std::pair<float, float> m_etaRange{-0.7, 0.7};
  std::pair<float, float> m_ptLeadRange{1, 100};
  std::pair<float, float> m_ptSubRange{1, 100};


  // easy function to check if a jet is within the detector acceptance, given the jet eta, the jet radius, and the z vertex position
  static bool jet_in_calo_acceptance( const float eta, const float R, const float z )
  {
      const float CALO_ABS_Z[3] = {130.23, 170.299, 301.683};
      const float CALO_RADIUS[3] = {93.5, 127.503, 225.87};
      float min_eta = -999, MAX_eta = 999;
      for (int i=0; i<3; i++) 
      {
        float min_z = -1.0*CALO_ABS_Z[i];
        float max_z = CALO_ABS_Z[i];
        float z_l = min_z - zvrtx;
        float z_h = max_z - zvrtx;
        float r = CALO_RADIUS[i];
        float eta_min = asinh( z_l / r );
        float eta_max = asinh( z_h / r );
        if ( eta_min > min_eta ) { min_eta = eta_min; }
        if ( eta_max < MAX_eta ) { MAX_eta = eta_max; }
      }
      min_eta += R;
      MAX_eta -= R;
      return (eta >= min_eta && eta <= MAX_eta);
  }

  

  


};

#endif  // DIJETQA_H
