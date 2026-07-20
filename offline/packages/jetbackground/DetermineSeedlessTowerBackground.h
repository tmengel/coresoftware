#ifndef JETBACKGROUND_DETERMINESEEDLESSTOWERBACKGROUND_H
#define JETBACKGROUND_DETERMINESEEDLESSTOWERBACKGROUND_H

//===========================================================
/// \file DetermineSeedlessTowerBackground.h
/// \brief alternative UE background calculator w/o seeds
/// \author Dennis V. Perepelitsa, Tanner Mengel
//===========================================================

#include <fun4all/SubsysReco.h>

// system includes
#include <jetbase/Jet.h>
#include <string>
#include <vector>
#include <array>

// forward declarations
class PHCompositeNode;

/// \class DetermineSeedlessTowerBackground
///
/// \brief alternative UE background calculator w/o seeds
///
/// This module constructs dE/deta vs. eta given the overall SumET in
/// different calorimeter layers
///
class DetermineSeedlessTowerBackground : public SubsysReco
{
 public:
  DetermineSeedlessTowerBackground(const std::string &name = "DetermineSeedlessTowerBackground");
  ~DetermineSeedlessTowerBackground() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void SetBackgroundOutputName(const std::string &name) { _backgroundName = name; }

  void set_towerNodePrefix(const std::string &prefix)
  {
    m_towerNodePrefix = prefix;
    return;
  }

  void set_ue_directpath( const std::string &name )  { m_ue_directpath = name; }

  float get_UE_estimate( int iz, int layer, int eta, float mbdQ_sum );

  void set_overlay_node( const std::string &name ) { m_overlay_node = name; } // overrides z and mbd
  

 private:

  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  std::vector<std::vector<float>> _UE;

  std::string m_ue_directpath;

  std::string m_overlay_node;

  static constexpr int k_n_zvertex_bins = 11;
  static constexpr std::array< double, k_n_zvertex_bins + 1> k_zvertex_bins = {
    -60, -50, -40, -30, -20, -10,
    10,  20,  30,  40,  50,  60
  };

  static constexpr int k_n_layers   = 3;
  static constexpr int k_n_eta_bins = 24;
  static constexpr int k_n_params   = 3; // p0, p1, p2

  static int get_zbin(const double zvertex)
  {
    for (int i = 0; i < k_n_zvertex_bins; ++i)
    {
      if (zvertex >= k_zvertex_bins[i] &&
          zvertex <  k_zvertex_bins[i + 1])
      {
        return i;
      }
    }
    return -1;
  }

  double m_parameters[k_n_zvertex_bins]
                     [k_n_layers]
                     [k_n_eta_bins]
                     [k_n_params];


  std::string _backgroundName{"TestTowerBackground"};
  
  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};
  std::string EMTowerName;
  std::string IHTowerName;
  std::string OHTowerName;

  // float _layer_sumet[3];
  // float _seedless_UE_parameters_LAYER_ETA_PARAM[3][24][3];

};

#endif
