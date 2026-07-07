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

  float get_UE_estimate( int layer, int eta , float layer_sumet );
  

 private:

  int CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  std::vector<std::vector<float> > _UE;

  std::string _backgroundName{"TestTowerBackground"};
  
  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};
  std::string EMTowerName;
  std::string IHTowerName;
  std::string OHTowerName;

  float _layer_sumet[3];
  float _seedless_UE_parameters_LAYER_ETA_PARAM[3][24][3];

};

#endif
