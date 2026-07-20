#ifndef JETBACKGROUND_AverageESub_H
#define JETBACKGROUND_AverageESub_H

//===========================================================
/// \file AverageESub.h
/// \brief alternative UE background calculator w/o seeds
/// \author Dennis V. Perepelitsa, Tanner Mengel
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>
#include <array>
#include <vector>

class PHCompositeNode;


class AverageESub : public SubsysReco
{
 public:

  AverageESub( const std::string &name = "AverageESub" ) : SubsysReco(name) {}
  ~AverageESub() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_ue_directpath( const std::string &name )  { m_ue_directpath = name; }
  void set_tower_node_prefix( const std::string &prefix ) { m_towernode_prefix = prefix; }
  void set_tower_sub_prefix( const std::string &prefix ) { m_tower_sub_prefix = prefix; }
 
  void set_towerbackground_name( const std::string &name ) { m_background_node = name; }

  void set_overlay_node( const std::string &name ) { m_overlay_node = name; } // overrides z and mbd

 private:

  int CreateNode(PHCompositeNode *topNode);

  std::string m_ue_directpath;
  std::string m_towernode_prefix{"TOWERINFO_CALIB"};
  std::string m_tower_sub_prefix{"AVG_SUB_TOWERINFO_CALIB"};
  std::string m_background_node{"TowerBackground_AvgE"};

  std::string m_overlay_node;

  double m_mbdQ = 0;
  double m_zvertex = 0;

  std::vector<float>* m_avg_tower_cemc_e = nullptr;
  std::vector<float>* m_avg_tower_hcalin_e = nullptr;
  std::vector<float>* m_avg_tower_hcalout_e = nullptr;

  static constexpr int k_n_layers       = 3;
  static constexpr int k_n_eta_bins     = 24;
  static constexpr int k_n_phi_bins     = 64;
  static constexpr int k_n_zvertex_bins = 11;
  static constexpr int k_n_mbd_bins     = 42;

  static constexpr std::array<double, k_n_zvertex_bins + 1> k_zvertex_bins = {
    -60, -50, -40, -30, -20, -10,
     10,  20,  30,  40,  50,  60
  };

  static constexpr std::array<double, k_n_mbd_bins + 1> k_mbd_bins = {
    -100, -50, 0,
     50, 100, 150,
    200, 250, 300,
    350, 400, 450,
    500, 550, 600,
    650, 700, 750,
    800, 850, 900,
    950, 1000, 1050, 1100,
    1150, 1200, 1250, 1300,
    1350, 1400, 1450, 1500,
    1550, 1600, 1650, 1700,
    1750, 1800, 1850, 1900,
    1950, 2000
  };

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

  static int get_mbd_bin(const double mbd)
  {
    for (int i = 0; i < k_n_mbd_bins; ++i)
    {
      if (mbd >= k_mbd_bins[i] &&
          mbd <  k_mbd_bins[i + 1])
      {
        return i;
      }
    }
    return -1;
  }

  static size_t get_index(const int ieta,
                          const int iphi,
                          const int izbin,
                          const int imbd)
  {
    return izbin * (k_n_mbd_bins * k_n_phi_bins * k_n_eta_bins)
         + imbd  * (k_n_phi_bins * k_n_eta_bins)
         + iphi  * k_n_eta_bins
         + ieta;
  }

  static void decode_index(int& ieta,
                           int& iphi,
                           int& izbin,
                           int& imbd,
                           size_t index)
  {
    izbin = index / (k_n_mbd_bins * k_n_phi_bins * k_n_eta_bins);
    index %= (k_n_mbd_bins * k_n_phi_bins * k_n_eta_bins);

    imbd = index / (k_n_phi_bins * k_n_eta_bins);
    index %= (k_n_phi_bins * k_n_eta_bins);

    iphi = index / k_n_eta_bins;
    ieta = index % k_n_eta_bins;
  }

};

#endif
