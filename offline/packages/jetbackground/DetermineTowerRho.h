/*!
 * \file DetermineTowerRho.h
 * \brief UE background rho calculator.
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 2.0.1 $
 * \date $Date: 02/01/2024. Revised 09/19/2024$
 */

#ifndef JETBACKGROUND_DETERMINETOWERHO_H
#define JETBACKGROUND_DETERMINETOWERHO_H

#include "TowerRho.h"

#include <fun4all/SubsysReco.h>

#include <iostream>
#include <string>
#include <vector>

class PHCompositeNode;
class JetInput;
class Jet;

namespace fastjet
{
  class PseudoJet;
  class Selector;
}  // namespace fastjet

class DetermineTowerRho : public SubsysReco
{
 public:
  DetermineTowerRho(const std::string &name = "DetermineTowerRho");
  ~DetermineTowerRho() override;

  // standard Fun4All methods
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  // add rho method (Area or Multiplicity)
  //
  // If jet_node is non-empty, the background jets that this method actually used
  // -- i.e. the ones surviving get_jet_selector() and contributing to the median
  // -- are additionally written to the node tree as a JetContainer under that
  // name, so they can be analysed downstream. The container also carries the
  // resulting rho via JetContainer::get_rho_median().
  void add_method(TowerRho::Method rho_method, std::string output_node = "",
                  std::string jet_node = "");

  // Choose which four-vector is stored for the saved background jets.
  //   false (default) : rebuilt by summing the ORIGINAL constituent momenta.
  //                     This matches what FastJetAlgoSub writes, so the saved
  //                     jets are directly comparable to seed containers built
  //                     by JetReco.
  //   true            : the raw fastjet axis, i.e. built from constituents after
  //                     negative-energy towers were flipped to +1 MeV for
  //                     clustering. Note this is the axis the jet eta acceptance
  //                     cut is applied to, so it is the one that determines which
  //                     jets entered rho -- useful for studying that acceptance.
  // The two differ appreciably only for soft jets; see get_jet_selector().
  void set_save_fastjet_axis(const bool b) { m_save_fastjet_axis = b; }
  bool get_save_fastjet_axis() const { return m_save_fastjet_axis; }

  // inputs for estimating background
  void add_input(JetInput *input) { m_inputs.push_back(input); }
  void add_tower_input(JetInput *input) { add_input(input); }  // for backwards compatibility

  // set the jet algorithm used to cluster background jets
  // default is KT
  void set_algo(const Jet::ALGO algo) { m_bkgd_jet_algo = algo; }
  Jet::ALGO get_algo() const { return m_bkgd_jet_algo; }

  // set the jet algorithm parameter for background jets
  // default is 0.4
  void set_par(const float val) { m_par = val; }
  float get_par() { return m_par; }

  // set the absolute eta range for tower acceptance
  // default is 1.1
  void set_abs_eta(const float val) { m_abs_input_eta_range = val; }
  float get_abs_eta() const { return m_abs_input_eta_range; }

  void set_tower_abs_eta(const float val) { set_abs_eta(val); }  // for backwards compatibility
  float get_tower_abs_eta() const { return get_abs_eta(); }      // for backwards compatibility

  // set the absolute eta range for jet acceptance
  // default is 1.1
  void set_jet_abs_eta(float abseta) { m_abs_jet_eta_range = abseta; }
  float get_jet_abs_eta() const { return m_abs_jet_eta_range; }

  // set the number of hardest towers to omit
  // default is 2
  void set_omit_nhardest(const unsigned int val) { m_omit_nhardest = val; }
  unsigned int get_omit_nhardest() const { return m_omit_nhardest; }

  // set the ghost area
  // default is 0.01
  void set_ghost_area(const float val) { m_ghost_area = val; }
  float get_ghost_area() const { return m_ghost_area; }

  // set the minimum pT for jets accepted in the background estimation
  // default is off (VOID_CUT)
  void set_jet_min_pT(const float val) { m_jet_min_pT = val; }
  float get_jet_min_pT() const { return m_jet_min_pT; }

  // print settings
  void print_settings(std::ostream &os = std::cout);

 private:
  // variables
  std::vector<JetInput *> m_inputs{};
  std::vector<std::string> m_output_nodes{};
  std::vector<std::string> m_jet_output_nodes{};  // empty string = do not save jets
  std::vector<TowerRho::Method> m_rho_methods{};

  bool m_save_fastjet_axis{false};

  Jet::ALGO m_bkgd_jet_algo{Jet::ALGO::KT};  // default is KT
  float m_par{0.4};                          // default is 0.4

  float m_abs_input_eta_range{1.1};  // default is 1.1
  unsigned int m_omit_nhardest{2};   // default is 2
  float m_ghost_area{0.01};          // default is 0.01

  const float VOID_CUT{-999.0};
  float m_jet_min_pT{-999.0};         // default is off
  float m_abs_jet_eta_range{-999.0};  // default is off

  // internal methods
  int CreateNodes(PHCompositeNode *topNode);

  // Write the background jets used for one rho method into its JetContainer node.
  // has_ghosts must be true only when the jets came from a ClusterSequenceArea,
  // since PseudoJet::is_pure_ghost() requires area information to be present.
  void FillJetContainer(PHCompositeNode *topNode, const std::string &node_name,
                        const std::vector<fastjet::PseudoJet> &fastjets,
                        const std::vector<float> &areas,
                        const std::vector<Jet *> &particles,
                        bool has_ghosts, float rho);

  static float CalcPercentile(const std::vector<float> &sorted_vec,
                              const float percentile, const float nempty);

  static void CalcMedianStd(const std::vector<float> &vec,
                            float n_empty_jets, float &median, float &std_dev);

  fastjet::Selector get_jet_selector() const;
};

#endif  // JETBACKGROUND_DETERMINETOWERHO_H