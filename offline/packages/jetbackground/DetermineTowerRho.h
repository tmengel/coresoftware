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
class JetContainer;
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
  // The background jets are clustered with fastjet and selected with
  // get_jet_selector() (eta, min pT). The surviving jets are then converted back
  // to sPHENIX Jets: the four-vector of each is rebuilt from the ORIGINAL
  // constituents, so the negative-energy towers -- which are flipped to +1 MeV for
  // the clustering only -- are reincorporated. The omit_nhardest hardest of these
  // sPHENIX jets are then dropped, and rho and sigma are estimated from the rest,
  // never from the fastjet jets.
  //
  // If jet_node is non-empty, every selected jet -- the omitted hardest ones
  // included -- is additionally written to the node tree as a JetContainer under
  // that name, so they can be analysed downstream. Each jet carries
  //   Jet::PROPERTY::prop_SeedItr : 1 if it was omitted as one of the hardest
  //                                 (a seed), 0 if it entered rho
  //   Jet::PROPERTY::prop_SeedD   : the signed scalar sum of its constituent pT
  //                                 (see set_use_signed_sum), always filled
  // The container also carries the resulting rho via
  // JetContainer::get_rho_median(). For the AREA method each jet carries its
  // fastjet area as Jet::PROPERTY::prop_area; jets with no area are stored but do
  // not enter the median, and pure-ghost jets (no real constituents) are stored
  // with zero momentum.
  void add_method(TowerRho::Method rho_method, std::string output_node = "",
                  std::string jet_node = "");

  // Choose which four-vector is stored for the saved background jets. This only
  // affects what is written to the container: rho is always estimated from the
  // sPHENIX jets built from the original constituents.
  //   false (default) : keep the sPHENIX jet as built from the ORIGINAL constituent
  //                     momenta. This is what rho is estimated from, and it matches
  //                     what FastJetAlgoSub writes, so the saved jets are directly
  //                     comparable to seed containers built by JetReco.
  //   true            : overwrite it with the raw fastjet axis, i.e. built from
  //                     constituents after negative-energy towers were flipped to
  //                     +1 MeV for clustering. Note this is the axis the jet eta
  //                     acceptance cut is applied to, so it is the one that
  //                     determines which jets entered rho -- useful for studying
  //                     that acceptance. The saved jets then no longer reproduce
  //                     rho by themselves.
  // The two differ appreciably only for soft jets; see get_jet_selector().
  void set_save_fastjet_axis(const bool b) { m_save_fastjet_axis = b; }

  // Choose the per-jet quantity the median is taken over.
  //   false (default) : |sum of the constituent pT vectors|, i.e. the jet pT. This is
  //                     the PPG14 estimator. Because it is a magnitude it is positive
  //                     definite: a jet whose towers sum to a net-negative energy
  //                     enters the median folded up as a positive value, which biases
  //                     rho high wherever such jets are common (peripheral events).
  //   true            : the signed scalar sum of the constituent pT, i.e.
  //                     sum_i sign(E_i) |pT_i|. Negative-energy towers then cancel
  //                     positive ones instead of being folded up, and a net-negative
  //                     jet contributes a negative value.
  // The two agree closely where all towers are positive, since a jet is narrow in phi
  // and the constituents are nearly collinear; they part company once negative towers
  // are a sizeable fraction. The quantity chosen here is also the one the
  // omit-n-hardest ranking uses, so that the jets removed are the hardest by the same
  // measure the median is taken over.
  void set_use_signed_sum(const bool b) { m_use_signed_sum = b; }
  bool get_use_signed_sum() const { return m_use_signed_sum; }
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

  // set the number of hardest jets to omit
  // This is applied after the fastjet jets have been converted to sPHENIX Jets,
  // to the jets that passed the eta (and min pT) cuts, and the jets are ranked by
  // the pT of the sPHENIX jet, i.e. with the negative-energy towers included.
  // Note that this is |pT|: a jet with a large net-negative energy is as "hard"
  // as a positive one of the same pT.
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

  // working container for the methods that were not given a jet_node: the
  // sPHENIX jets are built for every method so rho is always estimated from them
  JetContainer *m_scratch_jets{nullptr};

  bool m_save_fastjet_axis{false};
  bool m_use_signed_sum{false};

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

  // Convert the selected fastjet jets back to sPHENIX Jets, replacing the contents
  // of `jets`. Each four-vector is the sum of the ORIGINAL constituent momenta, so
  // negative-energy towers are reincorporated, and the constituents are recorded
  // as the jet's components. with_area must be true only when the jets came from a
  // ClusterSequenceArea, since PseudoJet::area() and is_pure_ghost() require area
  // information to be present; the area is then stored as Jet::prop_area.
  // signed_pt is filled, one entry per stored jet, with the signed scalar sum of the
  // constituent pT (see set_use_signed_sum). It cannot be recovered from the stored
  // four-vector, so it has to be accumulated here, while the constituents are at hand.
  void ConvertJets(JetContainer *jets, std::vector<float> &signed_pt,
                   const std::vector<fastjet::PseudoJet> &fastjets,
                   const std::vector<Jet *> &particles, bool with_area) const;

  // The per-jet quantity the median is taken over, and that the omit-n-hardest
  // ranking orders by: the jet pT, or the signed scalar sum (set_use_signed_sum).
  float JetPt(JetContainer *jets, const std::vector<float> &signed_pt,
              unsigned int ijet) const;

  // Drop the m_omit_nhardest hardest sPHENIX jets in `jets` (ranked by JetPt(),
  // ties keeping their original order). Returns the indices into `jets` of the
  // jets that remain, in their original order.
  std::vector<unsigned int> SelectJets(JetContainer *jets,
                                       const std::vector<float> &signed_pt) const;

  // Estimate rho and sigma from the sPHENIX jets of `jets` listed in `keep` (as
  // returned by SelectJets)
  void CalcRho(JetContainer *jets, const std::vector<float> &signed_pt,
               const std::vector<unsigned int> &keep,
               TowerRho::Method rho_method, float &rho, float &sigma) const;

  static float CalcPercentile(const std::vector<float> &sorted_vec,
                              const float percentile, const float nempty);

  static void CalcMedianStd(const std::vector<float> &vec,
                            float n_empty_jets, float &median, float &std_dev);

  // the fastjet-side selection: jet eta acceptance and, if set, the min pT. The
  // omit-n-hardest cut is not part of it, it is applied by SelectJets().
  fastjet::Selector get_jet_selector() const;
};

#endif  // JETBACKGROUND_DETERMINETOWERHO_H