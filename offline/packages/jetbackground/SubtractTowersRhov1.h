#ifndef JETBACKGROUND_SUBTRACTTOWERSRHOV1_H
#define JETBACKGROUND_SUBTRACTTOWERSRHOV1_H

#include <fun4all/SubsysReco.h>

#include "TowerRho.h"
#include <globalvertex/GlobalVertex.h>

#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;
class CDBTTree;

/// \class SubtractTowersRhov1
///
/// \brief creates new UE-subtracted towers, with an optional eta-shape
///        calibration on top of the flat per-event rho
///
/// Same flat-rho subtraction as SubtractTowersRho ( E' = E - rho*cosh(eta),
/// optionally area- or multiplicity-weighted ) but with an optional
/// multiplicative eta-shape correction w(layer,ieta;zvertex,mbdQ):
///
///     E' = E - rho * w(layer,ieta;zvertex,mbdQ) * cosh(eta)
///
/// w corrects for the facts that (i) the calorimeter response and (ii) the
/// underlying event dN/deta are not flat in eta, and that this shape shifts
/// with the collision z-vertex and depends on centrality (here proxied by
/// the raw MBD charge sum, to avoid a dependency on a separate centrality
/// calibration). w is built by rho_calib/make_rho_calib.C as a simple
/// binned lookup table and stored in CDBTTree format; see that macro for
/// details. If no calibration is configured, w === 1 and the behavior is
/// identical to SubtractTowersRho.
///
/// The calibration table can be loaded two ways -- whichever is set last
/// wins if both are called:
///   - set_etaCalib_directPath(path): read directly from a ROOT file on disk
///   - set_etaCalib_cdbTag(tag): resolve a path via CDBInterface, for use
///     once the payload has been committed as an actual CDB calibration
/// Both paths load through the exact same CDBTTree machinery, so switching
/// from one to the other is a one-line change.
class SubtractTowersRhov1 : public SubsysReco
{
 public:
  SubtractTowersRhov1(const std::string &name = "SubtractTowersRhov1") : SubsysReco(name) {}
  ~SubtractTowersRhov1() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_flowMod( const bool b ) { m_doFlowMod = b; }
  void set_rhoNode( const std::string & node ) { m_rhoNode = node; }
  void add_targetTowerNode( const std::string & node ) { m_targetTowerNodes.push_back(node); }
  void set_subSuffix( const std::string & node ) { m_subSuffix = node; }
  void set_globalVertexType( GlobalVertex::VTXTYPE type ) { m_vertex_type = type; }

  /// Point directly at a rho eta-shape calibration file on disk (CDBTTree
  /// format, as written by rho_calib/make_rho_calib.C).
  void set_etaCalib_directPath( const std::string & path ) { m_calib_direct_path = path; m_use_etaCalib = true; }

  /// Resolve the calibration through CDBInterface using a CDB global tag
  /// name, e.g. once the payload has been committed to the calibration
  /// database.
  void set_etaCalib_cdbTag( const std::string & tag ) { m_calib_cdb_tag = tag; m_use_etaCalib = true; }

  /// Node name for the MBD raw charge sum used to index the eta-shape
  /// calibration; only read from the node tree when a calibration is
  /// configured.
  void set_mbdNode( const std::string & node ) { m_mbd_node = node; }

 private:

  int CreateNode( PHCompositeNode *topNode );
  int grab_zvrtx( PHCompositeNode *topNode );
  int grab_mbdQ( PHCompositeNode *topNode );
  int LoadEtaCalib();
  float get_etaWeight( const int layer_index, const int ieta ) const;

  bool m_doFlowMod = false;
  std::string m_rhoNode = "";

  std::vector< std::string > m_targetTowerNodes {};

  std::string m_subSuffix = "";

  TowerRho::Method m_rho_method = TowerRho::Method::NONE;

  static std::string get_outputTowerNode( const std::string & input_node, const std::string & sub_suffix )
  {
    return  sub_suffix + "_" + input_node;
  }

  static int get_layer_index( const bool is_emcal, const bool is_ihcal, const bool is_ohcal )
  {
    if ( is_emcal ) { return 0; }
    if ( is_ihcal ) { return 1; }
    if ( is_ohcal ) { return 2; }
    return -1;
  }

  GlobalVertex::VTXTYPE m_vertex_type {GlobalVertex::UNDEFINED};
  float m_vtxz = 0.0;

  // -- optional eta-shape calibration --
  bool m_use_etaCalib = false;
  bool m_etaCalib_loaded = false;
  std::string m_calib_direct_path = "";
  std::string m_calib_cdb_tag = "";
  std::string m_mbd_node = "MbdOut";
  float m_mbdQ = 0.0;

  static constexpr int m_calib_n_eta_expected = 24;
  CDBTTree * m_calib_tree = nullptr;
  int m_calib_izbin = -1;
  int m_calib_imbd = -1;
  std::vector<float> m_calib_zvtx_edges {};
  std::vector<float> m_calib_mbdQ_edges {};
  static const std::vector<std::string> m_calib_layer_names;

  static int find_bin( const float val, const std::vector<float> & edges );
  static int encode_channel( const int ieta, const int izbin, const int imbd, const int n_mbd_bins )
  {
    return izbin * ( n_mbd_bins * m_calib_n_eta_expected ) + imbd * m_calib_n_eta_expected + ieta;
  }

};

#endif
