#include "SubtractTowersRhov1.h"

#include "TowerRho.h"
#include "TowerRhov1.h"

// sPHENIX includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv1.h>

#include <globalvertex/GlobalVertexv3.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <mbd/MbdOutV2.h>

#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TMath.h>
#include <TString.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <cassert>

const std::vector<std::string> SubtractTowersRhov1::m_calib_layer_names = {"w_cemc", "w_hcalin", "w_hcalout"};

SubtractTowersRhov1::~SubtractTowersRhov1()
{
  delete m_calib_tree;
}

int SubtractTowersRhov1::find_bin( const float val, const std::vector<float> & edges )
{
  if ( edges.size() < 2 || std::isnan(val) || val < edges.front() || val >= edges.back() )
  {
    return -1;
  }
  // edges.size() is small (O(10-20)) -- linear search is fine
  for ( size_t i = 0; i + 1 < edges.size(); ++i )
  {
    if ( val >= edges[i] && val < edges[i + 1] )
    {
      return static_cast<int>(i);
    }
  }
  return -1;
}

int SubtractTowersRhov1::LoadEtaCalib()
{
  if ( !m_use_etaCalib )
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::string calibpath = m_calib_direct_path;
  if ( calibpath.empty() && !m_calib_cdb_tag.empty() )
  {
    calibpath = CDBInterface::instance()->getUrl( m_calib_cdb_tag );
  }

  if ( calibpath.empty() )
  {
    std::cout << PHWHERE << "SubtractTowersRhov1::LoadEtaCalib - no calibration path resolved "
              << "(direct path and CDB tag both empty/unresolved). Falling back to flat rho "
              << "(w === 1 everywhere)." << std::endl;
    m_use_etaCalib = false;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_calib_tree = new CDBTTree( calibpath );
  m_calib_tree->LoadCalibrations();

  const int n_eta = m_calib_tree->GetSingleIntValue( "n_eta" );
  const int n_z = m_calib_tree->GetSingleIntValue( "n_zvtx_bins" );
  const int n_mbd = m_calib_tree->GetSingleIntValue( "n_mbdQ_bins" );
  if ( n_eta != m_calib_n_eta_expected || n_z <= 0 || n_mbd <= 0 )
  {
    std::cout << PHWHERE << "SubtractTowersRhov1::LoadEtaCalib - calibration file " << calibpath
              << " looks invalid (n_eta=" << n_eta << ", n_zvtx_bins=" << n_z
              << ", n_mbdQ_bins=" << n_mbd << "). Disabling eta-shape calibration." << std::endl;
    delete m_calib_tree;
    m_calib_tree = nullptr;
    m_use_etaCalib = false;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_calib_zvtx_edges.resize( n_z + 1 );
  for ( int i = 0; i <= n_z; ++i )
  {
    m_calib_zvtx_edges[i] = m_calib_tree->GetSingleFloatValue( Form( "zvtx_edge_%d", i ) );
  }
  m_calib_mbdQ_edges.resize( n_mbd + 1 );
  for ( int i = 0; i <= n_mbd; ++i )
  {
    m_calib_mbdQ_edges[i] = m_calib_tree->GetSingleFloatValue( Form( "mbdQ_edge_%d", i ) );
  }

  m_etaCalib_loaded = true;
  if ( Verbosity() > 0 )
  {
    std::cout << "SubtractTowersRhov1::LoadEtaCalib - loaded eta-shape calibration from " << calibpath
              << " (n_zvtx_bins=" << n_z << ", n_mbdQ_bins=" << n_mbd << ")" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

float SubtractTowersRhov1::get_etaWeight( const int layer_index, const int ieta ) const
{
  if ( !m_use_etaCalib || !m_etaCalib_loaded || !m_calib_tree )
  {
    return 1.0f;
  }
  if ( layer_index < 0 || layer_index >= static_cast<int>(m_calib_layer_names.size()) )
  {
    return 1.0f;
  }
  if ( m_calib_izbin < 0 || m_calib_imbd < 0 )
  {
    return 1.0f;  // event's (zvertex,mbdQ) falls outside the calibrated range
  }
  if ( ieta < 0 || ieta >= m_calib_n_eta_expected )
  {
    return 1.0f;
  }

  const int n_mbd_bins = static_cast<int>(m_calib_mbdQ_edges.size()) - 1;
  const int channel = encode_channel( ieta, m_calib_izbin, m_calib_imbd, n_mbd_bins );
  const float w = m_calib_tree->GetFloatValue( channel, m_calib_layer_names.at( layer_index ) );
  if ( !(w > 0) || std::isnan(w) )
  {
    return 1.0f;
  }
  return w;
}

int SubtractTowersRhov1::InitRun(PHCompositeNode *topNode)
{
  auto res = LoadEtaCalib();
  if ( res != Fun4AllReturnCodes::EVENT_OK )
  {
    return res;
  }
  return CreateNode(topNode);
}

int SubtractTowersRhov1::grab_zvrtx( PHCompositeNode *topNode )
{
  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowersRhov1::grab_zvrtx - starting grab_zvrtx with m_vertex_type = " << m_vertex_type << std::endl;
  }

  m_vtxz = 0;  // default to 0
  auto vertexmap = findNode::getClass< GlobalVertexMap >( topNode, "GlobalVertexMap" );
  if ( !vertexmap )
  {
    std::cout << "SubtractTowersRhov1::grab_zvrtx - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if ( vertexmap->empty() )
  {
    if (Verbosity() > 0 )
    {
      std::cout << "SubtractTowersRhov1::grab_zvrtx - empty vertex map, continuing as if zvtx = 0" << std::endl;
    }
    m_vtxz = 0;  // default to 0
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if ( m_vertex_type == GlobalVertex::UNDEFINED )
  {
    auto * vtx = vertexmap->begin()->second;
    if ( vtx )
    {
      m_vtxz = vtx->get_z();
    }
  }
  else
  {
    auto vertices = vertexmap -> get_gvtxs_with_type( { m_vertex_type } );
    if( !vertices.empty() && vertices.at(0) )
    {
      m_vtxz = vertices.at(0) -> get_z();
    }
  }

  if ( std::isnan(m_vtxz)  || std::abs(m_vtxz) > 1e3 )
  {
    static bool once = true;
    if (once)
    {
      once = false;
      std::cout << "SubtractTowersRhov1::grab_zvrtx - WARNING - vertex is " << m_vtxz << ". Continue with zvtx = 0 (further vertex warning will be suppressed)." << std::endl;
    }
    m_vtxz = 0;
  }

  if (Verbosity() > 1)
  {
    std::cout << "SubtractTowersRhov1::grab_zvrtx - finished grab_zvrtx with m_vtxz = " << m_vtxz << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowersRhov1::grab_mbdQ( PHCompositeNode *topNode )
{
  m_mbdQ = 0;
  auto * mbd_node = findNode::getClass<MbdOutV2>( topNode, m_mbd_node );
  if ( !mbd_node )
  {
    static bool once = true;
    if ( once )
    {
      once = false;
      std::cout << PHWHERE << "SubtractTowersRhov1::grab_mbdQ - WARNING - MBD node " << m_mbd_node
                << " not found, continuing with mbdQ = 0 (further warnings will be suppressed)." << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_mbdQ = mbd_node->get_q(0) + mbd_node->get_q(1);
  if ( std::isnan(m_mbdQ) )
  {
    m_mbdQ = 0;
  }

  if (Verbosity() > 1)
  {
    std::cout << "SubtractTowersRhov1::grab_mbdQ - finished grab_mbdQ with m_mbdQ = " << m_mbdQ << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowersRhov1::process_event(PHCompositeNode *topNode)
{

  if (Verbosity() > 1)
  {
    std::cout << "SubtractTowersRhov1::process_event - starting event processing" << std::endl;
  }

  auto * rho_node = findNode::getClass<TowerRhov1 >(topNode, m_rhoNode);
  if ( !rho_node )
  {
    std::cout << PHWHERE << "TowerRho node " << m_rhoNode << " not found, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto res = grab_zvrtx( topNode );
  if ( res != Fun4AllReturnCodes::EVENT_OK )
  {
    if ( Verbosity() > 0 )
    {
      std::cout << "SubtractTowersRhov1::process_event - grab_zvrtx failed, skipping tower subtraction for this event" << std::endl;
    }
    return res;
  }

  // resolve the (zvertex,mbdQ) calibration bin once per event -- ieta-independent
  m_calib_izbin = -1;
  m_calib_imbd = -1;
  if ( m_use_etaCalib && m_etaCalib_loaded )
  {
    grab_mbdQ( topNode );
    m_calib_izbin = find_bin( m_vtxz, m_calib_zvtx_edges );
    m_calib_imbd = find_bin( m_mbdQ, m_calib_mbdQ_edges );
    if (Verbosity() > 1)
    {
      std::cout << "SubtractTowersRhov1::process_event - eta-shape calibration bin: zvtx=" << m_vtxz
                << " -> izbin=" << m_calib_izbin << ", mbdQ=" << m_mbdQ << " -> imbd=" << m_calib_imbd << std::endl;
    }
  }

  auto rho_val = rho_node -> get_rho();
  // const double MULT_THRES_VAL = TMath::Sqrt(2 * 1.0 );
  const double MULT_THRES_VAL = 0.0;
  for ( const auto & target_node_name : m_targetTowerNodes )
  {

    auto * towerinfos = findNode::getClass<TowerInfoContainer>(topNode, target_node_name.c_str());
    if ( !towerinfos )
    {
      std::cout << PHWHERE << "TowerInfo node " << target_node_name << " not found, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    auto output_name = get_outputTowerNode(target_node_name, m_subSuffix);

    auto * sub_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, output_name.c_str());
    if ( !sub_towerinfos )
    {
      std::cout << PHWHERE << "Subtracted TowerInfo node " << output_name << " not found, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    bool is_emcal = (target_node_name.find("CEMC") != std::string::npos);
    bool is_ihcal = (target_node_name.find("HCALIN") != std::string::npos);
    bool is_ohcal = (target_node_name.find("HCALOUT") != std::string::npos);
    bool is_retowered = (target_node_name.find("RETOWER") != std::string::npos);
    const int layer_index = get_layer_index( is_emcal, is_ihcal, is_ohcal );

    std::string geo_node_name = is_emcal ? (is_retowered ? "TOWERGEOM_HCALIN" : "TOWERGEOM_CEMC") : (is_ihcal ? "TOWERGEOM_HCALIN" : (is_ohcal ? "TOWERGEOM_HCALOUT" : ""));
    RawTowerDefs::CalorimeterId calo_id = is_emcal ? ( is_retowered ? RawTowerDefs::CalorimeterId::HCALIN : RawTowerDefs::CalorimeterId::CEMC ) : (is_ihcal ? RawTowerDefs::CalorimeterId::HCALIN :  (is_ohcal ? RawTowerDefs::CalorimeterId::HCALOUT : RawTowerDefs::CalorimeterId::NONE) );

    auto * geom = findNode::getClass<RawTowerGeomContainer>(topNode, geo_node_name.c_str());
    if ( !geom )
    {
      std::cout << PHWHERE << "RawTowerGeomContainer node " << geo_node_name << " not found, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    auto n_eta = geom -> get_etabins();
    std::vector< double > dA_vals (n_eta, 0);
    for (auto ieta = 0; ieta < n_eta; ++ieta)
    {
      auto etabounds = geom -> get_etabounds(ieta);
      auto phibounds = geom -> get_phibounds(0);  // assume constant phi binning
      double dphi = phibounds.second - phibounds.first;
      double deta = etabounds.second - etabounds.first;
      dA_vals.at(ieta) = deta * dphi;
    }

    double calo_radius = 0;
    if ( is_emcal && is_retowered )
    {
      // this needs to load the emcal geo
      const RawTowerDefs::keytype EMCal_key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, 0);
      auto emcal_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
      if ( !emcal_geom )
      {
        std::cout << PHWHERE << "RawTowerGeomContainer node TOWERGEOM_CEMC not found, cannot determine calo radius for flow modulation, doing nothing." << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      auto emcal_tower_geom = emcal_geom -> get_tower_geometry(EMCal_key);
      assert(emcal_tower_geom);
      calo_radius = emcal_tower_geom -> get_center_radius();
    }
    else
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(calo_id, 0, 0);
      auto tower_geom = geom -> get_tower_geometry(key);
      assert(tower_geom);
      calo_radius = tower_geom -> get_center_radius();
    }

    for ( auto ich = 0; ich < (int)towerinfos -> size(); ++ich)
    {
      auto tower = towerinfos -> get_tower_at_channel(ich);
      assert(tower);

      auto tower_key = towerinfos -> encode_key(ich);
      auto ieta = towerinfos -> getTowerEtaBin(tower_key);
      auto iphi = towerinfos -> getTowerPhiBin(tower_key);
      const RawTowerDefs::keytype geo_key = RawTowerDefs::encode_towerid(calo_id, ieta, iphi);

      auto raw_E = tower -> get_energy();
      if ( !tower -> get_isGood() || std::isnan(raw_E) )
      {
        // don't waste time
        sub_towerinfos -> get_tower_at_channel(ich) -> set_energy(0);
        sub_towerinfos -> get_tower_at_channel(ich) -> set_time(tower -> get_time());
        sub_towerinfos -> get_tower_at_channel(ich) -> set_status(tower -> get_status());
        continue;
      }

      auto tower_geom = geom -> get_tower_geometry(geo_key);
      assert(tower_geom);

      auto eta0 = tower_geom -> get_eta();
      auto z0 =  sinh(eta0) * calo_radius;
      auto dz = z0 - m_vtxz;
      double eta = asinh( dz / calo_radius);

      const double w = get_etaWeight( layer_index, ieta );
      double UE = rho_val * cosh(eta) * w;
      if ( m_rho_method == TowerRho::MULT )
      {
        double eT = raw_E / cosh(eta);
        // we  don't apply the multiplicity threshold for now, since it is not clear that it is needed for the current rho calculation.
        if ( eT > MULT_THRES_VAL * rho_val * w ) // mult thres is 0
        {
          // don't subtract negative towers
          UE = 0;
        }
      }
      else if ( m_rho_method == TowerRho::AREA )
      {
        UE *= dA_vals.at(ieta);
      }

      double sub_E = raw_E - UE;

      sub_towerinfos -> get_tower_at_channel(ich) -> set_energy(sub_E);
      sub_towerinfos -> get_tower_at_channel(ich) -> set_time(tower -> get_time());
      sub_towerinfos -> get_tower_at_channel(ich) -> set_status(tower -> get_status());
      if (Verbosity() > 5)
      {
        std::cout << "SubtractTowersRhov1::process_event - tower " << ich << " at ieta / iphi = " << ieta << " / " << iphi << ", raw E = " << raw_E << ", w = " << w << ", UE = " << UE << ", subtracted E = " << sub_E << std::endl;
      }
    }

    if (Verbosity() > 1)
    {
      std::cout << "SubtractTowersRhov1::process_event - finished processing target node " << target_node_name << std::endl;
    }

  }

  if (Verbosity() > 1)
  {
    std::cout << "SubtractTowersRhov1::process_event - finished event processing" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;

}

int SubtractTowersRhov1::CreateNode( PHCompositeNode *topNode )
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if ( !dstNode )
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto * rho_node = findNode::getClass<TowerRhov1 >(topNode, m_rhoNode);
  if ( !rho_node )
  {
    std::cout << PHWHERE << "TowerRho node " << m_rhoNode << " not found, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if ( m_subSuffix.empty() )
  {
    m_subSuffix = rho_node -> get_method_string(rho_node -> get_method());
  }

  m_rho_method = rho_node -> get_method ();

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowersRhov1::CreateNode - using rho method " << m_rho_method << " with suffix " << m_subSuffix << std::endl;
  }

  // loop over target nodes and create new tower nodes for the subtracted towers, if they don't already exist
  for (const auto & target_node_name : m_targetTowerNodes)
  {

    bool is_emcal = (target_node_name.find("CEMC") != std::string::npos);
    bool is_ihcal = (target_node_name.find("HCALIN") != std::string::npos);
    bool is_ohcal = (target_node_name.find("HCALOUT") != std::string::npos);

    std::string target_node_query_str = is_emcal ? "CEMC" : (is_ihcal ? "HCALIN" : (is_ohcal ? "HCALOUT" : ""));
    if ( target_node_query_str.empty() )
    {
      std::cout << PHWHERE << "Target node name " << target_node_name << " does not match expected calorimeter types, skipping." << std::endl;
      continue;
    }

    auto * towerinfos = findNode::getClass<TowerInfoContainer>(topNode, target_node_name.c_str());
    if ( !towerinfos )
    {
      std::cout << PHWHERE << "TowerInfo node " << target_node_name << " not found, doing nothing." << std::endl;
      continue;
    }

    auto * calonode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", target_node_query_str.c_str()));
    if ( !calonode )
    {
      std::cout << PHWHERE << "Calorimeter Node " << target_node_query_str << " not found, doing nothing." << std::endl;
      continue;
    }

    auto output_name = get_outputTowerNode(target_node_name, m_subSuffix);

    auto * test_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, output_name.c_str());
    if ( !test_towerinfos )
    {
      if (Verbosity() > 0)
      {
        std::cout << "SubtractTowersRhov1::CreateNode : creating " << output_name << " node " << std::endl;
      }
      auto * new_towerinfos =  dynamic_cast<TowerInfoContainer *>( towerinfos->CloneMe() );
      auto * towerInfoNode = new PHIODataNode<PHObject>(new_towerinfos, output_name.c_str(), "PHObject");
      calonode -> addNode( towerInfoNode );
    }
    else
    {
      std::cout << "SubtractTowersRhov1::CreateNode : " << output_name << " already exists! " << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowersRhov1::CreateNode - finished creating nodes" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
