#include "AverageESub.h"
#include "TowerRhov1.h"
#include "TowerBackgroundv2.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <mbd/MbdOutV2.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TTree.h>


int AverageESub::InitRun( PHCompositeNode *topNode )
{

  auto * fin = new TFile( m_ue_directpath.c_str() , "READ" );
  if ( !fin || fin -> IsZombie() )
  {
    std::cerr << "Error: Could not open input file: " << m_ue_directpath << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  auto * ttree = dynamic_cast<TTree*>( fin -> Get( "T" ) );
  if ( !ttree )
  {
    std::cerr << "Error: Could not find TTree 'cdbobj' in input file: " << m_ue_directpath << std::endl;
    fin -> Close();
    return Fun4AllReturnCodes::ABORTRUN;
  }
  ttree -> SetBranchAddress( "avg_tower_cemc_e", &m_avg_tower_cemc_e );
  ttree -> SetBranchAddress( "avg_tower_hcalin_e", &m_avg_tower_hcalin_e );
  ttree -> SetBranchAddress( "avg_tower_hcalout_e", &m_avg_tower_hcalout_e );
  ttree -> GetEntry( 0 );
  if ( Verbosity() > 1 )
  {
    std::cout << "AverageESub::InitRun: loaded UE parameters from " << m_ue_directpath << std::endl;
  } 
  if ( Verbosity() > 2 )
  {
    std::cout << "AverageESub::InitRun: loaded avg_tower_cemc_e with size " << m_avg_tower_cemc_e -> size() << std::endl;
    std::cout << "AverageESub::InitRun: loaded avg_tower_hcalin_e with size " << m_avg_tower_hcalin_e -> size() << std::endl;
    std::cout << "AverageESub::InitRun: loaded avg_tower_hcalout_e with size " << m_avg_tower_hcalout_e -> size() << std::endl;
    for ( size_t i = 0; i < m_avg_tower_cemc_e -> size(); ++i )
    {
      std::cout << "AverageESub::InitRun: avg_tower_cemc_e[" << i << "] = " << m_avg_tower_cemc_e -> at(i) << std::endl;
    }
    for ( size_t i = 0; i < m_avg_tower_hcalin_e -> size(); ++i )
    {
      std::cout << "AverageESub::InitRun: avg_tower_hcalin_e[" << i << "] = " << m_avg_tower_hcalin_e -> at(i) << std::endl;
    }
    for ( size_t i = 0; i < m_avg_tower_hcalout_e -> size(); ++i )
    {
      std::cout << "AverageESub::InitRun: avg_tower_hcalout_e[" << i << "] = " << m_avg_tower_hcalout_e -> at(i) << std::endl;
    }
  }
  return CreateNode(topNode);
}
  
int AverageESub::process_event(PHCompositeNode *topNode)
{
  
  
  if ( Verbosity() > 0 )
  {
    std::cout << "AverageESub::process_event: start" << std::endl;
  }
  
  m_zvertex = 0;
  m_mbdQ = 0;
  
  
  GlobalVertex * vtx { nullptr };
  auto * vertexmap = findNode::getClass<GlobalVertexMap>( topNode, "GlobalVertexMap" );
  if ( !vertexmap  ) 
  {
    std::cout << PHWHERE << "GlobalVertexMap node missing, skipping event." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if ( vertexmap->empty() ) 
  {
    std::cout << PHWHERE << "GlobalVertexMap node empty, skipping event." << std::endl;
  }

  auto vertices = vertexmap -> get_gvtxs_with_type( { GlobalVertex::MBD } );
  if( !vertices.empty() )
  {
    vtx = vertices.at(0);
  }
  else 
  {
    vtx = vertexmap->begin()->second;
  }
  if ( vtx )
  {
    m_zvertex = vtx->get_z();
  }
    
  if ( std::isnan(m_zvertex) || std::abs(m_zvertex) > 1e3 )
  {
    static bool z_warning_once = true;
    if ( z_warning_once )
    {
        z_warning_once = false;
        std::cout << PHWHERE << " vertex z is " << m_zvertex << ", skipping event (further warnings will be suppressed)." << std::endl;
    }
    m_zvertex = 0;      
  }  

  if ( !m_overlay_node.empty() ) 
  { 

    auto * overlayinfo = findNode::getClass<TowerRhov1>( topNode, m_overlay_node );
    if ( !overlayinfo ) 
    {
        std::cout << PHWHERE << " Input node " << m_overlay_node << " Node missing, doing nothing." << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }
    
    auto m_overlay_zvrtx = overlayinfo -> get_sigma();
    auto m_overlay_mbd = overlayinfo -> get_rho();
    if ( Verbosity() > 1 ) 
    {
      std::cout << PHWHERE << " - overlay_zvrtx = " << m_overlay_zvrtx << std::endl;
      std::cout << PHWHERE << " - overlay_mbd_q_N = " << m_overlay_mbd << std::endl;
    }

    m_mbdQ += m_overlay_mbd;
    m_zvertex = m_overlay_zvrtx;

  }

  // still getmbdQ from MbdOut node if it exists (even w/ overlay node) to be consistent with DetermineSeedlessTowerBackground
  auto * mbd_node = findNode::getClass< MbdOutV2 >( topNode, "MbdOut" );
  if ( !mbd_node ) 
  {
    std::cout << PHWHERE << "MbdOut node missing, skipping event." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // m_mbdQ = mbd_node -> get_q(0) + mbd_node -> get_q(1);
  m_mbdQ += mbd_node -> get_q(0);
  m_mbdQ += mbd_node -> get_q(1);
  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - sum_mbdQ = " << m_mbdQ << std::endl;
    std::cout << PHWHERE << " - zvtx = " << m_zvertex << std::endl;
  }


  int z_bin = get_zbin( m_zvertex );
  int mbd_bin = get_mbd_bin( m_mbdQ );
  
  if ( z_bin < 0 || z_bin >= k_n_zvertex_bins )
  {
    std::cout << PHWHERE << " z_bin = " << z_bin << " out of range, skipping event." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if ( mbd_bin < 0 || mbd_bin >= k_n_mbd_bins )
  {
    std::cout << PHWHERE << " mbd_bin = " << mbd_bin << " out of range, skipping event." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  std::vector< float > this_avg_tower_e;
  int layer_idx = 0;

  auto * tower_background = findNode::getClass<TowerBackgroundv2>( topNode, m_background_node );
  if ( !tower_background )
  {
    std::cout << PHWHERE << " TowerBackgroundv2 node is missing, skipping." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }

  for (const auto& calo : {std::string("CEMC_RETOWER"),
                         std::string("HCALIN"),
                         std::string("HCALOUT")})
  { 
    if ( calo == "CEMC_RETOWER" )
    {
      this_avg_tower_e = *m_avg_tower_cemc_e;
      layer_idx = 0;
    }
    else if ( calo == "HCALIN" )
    {
      this_avg_tower_e = *m_avg_tower_hcalin_e;
      layer_idx = 1;
    }
    else if ( calo == "HCALOUT" )
    {
      this_avg_tower_e = *m_avg_tower_hcalout_e;
      layer_idx = 2;
    }
    else
    {
      std::cerr << "AverageESub::process_event: Unknown calorimeter: " << calo << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    std::string unsub_node = m_towernode_prefix + "_" + calo;
    std::string sub_node = m_tower_sub_prefix + "_" + calo;
    auto * unsub_towers = findNode::getClass<TowerInfoContainer>( topNode, unsub_node );
    if ( !unsub_towers )
    {
      std::cout << "AverageESub::process_event: Cannot find node " << unsub_node << std::endl;
      exit(1);
    }
    auto * sub_towers = findNode::getClass<TowerInfoContainer>( topNode, sub_node );
    if ( !sub_towers )
    {
      std::cout << "AverageESub::process_event: Cannot find node " << sub_node << std::endl;
      exit(1);
    }

    auto nch = unsub_towers -> size();
    for ( unsigned int ich = 0; ich < nch; ++ich )
    {
      auto * tower = unsub_towers -> get_tower_at_channel(ich);
      unsigned int towerkey = unsub_towers -> encode_key(ich);
      int ieta = unsub_towers -> getTowerEtaBin(towerkey);
      int iphi = unsub_towers -> getTowerPhiBin(towerkey);

      float unsub_energy = tower -> get_energy();
      float avg_energy = this_avg_tower_e.at( AverageESub::get_index(ieta, iphi, z_bin, mbd_bin) );
      
      float new_energy = unsub_energy - avg_energy;
      int is_used_for_bkg = 1;
      if ( !tower -> get_isGood() )
      {
        new_energy = 0;
        is_used_for_bkg = 0;
        avg_energy = 0;
      }

      tower_background -> set_UE_layer_eta_phi( layer_idx, ieta, iphi, avg_energy );
      tower_background -> set_is_used_for_bkg( layer_idx, ieta, iphi, is_used_for_bkg );

      // update the subtracted tower container with the new energy, time, and status
      unsub_towers -> get_tower_at_channel(ich) -> set_energy(new_energy);
      unsub_towers -> get_tower_at_channel(ich) -> set_time(tower -> get_time());
      // unsub_towers -> get_tower_at_channel(ich) -> set_status(tower -> get_status());
      sub_towers -> get_tower_at_channel(ich) -> set_energy(new_energy);
      sub_towers -> get_tower_at_channel(ich) -> set_time(tower -> get_time());
      
      if ( Verbosity() > 5 )
      {
        std::cout << "AverageESub::process_event: " << calo << " tower at ieta / iphi = " << ieta << " / " << iphi << ", pre-sub / after-sub E = " << unsub_energy << " / " << new_energy << std::endl;
      }

    } // end loop over channels
  }
  
  if (Verbosity() > 0)
  {
    std::cout << "AverageESub::process_event: finished" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AverageESub::CreateNode( PHCompositeNode *topNode )
{

  PHNodeIterator iter(topNode);
  auto * dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if ( !dstNode )
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto * bkgNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JETBACKGROUND"));
  if ( !bkgNode )
  {
    bkgNode = new PHCompositeNode("JETBACKGROUND");
    dstNode -> addNode(bkgNode);
  }

  auto * tower_bkg = findNode::getClass<TowerBackgroundv2>( topNode, m_background_node );
  if ( !tower_bkg )
  {
    tower_bkg = new TowerBackgroundv2();
    PHIODataNode<PHObject> *bkgDataNode = new PHIODataNode<PHObject>(tower_bkg, m_background_node, "PHObject");
    bkgNode -> addNode(bkgDataNode);
  }
  else
  {
    std::cout << PHWHERE << "::ERROR - " << m_background_node << " pre-exists, but should not" << std::endl;
  }

  for ( const auto& calo : {std::string("CEMC"),
                         std::string("HCALIN"),
                         std::string("HCALOUT")})
  {
    auto * caloNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", calo));
    if ( !caloNode )
    {
      std::cout << PHWHERE << calo << " Node missing, doing nothing." << std::endl;
      exit(1);
    }

    std::string unsub_node_name = m_towernode_prefix + "_" + calo;
    std::string sub_node_name = m_tower_sub_prefix + "_" + calo;
    if ( calo == "CEMC" )
    {
      unsub_node_name = m_towernode_prefix + "_CEMC_RETOWER";
      sub_node_name = m_tower_sub_prefix + "_CEMC_RETOWER";
    }

    auto * sub_tower_node = findNode::getClass<TowerInfoContainer>( topNode, sub_node_name );
    if ( !sub_tower_node )
    {
      if ( Verbosity() > 0 )
      {
        std::cout << "AverageESub::CreateNode: creating node " << sub_node_name << std::endl;
      }
      auto * unsub_tower_node = findNode::getClass<TowerInfoContainer>( topNode, unsub_node_name );
      if ( !unsub_tower_node )
      {
        std::cout << "AverageESub::CreateNode: Cannot find unsubtracted tower node " << unsub_node_name << std::endl;
        exit(1);
      }
      TowerInfoContainer * sub_towers = dynamic_cast<TowerInfoContainer *>( unsub_tower_node -> CloneMe() );
      PHIODataNode<PHObject> * subTowerNode = new PHIODataNode<PHObject>( sub_towers, sub_node_name, "PHObject" );
      caloNode -> addNode( subTowerNode );
    }
    else
    {
      std::cout << "AverageESub::CreateNode: node " << sub_node_name << " already exists!" << std::endl;
    }
  }


  return Fun4AllReturnCodes::EVENT_OK;
}

