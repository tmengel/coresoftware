#include "DetermineSeedlessTowerBackground.h"

#include <jetbackground/TowerBackground.h>
#include <jetbackground/TowerBackgroundv1.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <eventplaneinfo/Eventplaneinfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

#include <centrality/CentralityInfo.h>
#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TLorentzVector.h>

// standard includes
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <set>

DetermineSeedlessTowerBackground::DetermineSeedlessTowerBackground(const std::string &name)
  : SubsysReco(name)
{
  _UE.resize(3, std::vector<float>(1, 0));
}

int DetermineSeedlessTowerBackground::InitRun(PHCompositeNode *topNode)
{

  // probably this should go in a calibration method 

  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 0 ][ 0 ] = 0.2120965196;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 0 ][ 1 ] = 0.05031564211;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 0 ][ 2 ] = -1.113034585e-06;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 1 ][ 0 ] = 0.1817338304;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 1 ][ 1 ] = 0.0472055138;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 1 ][ 2 ] = -6.132067527e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 2 ][ 0 ] = 0.1736861722;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 2 ][ 1 ] = 0.04495695851;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 2 ][ 2 ] = -2.48724152e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 3 ][ 0 ] = 0.1837111248;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 3 ][ 1 ] = 0.04301007009;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 3 ][ 2 ] = 1.186183463e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 4 ][ 0 ] = 0.1799415174;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 4 ][ 1 ] = 0.04127921522;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 4 ][ 2 ] = 3.55174846e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 5 ][ 0 ] = 0.1546017354;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 5 ][ 1 ] = 0.03902324678;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 5 ][ 2 ] = 5.973436416e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 6 ][ 0 ] = 0.1319221946;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 6 ][ 1 ] = 0.03772326526;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 6 ][ 2 ] = 7.004122164e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 7 ][ 0 ] = 0.1128588244;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 7 ][ 1 ] = 0.03607995457;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 7 ][ 2 ] = 8.010000448e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 8 ][ 0 ] = 0.1265889922;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 8 ][ 1 ] = 0.03423885831;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 8 ][ 2 ] = 9.535080443e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 9 ][ 0 ] = 0.1130356916;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 9 ][ 1 ] = 0.03350577856;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 9 ][ 2 ] = 9.921697356e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 10 ][ 0 ] = 0.08728665167;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 10 ][ 1 ] = 0.03396326782;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 10 ][ 2 ] = 9.975970952e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 11 ][ 0 ] = 0.1091489055;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 11 ][ 1 ] = 0.03302033224;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 11 ][ 2 ] = 1.104297145e-06;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 12 ][ 0 ] = 0.104515202;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 12 ][ 1 ] = 0.03374247946;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 12 ][ 2 ] = 1.129639623e-06;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 13 ][ 0 ] = 0.08991978446;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 13 ][ 1 ] = 0.03447238522;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 13 ][ 2 ] = 1.106562314e-06;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 14 ][ 0 ] = 0.1125588615;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 14 ][ 1 ] = 0.03349849942;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 14 ][ 2 ] = 1.063508117e-06;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 15 ][ 0 ] = 0.1415612047;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 15 ][ 1 ] = 0.03380469268;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 15 ][ 2 ] = 1.050287122e-06;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 16 ][ 0 ] = 0.1207641481;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 16 ][ 1 ] = 0.03462512453;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 16 ][ 2 ] = 8.715866025e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 17 ][ 0 ] = 0.1459010516;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 17 ][ 1 ] = 0.03558933523;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 17 ][ 2 ] = 7.97888879e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 18 ][ 0 ] = 0.1725986651;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 18 ][ 1 ] = 0.03727176344;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 18 ][ 2 ] = 7.256919699e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 19 ][ 0 ] = 0.1919291268;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 19 ][ 1 ] = 0.03942669162;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 19 ][ 2 ] = 5.528969804e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 20 ][ 0 ] = 0.2030374424;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 20 ][ 1 ] = 0.04082207905;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 20 ][ 2 ] = 3.831317392e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 21 ][ 0 ] = 0.1691248798;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 21 ][ 1 ] = 0.04307476882;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 21 ][ 2 ] = -7.419568831e-08;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 22 ][ 0 ] = 0.1410257131;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 22 ][ 1 ] = 0.04615530662;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 22 ][ 2 ] = -4.879195626e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 23 ][ 0 ] = 0.1153054715;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 23 ][ 1 ] = 0.05032899584;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 23 ][ 2 ] = -9.595323109e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 0 ][ 0 ] = 0.1977771725;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 0 ][ 1 ] = 0.005918793556;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 0 ][ 2 ] = -7.259353269e-09;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 1 ][ 0 ] = 0.2208571546;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 1 ][ 1 ] = 0.003175225483;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 1 ][ 2 ] = 3.975835829e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 2 ][ 0 ] = 0.2161655328;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 2 ][ 1 ] = 0.002982590709;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 2 ][ 2 ] = 4.626212925e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 3 ][ 0 ] = 0.218479975;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 3 ][ 1 ] = 0.002557665388;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 3 ][ 2 ] = 4.946929979e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 4 ][ 0 ] = 0.2184736923;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 4 ][ 1 ] = 0.002230489696;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 4 ][ 2 ] = 5.098212872e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 5 ][ 0 ] = 0.2254862045;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 5 ][ 1 ] = 0.001944753547;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 5 ][ 2 ] = 5.209088531e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 6 ][ 0 ] = 0.224122218;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 6 ][ 1 ] = 0.001719860659;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 6 ][ 2 ] = 5.058672332e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 7 ][ 0 ] = 0.2271340041;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 7 ][ 1 ] = 0.001554593673;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 7 ][ 2 ] = 5.013579076e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 8 ][ 0 ] = 0.2288534658;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 8 ][ 1 ] = 0.00137120559;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 8 ][ 2 ] = 4.779423453e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 9 ][ 0 ] = 0.2276760881;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 9 ][ 1 ] = 0.001207441324;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 9 ][ 2 ] = 4.474249892e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 10 ][ 0 ] = 0.2236360497;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 10 ][ 1 ] = 0.001236441242;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 10 ][ 2 ] = 4.436846871e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 11 ][ 0 ] = 0.2287048618;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 11 ][ 1 ] = 0.001209408447;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 11 ][ 2 ] = 4.514057368e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 12 ][ 0 ] = 0.2283755648;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 12 ][ 1 ] = 0.001207903426;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 12 ][ 2 ] = 4.500283304e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 13 ][ 0 ] = 0.2233407769;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 13 ][ 1 ] = 0.001235190051;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 13 ][ 2 ] = 4.409357624e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 14 ][ 0 ] = 0.2272496812;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 14 ][ 1 ] = 0.001205374311;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 14 ][ 2 ] = 4.411122806e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 15 ][ 0 ] = 0.2281299548;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 15 ][ 1 ] = 0.00136635695;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 15 ][ 2 ] = 4.710385051e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 16 ][ 0 ] = 0.2288353895;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 16 ][ 1 ] = 0.001555709026;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 16 ][ 2 ] = 4.934539211e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 17 ][ 0 ] = 0.2266255838;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 17 ][ 1 ] = 0.001714542677;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 17 ][ 2 ] = 5.031838854e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 18 ][ 0 ] = 0.2268472562;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 18 ][ 1 ] = 0.001938627877;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 18 ][ 2 ] = 5.151298258e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 19 ][ 0 ] = 0.2199272726;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 19 ][ 1 ] = 0.002218715122;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 19 ][ 2 ] = 4.996123591e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 20 ][ 0 ] = 0.2184220561;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 20 ][ 1 ] = 0.002517138474;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 20 ][ 2 ] = 4.915608206e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 21 ][ 0 ] = 0.2166350649;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 21 ][ 1 ] = 0.002935611121;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 21 ][ 2 ] = 4.6219667e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 22 ][ 0 ] = 0.2206710453;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 22 ][ 1 ] = 0.003137362329;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 22 ][ 2 ] = 3.830413492e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 23 ][ 0 ] = 0.198504869;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 23 ][ 1 ] = 0.005713949191;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 23 ][ 2 ] = 1.057068586e-08;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 0 ][ 0 ] = 0.1464578278;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 0 ][ 1 ] = 0.02535797864;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 0 ][ 2 ] = -1.563768599e-06;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 1 ][ 0 ] = 0.1065649392;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 1 ][ 1 ] = 0.008761897524;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 1 ][ 2 ] = -3.57194577e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 2 ][ 0 ] = 0.100240904;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 2 ][ 1 ] = 0.00999668167;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 2 ][ 2 ] = -3.581594087e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 3 ][ 0 ] = 0.09576232891;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 3 ][ 1 ] = 0.007795809374;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 3 ][ 2 ] = -1.508583266e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 4 ][ 0 ] = 0.07285955413;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 4 ][ 1 ] = 0.008134762215;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 4 ][ 2 ] = 2.266406674e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 5 ][ 0 ] = 0.06576189358;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 5 ][ 1 ] = 0.006192823688;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 5 ][ 2 ] = 2.638589856e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 6 ][ 0 ] = 0.06618350589;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 6 ][ 1 ] = 0.005567834478;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 6 ][ 2 ] = 2.835233455e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 7 ][ 0 ] = 0.06117568515;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 7 ][ 1 ] = 0.00543203112;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 7 ][ 2 ] = 2.816557708e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 8 ][ 0 ] = 0.06011015376;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 8 ][ 1 ] = 0.0050126664;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 8 ][ 2 ] = 2.558503702e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 9 ][ 0 ] = 0.06744393651;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 9 ][ 1 ] = 0.004703943817;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 9 ][ 2 ] = 2.754483032e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 10 ][ 0 ] = 0.06388098862;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 10 ][ 1 ] = 0.004864377698;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 10 ][ 2 ] = 2.860649967e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 11 ][ 0 ] = 0.0658736226;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 11 ][ 1 ] = 0.004821619679;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 11 ][ 2 ] = 2.948491528e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 12 ][ 0 ] = 0.06555143733;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 12 ][ 1 ] = 0.004825562404;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 12 ][ 2 ] = 2.886793527e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 13 ][ 0 ] = 0.06720369769;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 13 ][ 1 ] = 0.004858058152;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 13 ][ 2 ] = 2.775937661e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 14 ][ 0 ] = 0.06877584293;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 14 ][ 1 ] = 0.004734190352;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 14 ][ 2 ] = 2.534302938e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 15 ][ 0 ] = 0.06356673181;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 15 ][ 1 ] = 0.005059240744;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 15 ][ 2 ] = 2.363434666e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 16 ][ 0 ] = 0.07405824355;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 16 ][ 1 ] = 0.005606746235;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 16 ][ 2 ] = 2.44604726e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 17 ][ 0 ] = 0.0791229347;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 17 ][ 1 ] = 0.005842182651;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 17 ][ 2 ] = 2.262605975e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 18 ][ 0 ] = 0.07646442502;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 18 ][ 1 ] = 0.006518746884;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 18 ][ 2 ] = 1.428283131e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 19 ][ 0 ] = 0.08810613024;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 19 ][ 1 ] = 0.008209039681;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 19 ][ 2 ] = 3.835134755e-08;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 20 ][ 0 ] = 0.09833221924;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 20 ][ 1 ] = 0.007951396113;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 20 ][ 2 ] = -1.804688501e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 21 ][ 0 ] = 0.105225752;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 21 ][ 1 ] = 0.01006848423;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 21 ][ 2 ] = -4.29247696e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 22 ][ 0 ] = 0.1135659493;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 22 ][ 1 ] = 0.009154409329;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 22 ][ 2 ] = -4.481355404e-07;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 23 ][ 0 ] = 0.1857438289;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 23 ][ 1 ] = 0.02782668769;
  _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 23 ][ 2 ] = -2.413541556e-06;

  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 0 ][ 0 ] = 0.0410237;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 0 ][ 1 ] = 0.0607689;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 0 ][ 2 ] = -8.72125e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 1 ][ 0 ] = 0.0342337;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 1 ][ 1 ] = 0.0543281;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 1 ][ 2 ] = -5.85177e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 2 ][ 0 ] = 0.027464;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 2 ][ 1 ] = 0.0509642;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 2 ][ 2 ] = -4.41797e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 3 ][ 0 ] = 0.0223333;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 3 ][ 1 ] = 0.0475657;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 3 ][ 2 ] = -3.56281e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 4 ][ 0 ] = 0.0190667;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 4 ][ 1 ] = 0.0457327;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 4 ][ 2 ] = -2.73986e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 5 ][ 0 ] = 0.0165371;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 5 ][ 1 ] = 0.0422734;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 5 ][ 2 ] = -1.63786e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 6 ][ 0 ] = 0.0143733;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 6 ][ 1 ] = 0.0408131;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 6 ][ 2 ] = -4.30828e-08;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 7 ][ 0 ] = 0.0129414;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 7 ][ 1 ] = 0.0387698;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 7 ][ 2 ] = 2.59592e-08;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 8 ][ 0 ] = 0.0106727;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 8 ][ 1 ] = 0.0373523;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 8 ][ 2 ] = 5.84453e-08;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 9 ][ 0 ] = 0.00581824;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 9 ][ 1 ] = 0.0356099;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 9 ][ 2 ] = 1.3337e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 10 ][ 0 ] = 0.00388854;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 10 ][ 1 ] = 0.0357452;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 10 ][ 2 ] = 1.74325e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 11 ][ 0 ] = -1.92154e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 11 ][ 1 ] = 0.0325906;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 11 ][ 2 ] = 2.13534e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 12 ][ 0 ] = 0.00328004;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 12 ][ 1 ] = 0.0315743;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 12 ][ 2 ] = 5.86105e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 13 ][ 0 ] = 0.00174584;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 13 ][ 1 ] = 0.0324987;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 13 ][ 2 ] = 5.60382e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 14 ][ 0 ] = -0.00143222;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 14 ][ 1 ] = 0.0332299;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 14 ][ 2 ] = 5.27528e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 15 ][ 0 ] = -0.000383052;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 15 ][ 1 ] = 0.0349354;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 15 ][ 2 ] = 4.42469e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 16 ][ 0 ] = -0.00149554;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 16 ][ 1 ] = 0.0357612;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 16 ][ 2 ] = 3.82722e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 17 ][ 0 ] = -0.00324126;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 17 ][ 1 ] = 0.0368239;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 17 ][ 2 ] = 2.78986e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 18 ][ 0 ] = -0.0080648;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 18 ][ 1 ] = 0.0392451;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 18 ][ 2 ] = 1.9644e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 19 ][ 0 ] = -0.0151704;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 19 ][ 1 ] = 0.0405672;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 19 ][ 2 ] = 8.13107e-08;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 20 ][ 0 ] = -0.0262827;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 20 ][ 1 ] = 0.0429652;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 20 ][ 2 ] = -2.68251e-08;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 21 ][ 0 ] = -0.0316013;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 21 ][ 1 ] = 0.0466257;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 21 ][ 2 ] = -1.48407e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 22 ][ 0 ] = -0.0412391;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 22 ][ 1 ] = 0.0498418;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 22 ][ 2 ] = -2.49831e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 23 ][ 0 ] = -0.0764219;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 23 ][ 1 ] = 0.0533237;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 0 ][ 23 ][ 2 ] = -4.21773e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 0 ][ 0 ] = 0.0108872;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 0 ][ 1 ] = 0.0974968;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 0 ][ 2 ] = -6.08305e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 1 ][ 0 ] = 0.000853813;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 1 ][ 1 ] = 0.0541406;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 1 ][ 2 ] = -8.8336e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 2 ][ 0 ] = -0.00131265;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 2 ][ 1 ] = 0.0514461;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 2 ][ 2 ] = 2.35087e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 3 ][ 0 ] = -0.00179621;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 3 ][ 1 ] = 0.04647;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 3 ][ 2 ] = 5.90116e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 4 ][ 0 ] = -0.00192955;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 4 ][ 1 ] = 0.0415699;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 4 ][ 2 ] = 8.4638e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 5 ][ 0 ] = -0.00190525;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 5 ][ 1 ] = 0.0379303;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 5 ][ 2 ] = 9.10542e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 6 ][ 0 ] = -0.00207737;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 6 ][ 1 ] = 0.0359686;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 6 ][ 2 ] = 1.0676e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 7 ][ 0 ] = -0.00196435;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 7 ][ 1 ] = 0.0319346;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 7 ][ 2 ] = 1.061e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 8 ][ 0 ] = -0.00196997;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 8 ][ 1 ] = 0.0302201;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 8 ][ 2 ] = 1.09495e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 9 ][ 0 ] = -0.00189448;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 9 ][ 1 ] = 0.0272469;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 9 ][ 2 ] = 1.13289e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 10 ][ 0 ] = -0.00272024;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 10 ][ 1 ] = 0.0303318;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 10 ][ 2 ] = 1.48984e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 11 ][ 0 ] = -0.00190785;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 11 ][ 1 ] = 0.0272064;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 11 ][ 2 ] = 1.17566e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 12 ][ 0 ] = -0.001932;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 12 ][ 1 ] = 0.0269148;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 12 ][ 2 ] = 1.16099e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 13 ][ 0 ] = -0.00198508;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 13 ][ 1 ] = 0.0268362;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 13 ][ 2 ] = 1.18974e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 14 ][ 0 ] = -0.00195361;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 14 ][ 1 ] = 0.0267548;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 14 ][ 2 ] = 1.1398e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 15 ][ 0 ] = -0.00204988;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 15 ][ 1 ] = 0.0295945;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 15 ][ 2 ] = 1.11054e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 16 ][ 0 ] = -0.00189532;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 16 ][ 1 ] = 0.0314371;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 16 ][ 2 ] = 8.83787e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 17 ][ 0 ] = -0.00203241;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 17 ][ 1 ] = 0.0345381;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 17 ][ 2 ] = 7.98156e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 18 ][ 0 ] = -0.00184849;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 18 ][ 1 ] = 0.036529;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 18 ][ 2 ] = 6.27554e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 19 ][ 0 ] = -0.00185812;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 19 ][ 1 ] = 0.0408127;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 19 ][ 2 ] = 4.97121e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 20 ][ 0 ] = -0.00160055;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 20 ][ 1 ] = 0.0443868;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 20 ][ 2 ] = 3.35423e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 21 ][ 0 ] = -0.000742139;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 21 ][ 1 ] = 0.0493966;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 21 ][ 2 ] = -3.6575e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 22 ][ 0 ] = 0.00217601;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 22 ][ 1 ] = 0.0514536;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 22 ][ 2 ] = -1.87708e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 23 ][ 0 ] = 0.0143985;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 23 ][ 1 ] = 0.089756;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 1 ][ 23 ][ 2 ] = -8.42829e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 0 ][ 0 ] = 0.0421447;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 0 ][ 1 ] = 0.116571;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 0 ][ 2 ] = -8.90003e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 1 ][ 0 ] = 0.0120358;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 1 ][ 1 ] = 0.0491275;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 1 ][ 2 ] = -1.93597e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 2 ][ 0 ] = 0.00732104;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 2 ][ 1 ] = 0.0532094;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 2 ][ 2 ] = -1.08101e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 3 ][ 0 ] = -0.000123433;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 3 ][ 1 ] = 0.0411288;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 3 ][ 2 ] = 6.97794e-07;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 4 ][ 0 ] = -0.00735767;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 4 ][ 1 ] = 0.0442175;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 4 ][ 2 ] = 7.52332e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 5 ][ 0 ] = -0.0125412;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 5 ][ 1 ] = 0.0331119;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 5 ][ 2 ] = 7.18691e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 6 ][ 0 ] = -0.0140661;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 6 ][ 1 ] = 0.0289481;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 6 ][ 2 ] = 6.948e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 7 ][ 0 ] = -0.0158358;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 7 ][ 1 ] = 0.0285886;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 7 ][ 2 ] = 7.39142e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 8 ][ 0 ] = -0.0151552;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 8 ][ 1 ] = 0.027205;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 8 ][ 2 ] = 6.99881e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 9 ][ 0 ] = -0.0143379;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 9 ][ 1 ] = 0.0253338;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 9 ][ 2 ] = 6.74075e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 10 ][ 0 ] = -0.0149106;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 10 ][ 1 ] = 0.0249937;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 10 ][ 2 ] = 6.7351e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 11 ][ 0 ] = -0.0150088;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 11 ][ 1 ] = 0.0243681;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 11 ][ 2 ] = 6.865e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 12 ][ 0 ] = -0.014241;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 12 ][ 1 ] = 0.0240858;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 12 ][ 2 ] = 6.60499e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 13 ][ 0 ] = -0.01437;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 13 ][ 1 ] = 0.0248938;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 13 ][ 2 ] = 6.53552e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 14 ][ 0 ] = -0.0141026;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 14 ][ 1 ] = 0.0249408;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 14 ][ 2 ] = 6.13532e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 15 ][ 0 ] = -0.0142472;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 15 ][ 1 ] = 0.0271072;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 15 ][ 2 ] = 6.06314e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 16 ][ 0 ] = -0.0114381;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 16 ][ 1 ] = 0.0290995;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 16 ][ 2 ] = 6.22165e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 17 ][ 0 ] = -0.00934935;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 17 ][ 1 ] = 0.0298785;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 17 ][ 2 ] = 4.80892e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 18 ][ 0 ] = -0.0068252;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 18 ][ 1 ] = 0.034049;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 18 ][ 2 ] = 2.90508e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 19 ][ 0 ] = 0.00165265;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 19 ][ 1 ] = 0.0447341;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 19 ][ 2 ] = -1.30743e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 20 ][ 0 ] = 0.00804001;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 20 ][ 1 ] = 0.0423475;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 20 ][ 2 ] = -5.68899e-06;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 21 ][ 0 ] = 0.0195586;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 21 ][ 1 ] = 0.0516363;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 21 ][ 2 ] = -1.21752e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 22 ][ 0 ] = 0.0280568;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 22 ][ 1 ] = 0.0513265;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 22 ][ 2 ] = -1.55342e-05;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 23 ][ 0 ] = 0.0788104;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 23 ][ 1 ] = 0.118279;
  // _seedless_UE_parameters_LAYER_ETA_PARAM[ 2 ][ 23 ][ 2 ] = -4.62318e-05;

  
  return CreateNode(topNode);
}

float DetermineSeedlessTowerBackground::get_UE_estimate( int layer, int eta , float layer_sumet ) {

  float p0 = _seedless_UE_parameters_LAYER_ETA_PARAM[ layer ][ eta ][ 0 ];
  float p1 = _seedless_UE_parameters_LAYER_ETA_PARAM[ layer ][ eta ][ 1 ];
  float p2 = _seedless_UE_parameters_LAYER_ETA_PARAM[ layer ][ eta ][ 2 ];

  // note: no checks on good bounds for layer_sumet
  float estimate = p0 + p1 * layer_sumet + p2 * pow( layer_sumet , 2 );
  
  if ( Verbosity() >= 10 ) {
    std::cout << "DetermineSeedlessTowerBackground::get_UE_estimate called on layer / eta = " << layer << " / " << eta << " with SumET = " << layer_sumet << ", returning " << estimate << std::endl;    
  }

  return estimate;
  
}
  
int DetermineSeedlessTowerBackground::process_event(PHCompositeNode *topNode)
{
  
  
  if ( Verbosity() > 0 )
  {
    std::cout << "DetermineSeedlessTowerBackground::process_event: start" << std::endl;
  }

  
  // pull out the tower containers and geometry objects at the start
  EMTowerName = m_towerNodePrefix + "_CEMC_RETOWER";
  IHTowerName = m_towerNodePrefix + "_HCALIN";
  OHTowerName = m_towerNodePrefix + "_HCALOUT";
  auto * towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
  auto * towerinfosIH3 = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  auto * towerinfosOH3 = findNode::getClass<TowerInfoContainer>(topNode, OHTowerName);
  if ( !towerinfosEM3 )
  {
    std::cout << "DetermineSeedlessTowerBackground::process_event: Cannot find node " << EMTowerName << std::endl;
    exit(1);
  }
  if ( !towerinfosIH3 )
  {
    std::cout << "DetermineSeedlessTowerBackground::process_event: Cannot find node " << IHTowerName << std::endl;
    exit(1);
  }
  if ( !towerinfosOH3 )
  {
    std::cout << "DetermineSeedlessTowerBackground::process_event: Cannot find node " << OHTowerName << std::endl;
    exit(1);
  }

  // reset all maps
  _UE.assign(3, std::vector<float>(24, 0));

  // determine total SumET
   
  _layer_sumet[0] = 0;
  for (long unsigned int ch = 0; ch < towerinfosEM3->size(); ch++)
    {
      _layer_sumet[0] += towerinfosEM3->get_tower_at_channel( ch )->get_energy();
    }
  _layer_sumet[1] = 0;
  for (long unsigned int ch = 0; ch < towerinfosIH3->size(); ch++)
    {
      _layer_sumet[1] += towerinfosIH3->get_tower_at_channel( ch )->get_energy();
    }
  _layer_sumet[2] = 0;
  for (long unsigned int ch = 0; ch < towerinfosOH3->size(); ch++)
    {
      _layer_sumet[2] += towerinfosOH3->get_tower_at_channel( ch )->get_energy();
    }

  if ( Verbosity() >= 1 ) {
    std::cout << "DetermineSeedlessTowerBackground::process_event layer SumET = " << _layer_sumet[0] << " / " << _layer_sumet[1] << " / " << _layer_sumet[2] << std::endl;
  }
  
  // fill UE vectors

  for (int layer = 0; layer < 3; layer++) {
    for (int eta = 0; eta < 24; eta++) {
      
      float UE_estimate = get_UE_estimate( layer, eta , _layer_sumet[ layer ] );

      // note: may need some special casing here for, e.g., eta region with fully disabled towers
      UE_estimate /= 64.0;
      
      _UE[layer][eta] = UE_estimate ;
    }
  }
  
  // that's it! 

  // end main code
  
  if (Verbosity() > 0)
  {
    for (int layer = 0; layer < 3; layer++)
    {
      std::cout << "DetermineSeedlessTowerBackground::process_event: summary UE in layer " << layer << " : ";
      for (int eta = 0; eta < 24; eta++)
      {
        std::cout << _UE[layer].at(eta) << " , ";
      }
      std::cout << std::endl;
    }
  }

  FillNode(topNode);

  if (Verbosity() > 0)
  {
    std::cout << "DetermineSeedlessTowerBackground::process_event: exiting" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int DetermineSeedlessTowerBackground::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the jet background stuff under a sub-node directory
  PHCompositeNode *bkgNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JETBACKGROUND"));
  if (!bkgNode)
  {
    bkgNode = new PHCompositeNode("JETBACKGROUND");
    dstNode->addNode(bkgNode);
  }

  // create the TowerBackground node...
  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, _backgroundName);
  if (!towerbackground)
  {
    towerbackground = new TowerBackgroundv1();
    PHIODataNode<PHObject> *bkgDataNode = new PHIODataNode<PHObject>(towerbackground, _backgroundName, "PHObject");
    bkgNode->addNode(bkgDataNode);
  }
  else
  {
    std::cout << PHWHERE << "::ERROR - " << _backgroundName << " pre-exists, but should not" << std::endl;
    // exit(-1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void DetermineSeedlessTowerBackground::FillNode(PHCompositeNode *topNode)
{
  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, _backgroundName);
  if (!towerbackground)
  {
    std::cout << " ERROR -- can't find TowerBackground node after it should have been created" << std::endl;
    return;
  }

  towerbackground->set_UE(0, _UE[0]);
  towerbackground->set_UE(1, _UE[1]);
  towerbackground->set_UE(2, _UE[2]);

  // do not use these in seedless
  towerbackground->set_v2( 0 );
  towerbackground->set_Psi2( 0 );
  towerbackground->set_nStripsUsedForFlow( 0 );
  towerbackground->set_nTowersUsedForBkg( 0 );
  towerbackground->set_flow_failure_flag( 0 );

  return;
}
