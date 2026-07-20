#include "TowerBackgroundv2.h"

#include <iostream>

TowerBackgroundv2::TowerBackgroundv2()
{
  _UE.resize(3, std::vector<float>(1, 0));
  _UE_layer_eta_phi.resize(3, std::vector<std::vector<float> >(24, std::vector<float>(64, 0)));
  _is_used_for_bkg.resize(3, std::vector<std::vector<int> >(24, std::vector<int>(64, 0)));

}

void TowerBackgroundv2::identify(std::ostream& os) const
{
  os << "TowerBackground: " << std::endl;
  for (int n = 0; n < 3; n++)
  {
    os << " layer " << n << " : UE in " << _UE[n].size() << " eta bins: ";
    for (float eta : _UE[n])
    {
      os << eta << " ";
    }
    os << std::endl;
  }

  os << " v2 = " << _v2 << ", Psi2 = " << _Psi2
     << ", # towers used for bkg = " << _nTowers
     << " , # strips used for flow = " << _nStrips
     << " , flow failure flag " << _flow_failure_flag << std::endl;

  return;
}
