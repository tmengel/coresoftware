#ifndef JETBACKGROUND_RHOUEPROFILE_H
#define JETBACKGROUND_RHOUEPROFILE_H

#include <cmath>

/// Shared arithmetic for the rho-based underlying-event profile.
///
/// Two independent paths subtract the same UE from the same towers:
///   SubtractTowersRhov1                         (computes UE per tower)
///   DetermineTowerBackgroundv1 -> SubtractTowers (stores UE per eta strip)
/// The second path stores the profile in TowerBackground as std::vector<float>,
/// so it can never carry more than float precision. For the two to produce
/// bit-identical towers -- and therefore identical jets -- every intermediate
/// must be computed with the same types, in the same order, by the same code.
/// Both classes call these functions; do not reimplement them inline.
namespace RhoUEProfile
{
  /// Tower eta corrected for the event vertex, rounded to float (the precision
  /// TowerBackground carries).
  inline float corrected_eta(const double tower_eta, const float radius, const float vtxz)
  {
    const double r = radius;
    const double z = (std::sinh(tower_eta) * r) - static_cast<double>(vtxz);
    return static_cast<float>(std::asinh(z / r));
  }

  /// UE energy per tower: rho * cosh(eta) * w, rounded to float.
  inline float ue(const float rho, const float eta_corr, const float w)
  {
    return static_cast<float>(static_cast<double>(rho) * std::cosh(static_cast<double>(eta_corr)) *
                              static_cast<double>(w));
  }
}  // namespace RhoUEProfile

#endif
