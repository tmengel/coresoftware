/**
 * @file trackbase/RawHitSet.cc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Implementation of RawHitSet
 */
#include "RawHitSet.h"

namespace
{
  RawHitSet::Vector dummy_vector;
  RawHitSet::VectorTpc2D dummy_vector_2d;
}  // namespace

void RawHitSet::addHit(RawHit* /*unused*/)
{
  return;
}

// void RawHitSet::addTpcHit(unsigned short, RawHit*)
//{ return; }

void RawHitSet::setTpcPhiBins(unsigned short /*unused*/)
{  // std::cout << "Deprecated settpcphibins " << phibins << std::endl;
  return;
}

RawHitSet::ConstRange RawHitSet::getHits() const
{
  return std::make_pair(dummy_vector.cbegin(), dummy_vector.cend());
}
RawHitSet::ConstVecIterator RawHitSet::getHits(int /*phibin*/)
{
  return dummy_vector_2d.begin();
}

// RawHitSet::ConstRange RawHitSet::getTpcHits(unsigned short phibins) const
//{ std::cout << "Deprecated settpcphibins " << phibins << std::endl; return std::make_pair( dummy_vector.cbegin(), dummy_vector.cend() ); }
