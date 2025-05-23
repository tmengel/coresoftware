/**
 * @file trackbase/TpcDefs.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Utility functions for TPC
 */
#ifndef TPC_TPCDEFS_H
#define TPC_TPCDEFS_H

#include "TrkrDefs.h"

#include <cstdint>  // for uint8_t, uint16_t, uint32_t

/**
 * @brief Utility functions for TPC
 *
 * Contains the functions for manipulating the various keys
 * used by the tpc for hits, hit sets, and clusters
 */
namespace TpcDefs
{
  constexpr int NSides = 2;
  constexpr int NSectors = 12;
  constexpr int NRSectors = 3;

  //! in memory representation of TPC ADC data: 10bit ADC value as 16bit signed integer.
  // This is signed to allow pedestal subtraction when needed
  using ADCDataType = int16_t;

  //! in memory representation of BCO clock using standard 64bit sPHENIX time stamp precision
  using BCODataType = uint64_t;

  // max values for pad and time bin
  static constexpr uint16_t MAXPAD __attribute__((unused)) = 1024;
  static constexpr uint16_t MAXTBIN __attribute__((unused)) = 512;

  /**
   * @brief Get the sector id from hitsetkey
   * @param[in] hitsetkey
   * @param[out] sector id
   */
  uint8_t getSectorId(TrkrDefs::hitsetkey key);

  /**
   * @brief Get the sector id from cluskey
   * @param[in] cluskey
   * @param[out] sector id
   */
  uint8_t getSectorId(TrkrDefs::cluskey key);

  /**
   * @brief Get the side from hitsetkey
   * @param[in] hitsetkey
   * @param[out] side
   */
  uint8_t getSide(TrkrDefs::hitsetkey key);

  /**
   * @brief Get the side id from cluskey
   * @param[in] cluskey
   * @param[out] side id
   */
  uint8_t getSide(TrkrDefs::cluskey key);

  /**
   * @brief Get the pad index from hitkey
   * @param[in] hitkey
   * @param[out] pad index
   */
  uint16_t getPad(TrkrDefs::hitkey key);

  /**
   * @brief Get the time bin from hitkey
   * @param[in] hitkey
   * @param[out] time bin
   */
  uint16_t getTBin(TrkrDefs::hitkey key);

  /**
   * @brief Generate a hitkey from a pad index and time bin
   * @param[in] pad Pad index
   * @param[in] tbin Time bin
   * @param[out] hitkey
   */
  TrkrDefs::hitkey genHitKey(const uint16_t pad, const uint16_t tbin);

  /**
   * @brief Generate a hitsetkey for the tpc
   * @param[in] lyr Layer index
   * @param[in] sector Sector index
   * @param[in] side Side index
   * @param[out] hitsetkey
   *
   * Generate a hitsetkey for the tpc. The tracker id is known
   * implicitly and used in the function.
   */
  TrkrDefs::hitsetkey genHitSetKey(const uint8_t lyr, const uint8_t sector, const uint8_t side);

  /**
   * @brief Generate a cluster key from indeces
   * @param[in] lyr Layer index
   * @param[in] sector Sector index
   * @param[in] side Side index
   * @param[in] clusid Cluster id
   * @param[out] cluskey
   */
  TrkrDefs::cluskey genClusKey(const uint8_t lyr, const uint8_t sector, const uint8_t side, const uint32_t clusid);

}  // namespace TpcDefs

#endif  // TPC_TPCDEFS_H
