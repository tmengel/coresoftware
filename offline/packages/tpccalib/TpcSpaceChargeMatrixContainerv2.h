#ifndef TPCCALIB_TpcSpaceChargeMatrixContainerv2_H
#define TPCCALIB_TpcSpaceChargeMatrixContainerv2_H

/**
 * @file tpccalib/TpcSpaceChargeMatrixContainer.h
 * @author Hugo Pereira Da Costa
 * @date June 2018
 * @brief Contains matrices needed for space charge trackbase reconstruction
 */

#include "TpcSpaceChargeMatrixContainer.h"

#include <array>

/**
 * @brief Cluster container object
 */
class TpcSpaceChargeMatrixContainerv2 : public TpcSpaceChargeMatrixContainer
{
 public:
  /// constructor
  TpcSpaceChargeMatrixContainerv2();

  /// destructor
  ~TpcSpaceChargeMatrixContainerv2() override = default;

  ///@name accessors
  //@{

  /// identify object
  void identify(std::ostream& os = std::cout) const override;

  /// get grid dimensions
  void get_grid_dimensions(int& phibins, int& rbins, int& zbins) const override;

  /// get grid size
  int get_grid_size() const override;

  /// get grid index for given sub-indexes
  int get_cell_index(int iphibin, int irbin, int izbin) const override;

  /// get entries for a given cell
  int get_entries(int cell_index) const override;

  /// get left hand side
  float get_lhs(int cell_index, int i, int j) const override;

  /// get right hand side
  float get_rhs(int cell_index, int i) const override;

  /// get reduced rphi left hand side
  float get_lhs_rphi(int cell_index, int i, int j) const override;

  /// get reduced rphi right hand side
  float get_rhs_rphi(int cell_index, int i) const override;

  /// get reduced z left hand side
  float get_lhs_z(int cell_index, int i, int j) const override;

  /// get reduced z right hand side
  float get_rhs_z(int cell_index, int i) const override;

  //@}

  ///@name modifiers
  //@{

  /// reset method
  void Reset() override;

  /// set grid dimensions
  /**
  \param phibins the number of bins in the azimuth direction
  \param zbins the number of bins along z
  */
  void set_grid_dimensions(int phibins, int rbins, int zbins) override;

  /// increment cell entries
  void add_to_entries(int cell_index) override
  {
    add_to_entries(cell_index, 1);
  }

  /// increment cell entries
  void add_to_entries(int cell_index, int value) override;

  /// increment left hand side matrix
  void add_to_lhs(int cell_index, int i, int j, float value) override;

  /// increment right hand side column
  void add_to_rhs(int cell_index, int i, float value) override;

  /// increment left hand side reduced rphi matrix
  void add_to_lhs_rphi(int cell_index, int i, int j, float value) override;

  /// increment right hand side reduced rphi column
  void add_to_rhs_rphi(int cell_index, int i, float value) override;

  /// increment left hand side reduced rphi matrix
  void add_to_lhs_z(int cell_index, int i, int j, float value) override;

  /// increment right hand side reduced rphi column
  void add_to_rhs_z(int cell_index, int i, float value) override;

  /// add content from other container
  bool add(const TpcSpaceChargeMatrixContainer& other) override;

  //@}

 private:
  /// boundary check
  bool bound_check(int cell_index) const;

  /// boundary check
  bool bound_check(int cell_index, int i) const;

  /// boundary check
  bool bound_check(int cell_index, int i, int j) const;

  /// map matrix index to flat array
  int get_flat_index(int i, int j) const;

  /// boundary check
  bool bound_check_reduced(int cell_index, int i) const;

  /// boundary check
  bool bound_check_reduced(int cell_index, int i, int j) const;

  /// map matrix index to flat array
  int get_flat_index_reduced(int i, int j) const;

  ///@name grid size
  //@{
  int m_phibins = 36;
  int m_rbins = 16;
  int m_zbins = 80;
  //@}

  //@name full 3D matrices (drphi, dz and dr)
  //@{
  //! number of coordinates for full matrix (drphi, dz and dr)
  static constexpr int m_ncoord = 3;

  /// internal matrix representation
  using matrix_t = std::array<float, m_ncoord * m_ncoord>;
  using column_t = std::array<float, m_ncoord>;

  /// left hand side matrices for distortion inversions
  std::vector<matrix_t> m_lhs;

  /// right hand side matrices for distortion inversions
  std::vector<column_t> m_rhs;
  //@}


  //@name reduced 2D matrices (drphi, dr) and (dz and dr)
  //@{
  //! number of coordinates for reduce matrix (drphi, dr) and (dz and dr)
  static constexpr int m_ncoord_reduced = 2;

  /// internal matrix representation
  using reduced_matrix_t = std::array<float, m_ncoord_reduced * m_ncoord_reduced>;
  using reduced_column_t = std::array<float, m_ncoord_reduced>;

  /// left hand side matrices for reduced rphi distortion inversions
  std::vector<reduced_matrix_t> m_lhs_rphi;

  /// right hand side matrices for distortion rphi inversions
  std::vector<reduced_column_t> m_rhs_rphi;

  /// left hand side matrices for reduced z distortion inversions
  std::vector<reduced_matrix_t> m_lhs_z;

  /// right hand side matrices for reduced z distortion inversions
  std::vector<reduced_column_t> m_rhs_z;
  //@}

  /// keep track of how many entries are used per cells
  std::vector<int> m_entries;

  ClassDefOverride(TpcSpaceChargeMatrixContainerv2, 1)
};

#endif
