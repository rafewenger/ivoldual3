/// @file isodual.cpp
/// @brief Marching cubes/hypercubes isosurface generation
/// - Version 0.3.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2026 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/



#include "ijkisocoord.tpp"
#include "ijkisopoly.tpp"
#include "ijkmesh.tpp"
#include "ijktime.tpp"

#include "ijkdual.tpp"
#include "ijkdual_extract.tpp"

#include "isodual.h"
#include "isodual_datastruct.h"

using namespace IJK;
using namespace IJKDUAL;
using namespace ISODUAL;


// **************************************************
// DUAL CONTOURING (HYPERCUBES)
// **************************************************

/// Dual Contouring Algorithm.
void ISODUAL::dual_contouring
(const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue, 
 DUAL_ISOSURFACE & dual_isosurface, DUALISO_INFO & dualiso_info)
{
  const int dimension = dualiso_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = dualiso_data.ScalarGrid().AxisSize();
  const bool allow_multiple_isov = dualiso_data.AllowMultipleIsoVertices();
  const bool flag_collapse = dualiso_data.CollapseFlag();
  const VERTEX_POSITION_METHOD vpos_method = 
    dualiso_data.VertexPositionMethod();
  
  float merge_time = 0.0;
  PROCEDURE_ERROR error("dual_contouring");

  clock_t t_start = clock();

  if (!dualiso_data.Check(error)) { throw error; };

  dual_isosurface.Clear();
  dualiso_info.time.Clear();

  ISO_MERGE_DATA merge_data(dimension, axis_size);

  if (allow_multiple_isov) {
    dual_contouring_multi_isov
      (dualiso_data.ScalarGrid(), isovalue, dualiso_data,
       dual_isosurface.isopoly_vert, dual_isosurface.isopoly_info, 
       dual_isosurface.vertex_coord, merge_data, dualiso_info);
  }
  else {
    dual_contouring_single_isov
      (dualiso_data.ScalarGrid(), isovalue, vpos_method, 
       dual_isosurface.isopoly_vert, dual_isosurface.isopoly_info,
       dual_isosurface.vertex_coord, merge_data, dualiso_info);
  }

  // store times
  clock_t t_end = clock();
  clock2seconds(t_end-t_start, dualiso_info.time.total);
}


// **************************************************
// DUAL CONTOURING (HYPERCUBES)
// **************************************************

// Extract isosurface using Dual Contouring algorithm.
// Allow multiple isosurface vertices per grid cube.
// Returns list of isosurface polytope vertices
//   and list of isosurface vertex coordinates
void ISODUAL::dual_contouring_multi_isov
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, 
 const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
 const DUALISO_DATA_FLAGS & param,
 std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
 GRID_EDGE_ARRAY & dual_edge,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 DUALISO_INFO & dualiso_info)
{
  std::vector<DUAL_ISOVERT> iso_vlist;

  dual_contouring_multi_isov
    (scalar_grid, isovalue, isodual_table, param, isopoly_vert,
     dual_edge, iso_vlist, vertex_coord, merge_data, dualiso_info);
}


// Extract isosurface using Dual Contouring algorithm.
// Allow multiple isosurface vertices per grid cube.
// Version which creates isodual_table.
void ISODUAL::dual_contouring_multi_isov
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, 
 const DUALISO_DATA_FLAGS & param,
 std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
 GRID_EDGE_ARRAY & dual_edge,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 DUALISO_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const bool flag_separate_neg = param.SeparateNegFlag();
  const VERTEX_POSITION_METHOD vpos_method = param.VertexPositionMethod();

  bool flag_always_separate_opposite(true);

  IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG
    isodual_table(dimension, flag_separate_neg, 
                  flag_always_separate_opposite);

  dual_contouring_multi_isov
    (scalar_grid, isovalue, isodual_table, param, isopoly_vert, 
     dual_edge, vertex_coord, merge_data, dualiso_info);
}


// Extract isosurface using Dual Contouring algorithm.
// Single isosurface vertex per grid cube.
// Returns list of isosurface polytope vertices
//   and list of isosurface vertex coordinates
void ISODUAL::dual_contouring_single_isov
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue, 
 const VERTEX_POSITION_METHOD vertex_position_method,
 std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
 GRID_EDGE_ARRAY & dual_edge,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 DUALISO_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const COORD_TYPE center_offset = 0.1;
  PROCEDURE_ERROR error("dual_contouring");
  clock_t t0, t1, t2, t3;

  t0 = clock();

  isopoly_vert.clear();
  vertex_coord.clear();
  dualiso_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> isopoly;
  extract_dual_isopoly
    (scalar_grid, isovalue, isopoly, dual_edge, dualiso_info);

  t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical(isopoly, iso_vlist, isopoly_vert, merge_data);
  t2 = clock();

  if (vertex_position_method == CUBE_CENTER) {
    position_all_dual_isovertices_cube_center
      (scalar_grid, iso_vlist, vertex_coord);
  }
  else {
    // default
    position_all_dual_isovertices_centroid
      (scalar_grid, isovalue, iso_vlist, vertex_coord);
  }
  t3 = clock();

  // store times
  clock2seconds(t1-t0, dualiso_info.time.extract);
  clock2seconds(t2-t1, dualiso_info.time.merge);
  clock2seconds(t3-t2, dualiso_info.time.position);
  clock2seconds(t3-t0, dualiso_info.time.total);
}
