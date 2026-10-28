/// @file ivoldual_query.h
/// @brief Query routines for ivoldual.
/// - Version 0.3.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2026 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 3.0 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IVOLDUAL_QUERY_
#define _IVOLDUAL_QUERY_

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"
#include "ivoldualtable.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

  // *****************************************************************
  // QUERY IN BOX OR PSEUDOBOX
  // *****************************************************************

  /// Return true if vertex is lower left vertex of an isosurface box.
  /// @param[out] box_corners[] Interval volume vertices which form the
  ///    8 corners of the isosurface box.
  bool is_ivolv_lower_left_vertex_of_isosurface_box
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const IVOL_VERTEX_INDEX ivolv,
   IVOL_VERTEX_INDEX box_corner[8]);

  /// Return true if upper, rightmost vertex of cube is surrounded
  ///   by an isosurface pseudobox.
  /// @param[out] cube_list_index[] Index in cube_list 
  ///   of eight cubes containing the grid vertex surrounded by the pseudobox.
  bool is_upper_right_grid_vertex_in_isosurface_pseudobox
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const std::vector<ACTIVE_GRID_CUBE> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const VERTEX_INDEX cube0_list_index,
   VERTEX_INDEX cube_list_index[8]);

  /// Return true if vertex iv equals w.sep_vert for some interval volume 
  ///   vertex w in cube cube_ivolv.
  bool is_vertex_sep_vertex
  (const VERTEX_INDEX iv,
   const ACTIVE_GRID_CUBE & cube_ivolv,
   const DUAL_IVOLVERT_ARRAY & ivolv_list);


  /// Return true if cube contains doubly connected interval volume vertex.
  bool does_cube_contain_doubly_connected_ivol_vertex
  (const ACTIVE_GRID_CUBE & cube_ivolv,
   const DUAL_IVOLVERT_ARRAY & ivolv_list);


  // *****************************************************************
  // GET FUNCTIONS
  // *****************************************************************

  /// Get the index in cube_list of the seven cubes sharing the
  ///   rightmost/upper vertex with cube0_list_index.
  /// - Gets the cubes only if there are interval volume edges
  ///   between interval volume vertices in each pair of cubes.
  /// - Return true if all 7 cubes found.
  /// @param[out] cube_list_index[] Cubes cube0_list_index
  ///     and the seven other cubes adjacent to cube0_list_index.
  bool get_seven_adjacent_cubes
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const std::vector<ACTIVE_GRID_CUBE> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const int cube0_list_index,
   int cube_list_index[8]);


  /// Get the index in cube_list of the adjacent cube in the oriented direction
  ///  if it shares an interval volume edge with cube0_list_index.
  bool get_adjacent_cube_in_oriented_direction
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const std::vector<ACTIVE_GRID_CUBE> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const int cubeA_list_index,
   const int direction,
   const int orientation,
   int & cubeB_list_index);

}

#endif
