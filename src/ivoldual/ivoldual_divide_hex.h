/// @file ivoldual_divide_hex.h
/// @brief Subdivide polytopes.
/// - Version 0.3.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2026 Rephael Wenger

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

#ifndef _IVOLDUAL_DIVIDE_HEX_
#define _IVOLDUAL_DIVIDE_HEX_

#include "ijk.tpp"

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"
#include "ivoldualtable.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

  /// Split hexahedron.
  void split_hex
  (std::vector<VERTEX_INDEX> & ivolpoly_vert,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord, 
   COORD_TYPE jacobian_limit);

  /// Collapse hexahedron.
  void collapse_hex
  (std::vector<VERTEX_INDEX> & ivolpoly_vert,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord, 
   COORD_TYPE jacobian_limit);

  void subdivide_hex_to_four
  (std::vector<VERTEX_INDEX> & ivolpoly_vert,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   std::vector<std::pair<VERTEX_INDEX, VERTEX_INDEX>> &
   	vertex_subdivide_list,
   COORD_ARRAY & vertex_coord);

  void generate_children_polytope
  (std::vector<VERTEX_INDEX> & ivolpoly_vert,
	int iv[],	int iw[], int corner[]);
}

#endif
