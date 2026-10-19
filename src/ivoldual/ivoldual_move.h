/// @file ivoldual_move.h
/// @brief Routines for moving/repositioning interval volume vertices.
/// - Version 0.3.0


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

#ifndef _IVOLDUAL_MOVE_
#define _IVOLDUAL_MOVE_

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"
#include "ivoldualtable.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

  // *****************************************************************
  // MOVE VERTEX TO CUBE CENTER
  // *****************************************************************

  /// Move interval volume vertex to cube center.
  void move_ivol_vertex_to_cube_center
  (const std::vector<GRID_CUBE_DATA> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const IVOL_VERTEX_INDEX ivolv,
   COORD_ARRAY & vertex_coord);


  // *****************************************************************
  // MOVE POINT/VERTEX AWAY FROM CUBE FACETS
  // *****************************************************************

  /// Move point so it is at least distance from a given cube facet.
  /// @pre 0 <= distance <= 0.5.
  void move_point_away_from_cube_facet_3D
  (const GRID_COORD_TYPE cube_coord[], const int jfacet,
   const COORD_TYPE distance, COORD_TYPE point_coord[], bool & flag_moved);

  /// Move interval volume vertex so it is at least distance 
  ///   from grid cube facet jfacet.
  void move_ivol_vertex_away_from_cube_facet
  (const std::vector<GRID_CUBE_DATA> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const IVOL_VERTEX_INDEX ivolv, const int jfacet, const COORD_TYPE distance, 
   COORD_ARRAY & vertex_coord, bool & flag_moved);

  /// Move point so it is at least distance from all cube facets.
  /// @param[out] flag_moved True, if vertex coordinates are changed.
  /// @pre 0 <= distance <= 0.5.
  void move_point_away_from_all_cube_facets_3D
  (const GRID_COORD_TYPE cube_coord[], const COORD_TYPE distance, 
   COORD_TYPE point_coord[], bool & flag_moved);

  /// Move interval volume vertex so it is at least distance 
  ///   from all grid cube facets.
  void move_ivol_vertex_away_from_all_cube_facets
  (const std::vector<GRID_CUBE_DATA> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const IVOL_VERTEX_INDEX ivolv, const COORD_TYPE distance, 
   COORD_ARRAY & vertex_coord, bool & flag_moved);


  // *****************************************************************
  // MOVE POINT/VERTEX IN BOUNDARY CUBE AWAY FROM GRID INTERIOR
  // *****************************************************************

  /// Move point so that it is at least distance from cubes in grid interior.
  /// @pre 0 <= distance <= 1.
  void move_point_away_from_grid_interior_3D
  (const DUALISO_GRID & grid, 
   const GRID_COORD_TYPE cube_coord[], const COORD_TYPE distance, 
   COORD_TYPE point_coord[], bool & flag_moved);

  /// Move interval volume vertex so that it is at least distance 
  ///   from cubes in grid interior.
  /// @pre Vertex in in grid boundary cube.
  void move_ivol_vertex_away_from_grid_interior
  (const DUALISO_GRID & grid, 
   const std::vector<GRID_CUBE_DATA> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   const IVOL_VERTEX_INDEX ivolv, const COORD_TYPE distance, 
   COORD_ARRAY & vertex_coord, bool & flag_moved);


  // *****************************************************************
  // EXPAND THIN REGIONS
  // *****************************************************************

  /// Expand thin regions.
  /// Move isosurface vertices in thin regions away from grid cube facets.
  void expand_thin_regions
  (const DUALISO_GRID & grid, const std::vector<GRID_CUBE_DATA> & cube_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list, 
   const COORD_TYPE separation_distance,
   COORD_ARRAY & vertex_coord,
   int & num_moved);

}

#endif

