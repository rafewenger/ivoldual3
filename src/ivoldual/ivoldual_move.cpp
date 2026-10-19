/// @file ivoldual_move.cpp
/// @brief Routines for moving/repositioning interval volume vertices.
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


#include "ivoldual_move.h"

#include "ijkdual_query.tpp"

// *****************************************************************
// MOVE VERTEX TO CUBE CENTER
// *****************************************************************

// Move interval volume vertex to cube center.
void IVOLDUAL::move_ivol_vertex_to_cube_center
(const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const IVOL_VERTEX_INDEX ivolv, 
 COORD_ARRAY & vertex_coord)
{
  const int DIM3(3);
  const GRID_COORD_TYPE * cube_coord =
    IJKDUAL::get_coord_of_cube_containing_isov(cube_list, ivolv_list, ivolv);

  for (int d = 0; d < DIM3; d++) 
    { vertex_coord[ivolv*DIM3+d] = cube_coord[d] + 0.5; }
}


// *****************************************************************
// MOVE POINT/VERTEX AWAY FROM CUBE FACETS
// *****************************************************************

// Move point so it is at least distance from cube facet jfacet.
// @pre 0 <= distance <= 0.5.
void IVOLDUAL::move_point_away_from_cube_facet_3D
(const GRID_COORD_TYPE cube_coord[], const int jfacet,
 const COORD_TYPE distance, COORD_TYPE point_coord[], bool & flag_moved)
{
  const int DIM3(3);
  const int d = jfacet%DIM3;
  const GRID_COORD_TYPE facet_coord = cube_coord[d];
  const COORD_TYPE c = point_coord[d];

  flag_moved = false;
  if (jfacet < DIM3) {
    if (c < facet_coord+distance) {
      point_coord[d] = facet_coord+distance;
      flag_moved = true;
    }
  }
  else {
    if (c > facet_coord+(1-distance)) {
      point_coord[d] = facet_coord+(1-distance);
      flag_moved = true;
    }
  }
}


// Move interval volume vertex so it is at least distance 
//   from grid cube facet jfacet.
void IVOLDUAL::move_ivol_vertex_away_from_cube_facet
(const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const IVOL_VERTEX_INDEX ivolv, const int jfacet, const COORD_TYPE distance, 
 COORD_ARRAY & vertex_coord, bool & flag_moved)
{
  const int DIM3(3);
  const GRID_COORD_TYPE * cube_coord =
    IJKDUAL::get_coord_of_cube_containing_isov(cube_list, ivolv_list, ivolv);

  move_point_away_from_cube_facet_3D
    (cube_coord, jfacet, distance, &(vertex_coord[DIM3*ivolv]), flag_moved);
}


// Move point so it is at least distance from all cube facets.
// @pre 0 <= distance <= 0.5.
void IVOLDUAL::move_point_away_from_all_cube_facets_3D
(const GRID_COORD_TYPE cube_coord[], const COORD_TYPE distance, 
 COORD_TYPE point_coord[], bool & flag_moved)
{
  const int DIM3(3);

  flag_moved = false;
  for (int d = 0; d < DIM3; d++) {
    GRID_COORD_TYPE facet_coord = cube_coord[d];
    COORD_TYPE c = point_coord[d];
    if (c < facet_coord+distance) {
      point_coord[d] = facet_coord+distance;
      flag_moved = true;
    }
    else if (c > facet_coord+(1-distance)) {
      point_coord[d] = facet_coord+(1-distance);
      flag_moved = true;
    }
  }
}


// Move interval volume vertex so it is at least distance 
//   from all grid cube facets.
void IVOLDUAL::move_ivol_vertex_away_from_all_cube_facets
(const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const IVOL_VERTEX_INDEX ivolv, const COORD_TYPE distance, 
 COORD_ARRAY & vertex_coord, bool & flag_moved)
{
  const int DIM3(3);
  const GRID_COORD_TYPE * cube_coord =
    IJKDUAL::get_coord_of_cube_containing_isov(cube_list, ivolv_list, ivolv);

  move_point_away_from_all_cube_facets_3D
    (cube_coord, distance, &(vertex_coord[DIM3*ivolv]), flag_moved);
}


// *****************************************************************
// MOVE POINT/VERTEX IN BOUNDARY CUBE AWAY FROM GRID INTERIOR
// *****************************************************************

// Move point away from cubes in grid interior.
// @pre 0 <= distance <= 1.
void IVOLDUAL::move_point_away_from_grid_interior_3D
(const DUALISO_GRID & grid, 
 const GRID_COORD_TYPE cube_coord[], const COORD_TYPE distance, 
 COORD_TYPE point_coord[], bool & flag_moved)
{
  const int DIM3(3);

  for (int d = 0; d < DIM3; d++) {

    GRID_COORD_TYPE facet_coord = cube_coord[d];
    COORD_TYPE c = point_coord[d];
    if (facet_coord == 0 && c > (1-distance)) {
      point_coord[d] = 1-distance;
      flag_moved = true;
    }
    else if (facet_coord+2 == grid.AxisSize(d) &&
             c < facet_coord+distance) {
      point_coord[d] = facet_coord+distance;
      flag_moved = true;
    }
  }
}


// Move interval volume vertex away from cubes in grid interior.
// @pre Vertex in in grid boundary cube.
void IVOLDUAL::move_ivol_vertex_away_from_grid_interior
(const DUALISO_GRID & grid, 
 const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const IVOL_VERTEX_INDEX ivolv, const COORD_TYPE distance, 
 COORD_ARRAY & vertex_coord, bool & flag_moved)
{
  const int DIM3(3);
  const GRID_COORD_TYPE * cube_coord =
    IJKDUAL::get_coord_of_cube_containing_isov(cube_list, ivolv_list, ivolv);

  move_point_away_from_grid_interior_3D
    (grid, cube_coord, distance, &(vertex_coord[ivolv*DIM3]), flag_moved);
}


// *****************************************************************
// EXPAND THIN REGIONS
// *****************************************************************

// Expand thin regions.
// Move isosurface vertices in thin regions away from grid cube facets.
void IVOLDUAL::expand_thin_regions
(const DUALISO_GRID & grid,
 const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const COORD_TYPE separation_distance,
 COORD_ARRAY & vertex_coord,
 int & num_moved)
{
  typedef DUAL_IVOLVERT::DIR_BITS_TYPE DIR_BITS_TYPE;

  const int DIM3(3);
  DIR_BITS_TYPE mask[DIM3] = 
    { DIR_BITS_TYPE(1), DIR_BITS_TYPE(2), DIR_BITS_TYPE(4) };
  bool flag_moved1, flag_moved2;

  num_moved = 0;

  for (IVOL_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) {

    if (!ivolv_list[ivolv].in_thin_region) { continue; }
    if (ivolv_list[ivolv].NumIncidentIsoQuad() != 4) { continue; }
    if (ivolv_list[ivolv].separation_vertex >= grid.NumVertices())
      { continue; }

    const int edge_dir = ivolv_list[ivolv].separation_edge_direction;

    const int d1 = (edge_dir+1)%DIM3;
    const int d2 = (edge_dir+2)%DIM3;
    const DIR_BITS_TYPE connect_dir = ivolv_list[ivolv].connect_dir;

    if ((mask[d1] & connect_dir) == 0) {
      move_ivol_vertex_away_from_cube_facet
        (cube_list, ivolv_list, ivolv, d1+DIM3, separation_distance, 
         vertex_coord, flag_moved1);
    }
    else {
      move_ivol_vertex_away_from_cube_facet
        (cube_list, ivolv_list, ivolv, d1, separation_distance, 
         vertex_coord, flag_moved1);
    }

    if ((mask[d2] & connect_dir) == 0) {
      move_ivol_vertex_away_from_cube_facet
        (cube_list, ivolv_list, ivolv, d2+DIM3, separation_distance, 
         vertex_coord, flag_moved2);
    }
    else {
      move_ivol_vertex_away_from_cube_facet
        (cube_list, ivolv_list, ivolv, d2, separation_distance, 
         vertex_coord, flag_moved2);
    }

    if (flag_moved1 || flag_moved2) { num_moved++; }
  }

}
