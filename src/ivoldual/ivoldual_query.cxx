/// \file ivoldual_query.cxx
/// Query routines for ivoldual.


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017 Rephael Wenger

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


#include "ivoldual_query.h"


// *****************************************************************
// QUERY IN BOX OR PSEUDOBOX
// *****************************************************************

// Return true if vertex is lower left vertex of an isosurface box.
bool IVOLDUAL::is_ivolv_lower_left_vertex_of_isosurface_box
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const IVOL_VERTEX_INDEX ivolv0,
 IVOL_VERTEX_INDEX box_corner[8])
{
  const int DIM3(3);
  const int NUM_CUBE_VERTICES(8);
  const int NUM_CUBE_FACET_VERTICES(NUM_CUBE_VERTICES/2);
  const VERTEX_INDEX separation_vertex = 
    ivolv_list[ivolv0].SeparationVertex();

  if (ivolv_list[ivolv0].NumIncidentIsoQuad() != 3) { return(false); }

  box_corner[0] = ivolv0;
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (ivolv0, 0, 1, box_corner[1])) { return(false); }
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (ivolv0, 1, 1, box_corner[2])) { return(false); }
  if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
      (box_corner[1], 1, 1, box_corner[3])) { return(false); }

  for (int j = 0; j < NUM_CUBE_FACET_VERTICES; j++) {
    if (!vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
        (box_corner[j], 2, 1, box_corner[4+j])) { return(false); }
  }

  for (int i = 1; i < NUM_CUBE_VERTICES; i++) {
    const VERTEX_INDEX ivolv = box_corner[i];
    if (ivolv_list[ivolv].separation_vertex != separation_vertex)
      { return(false); }
    if (ivolv_list[ivolv].NumIncidentIsoQuad() != 3)
      { return(false); }
  }

  return(true);
}


// Return true if upper, rightmost vertex of cube is surrounded
//   by an isosurface pseudobox.
bool IVOLDUAL::is_upper_right_grid_vertex_in_isosurface_pseudobox
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const VERTEX_INDEX cube0_list_index,
 VERTEX_INDEX cube_list_index[8])
{
  typedef IVOLDUAL_TABLE_VERTEX_INFO::CUBE_VERTEX_TYPE CUBE_VERTEX_TYPE;

  const int DIM3(3);
  const int NUM_CUBE_VERTICES(8);
  const int NUM_CUBE_FACET_VERTICES(NUM_CUBE_VERTICES/2);
  const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
    IVOLDUAL_TABLE_VERTEX_INFO::UNDEFINED_CUBE_VERTEX;
  static const CUBE_FACE_INFO cube(DIM3);

  if (!get_seven_adjacent_cubes
      (vertex_adjacency_list, cube_list, ivolv_list, cube0_list_index,
       cube_list_index))
    { return(false); }

  // The surrounded vertex 
  const VERTEX_INDEX iv0 = cube_list[cube_list_index[7]].cube_index;

  // *** NEED TO ALSO CHECK THAT ALL GRID EDGES INCIDENT ON iv0
  //     ARE INTERSECTED BY THE APPROPRIATE ISOSURFACE. ***

  for (int  d = 0; d < DIM3; d++) {

    // flag_connects is true if parallel quads with orth dir d
    //   are connected by some interval volume edge.
    bool flag_connects = false;
    for (int k = 0; k < NUM_CUBE_FACET_VERTICES; k++) {

      CUBE_VERTEX_TYPE icorner0 = cube.FacetVertex(d, k);
      CUBE_VERTEX_TYPE icorner1 = cube.VertexNeighbor(icorner0, d);
      VERTEX_INDEX i0 = cube_list_index[icorner0];
      VERTEX_INDEX i1 = cube_list_index[icorner1];

      bool flag0 = 
        (is_vertex_sep_vertex(iv0, cube_list[i0], ivolv_list) ||
         does_cube_contain_doubly_connected_ivol_vertex
         (cube_list[i0], ivolv_list));
      bool flag1 = 
        (is_vertex_sep_vertex(iv0, cube_list[i1], ivolv_list) ||
         does_cube_contain_doubly_connected_ivol_vertex
         (cube_list[i1], ivolv_list));

      if (flag0 && flag1) {
        flag_connects = true;
        break;
      }

    }

    if (!flag_connects) { return(false); }
  }

  return(true);
}


/// Return true if vertex iv equals w.sep_vert for some interval volume 
///   vertex w in cube cube_ivolv.
bool IVOLDUAL::is_vertex_sep_vertex
(const VERTEX_INDEX iv,
 const GRID_CUBE_DATA & cube_ivolv,
 const DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  for (int i = 0; i < cube_ivolv.num_isov; i++) {
    int k = cube_ivolv.first_isov + i;
    if ((ivolv_list[k].NumIncidentIsoQuad() == 3) &&
        (ivolv_list[k].separation_vertex == iv))
      { return(true); }
  }

  return(false);
}


// Return true if cube contains doubly connected interval volume vertex.
bool IVOLDUAL::does_cube_contain_doubly_connected_ivol_vertex
(const GRID_CUBE_DATA & cube_ivolv,
 const DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  for (int i = 0; i < cube_ivolv.num_isov; i++) {
    int k = cube_ivolv.first_isov + i;
    if (ivolv_list[k].is_doubly_connected) { return(true); }
  }

  return(false);
}


// *****************************************************************
// GET FUNCTIONS
// *****************************************************************

// Get the index in cube_list of the seven cubes sharing the
//   rightmost/upper vertex with cube0_list_index.
// - Gets the cubes only if there are interval volume edges
//   between interval volume vertices in each pair of cubes.
// - Return true if all 7 cubes found.
// @param[out] cube_list_index[] Cubes cube0_list_index
//     and the seven other cubes adjacent to cube0_list_index.
bool IVOLDUAL::get_seven_adjacent_cubes
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const int cube0_list_index,
 VERTEX_INDEX cube_list_index[8])
{
  const int NUM_CUBE_VERTICES(8);
  const int NUM_CUBE_FACET_VERTICES(NUM_CUBE_VERTICES/2);

  int list_index;

  cube_list_index[0] = cube0_list_index;

  if (!get_adjacent_cube_in_oriented_direction
      (vertex_adjacency_list, cube_list, ivolv_list,
       cube0_list_index, 0, 1, cube_list_index[1]))
    { return(false); }
  if (!get_adjacent_cube_in_oriented_direction
      (vertex_adjacency_list, cube_list, ivolv_list,
       cube0_list_index, 1, 1, cube_list_index[2]))
    { return(false); }
  if (!get_adjacent_cube_in_oriented_direction
      (vertex_adjacency_list, cube_list, ivolv_list,
       cube_list_index[1], 1, 1, cube_list_index[3]))
    { return(false); }

  for (int j = 0; j < NUM_CUBE_FACET_VERTICES; j++) {
    if (!get_adjacent_cube_in_oriented_direction
      (vertex_adjacency_list, cube_list, ivolv_list,
       cube_list_index[j], 2, 1, cube_list_index[4+j]))
      { return(false); }
  }

  return(true);
}


/// Get the index in cube_list of the adjacent cube in the oriented direction
///  if it shares an interval volume edge with cube0_list_index.
bool IVOLDUAL::get_adjacent_cube_in_oriented_direction
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const std::vector<GRID_CUBE_DATA> & cube_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 const int cubeA_list_index,
 const int direction,
 const int orientation,
 int & cubeB_list_index)
{
  ISO_VERTEX_INDEX ivolvB;

  for (int i = 0; i < cube_list[cubeA_list_index].num_isov; i++) {

    const ISO_VERTEX_INDEX ivolvA = 
      cube_list[cubeA_list_index].first_isov + i;

    if (vertex_adjacency_list.GetAdjacentVertexInOrientedDirection
        (ivolvA, direction, orientation, ivolvB)) {

      cubeB_list_index = ivolv_list[ivolvB].cube_list_index;
      return(true);
    }
  }

  return(false);
}

