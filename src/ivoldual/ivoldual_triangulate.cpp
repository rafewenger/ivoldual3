/// @file ivoldual_triangulate.cpp
/// @brief Triangulate (tetrahedralize) interval volume hexahedra.
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

#include "ijktriangulate.txx"

#include "ijkdual_triangulate.tpp"

#include "ivoldual_triangulate.h"

using namespace IVOLDUAL;


// **************************************************
// Triangulate interval volume hexahedra
// **************************************************

// Triangulate the interval volume.
void IVOLDUAL::triangulate_interval_volume
(const IVOLDUAL_DATA & ivoldual_data, 
 DUAL_INTERVAL_VOLUME & dual_interval_volume)
{
  if (ivoldual_data.flag_add_isov_dual_to_hexahedra) {
    add_vertex_coord_at_hexahedra_centroids
      (ivoldual_data, dual_interval_volume);

    /* NOT PROPERLY WORKING.
    triangulate_ivol_hexahedra_tri12_diagonal07
      (ivoldual_data.ScalarGrid(), dual_interval_volume.isopoly_vert,
       dual_interval_volume.isopoly_info, 
       dual_interval_volume.first_isov_dual_to_iso_poly,
       dual_interval_volume.tri_vert);
    */

    // *** FOR NOW, USE TRIANGULATION INTO 6 TETRAHEDRA ***
    triangulate_ivol_hexahedra_diagonal
      (dual_interval_volume.isopoly_vert,
       dual_interval_volume.tri_vert);
  }
  else {
    triangulate_ivol_hexahedra_diagonal
      (dual_interval_volume.isopoly_vert,
       dual_interval_volume.tri_vert);
  }
}



// Triangulate the interval volume hexahedra using hexahedra diagonals.
void IVOLDUAL::triangulate_ivol_hexahedra_diagonal
(const VERTEX_INDEX_ARRAY & ivolpoly_vert,
 VERTEX_INDEX_ARRAY & tri_vert)
{
  const int NUM_VERT_PER_HEXAHEDRON(8);
  IJK::POLYMESH_DATA<VERTEX_INDEX,int,
    IJK::HEX_TRIANGULATION_INFO<char,char>> hex_data;

  hex_data.AddPolytopes(ivolpoly_vert, NUM_VERT_PER_HEXAHEDRON);
  IJK::triangulate_hex_meshX(hex_data, tri_vert);
}


// Uniformly triangulate the interval volume hexahedra 
//   using an additional vertex inside each hexahedron.
void IVOLDUAL::triangulate_ivol_hexahedra_tri12_diagonal07
(const DUALISO_GRID & grid,
 const VERTEX_INDEX_ARRAY & ivolpoly_vert,
 const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 const VERTEX_INDEX first_isov_dual_to_iso_poly,
 VERTEX_INDEX_ARRAY & tri_vert)
{
  const int NUM_VERT_PER_HEXAHEDRON(8);
  const int num_hex = ivolpoly_vert.size()/NUM_VERT_PER_HEXAHEDRON;

  for (int ihex = 0; ihex < num_hex; ihex++) {

    const VERTEX_INDEX * hex_vert = 
      &(ivolpoly_vert[ihex*NUM_VERT_PER_HEXAHEDRON]);
    const VERTEX_INDEX ivert_in_hex = first_isov_dual_to_iso_poly+ihex;

    triangulate_hexahedron_tri12_diagonal07(hex_vert, ivert_in_hex, tri_vert);
  }

}


// Triangulate hexahedron using an additional vertex w0
// into 12 hexahedra.
// Hex facet triangles should match hex diagonal07 triangulation.
void IVOLDUAL::triangulate_hexahedron_tri12_diagonal07
(const VERTEX_INDEX hex_vert[], const VERTEX_INDEX w0,
 std::vector<VERTEX_INDEX> & tet_vert_list)
{
  using IJK::add_tetrahedron_vertices;

  add_tetrahedron_vertices
    (w0, hex_vert[7], hex_vert[3], hex_vert[2], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[3], hex_vert[2], hex_vert[0], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[7], hex_vert[2], hex_vert[6], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[2], hex_vert[6], hex_vert[0], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[7], hex_vert[6], hex_vert[4], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[6], hex_vert[4], hex_vert[0], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[7], hex_vert[4], hex_vert[5], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[4], hex_vert[5], hex_vert[0], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[7], hex_vert[5], hex_vert[1], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[5], hex_vert[1], hex_vert[0], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[7], hex_vert[1], hex_vert[3], tet_vert_list);
  add_tetrahedron_vertices
    (w0, hex_vert[1], hex_vert[3], hex_vert[0], tet_vert_list);
}
