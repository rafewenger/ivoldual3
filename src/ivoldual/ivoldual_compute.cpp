/// @file ivoldual_compute.cpp
/// @brief Compute routines for ivoldual.
/// - Version 0.3.1


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


#include "ijkcoord.tpp"
#include "ijkcube.tpp"

#include "ivoldual_compute.h"

// Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
void IVOLDUAL::compute_min_max_hexahedron_Jacobian_determinant
(const std::vector<VERTEX_INDEX> & hex_vert,
 const int ihex,
 const std::vector<COORD_TYPE> & vertex_coord,
 COORD_TYPE & min_Jacobian_determinant,
 COORD_TYPE & max_Jacobian_determinant)
{
  const int DIM3(3);
  const int POSITIVE_ORIENTATION(1);
  IJK::CUBE_FACE_INFO<int,int,int> cube(DIM3);

  IJK::compute_min_max_hexahedron_Jacobian_determinant_3D
    (hex_vert, ihex, POSITIVE_ORIENTATION, vertex_coord, cube,
     min_Jacobian_determinant, max_Jacobian_determinant);
}


// Compute the Jacobian matrix determinant of a hexahedron
//   at a given corner.
void IVOLDUAL::compute_hexahedron_Jacobian_determinant
(const std::vector<VERTEX_INDEX> & hex_vert,
 const int ihex,
 const std::vector<COORD_TYPE> & vertex_coord,
 const int icorner,
 COORD_TYPE & Jacobian_determinant)
{
  const int DIM3(3);
  const int POSITIVE_ORIENTATION(1);
  IJK::CUBE_FACE_INFO<int,int,int> cube(DIM3);

  IJK::compute_Jacobian_determinant_at_hex_vertex_3D
    (hex_vert, ihex, POSITIVE_ORIENTATION, vertex_coord, cube, icorner,
     Jacobian_determinant);
}

// Compute the normalized Jacobian matrix determinant of a hexahedron
//   at a given corner.
void IVOLDUAL::compute_hexahedron_normalized_Jacobian_determinant
(const std::vector<VERTEX_INDEX> & hex_vert,
 const int ihex,
 const std::vector<COORD_TYPE> & vertex_coord,
 const int icorner,
 COORD_TYPE & Jacobian_determinant)
{
  const int DIM3(3);
  const int POSITIVE_ORIENTATION(1);
  const int NUM_VERT_PER_HEX(8);
  const VERTEX_INDEX * hex_i_vert = &(hex_vert[ihex*NUM_VERT_PER_HEX]);
  const COORD_TYPE * vcoord = IJK::vector2pointer(vertex_coord);

  COORD_TYPE max_small_magnitude(0.0);
  bool flag_zero;
  IJK::CUBE_FACE_INFO<int,int,int> cube(DIM3);

  IJK::compute_normalized_Jacobian_determinant_at_hex_vertex_3D
    (hex_i_vert, POSITIVE_ORIENTATION, vcoord, cube, icorner,
      max_small_magnitude, Jacobian_determinant, flag_zero);
}

