/// @file ivoldual_compute.h
/// @brief Compute routines for ivoldual.
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

#ifndef _IVOLDUAL_COMPUTE_
#define _IVOLDUAL_COMPUTE_

#include "ivoldual_types.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

  /// Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// - Version with C++ STL vector hex_vert[] containing 
  ///     an array of vertices of multiple hexahedra.
  /// @param ihex Hexahedron index.
  void compute_min_max_hexahedron_Jacobian_determinant
  (const std::vector<VERTEX_INDEX> & hex_vert,
   const int ihex,
   const std::vector<COORD_TYPE> & vertex_coord,
   COORD_TYPE & min_Jacobian_determinant,
   COORD_TYPE & max_Jacobian_determinant);

  /// Compute the Jacobian matrix determinants of a hexahedron
  ///   at a given corner.
  void compute_hexahedron_Jacobian_determinant
  (const std::vector<VERTEX_INDEX> & hex_vert,
   const int ihex,
   const std::vector<COORD_TYPE> & vertex_coord,
   const int icorner,
   COORD_TYPE & Jacobian_determinant);

  /// Compute the normalized Jacobian matrix determinants of 
  /// a hexahedron at a given corner.
  void compute_hexahedron_normalized_Jacobian_determinant
  (const std::vector<VERTEX_INDEX> & hex_vert,
   const int ihex,
   const std::vector<COORD_TYPE> & vertex_coord,
   const int icorner,
   COORD_TYPE & Jacobian_determinant);

}

#endif
