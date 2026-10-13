/// @file ivoldual_triangulate.h
/// @brief Triangulate (tetrahedralize) interval volume hexahedra.
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


#ifndef _IVOLDUAL_TRIANGULATE_
#define _IVOLDUAL_TRIANGULATE_

#include "ivoldual_datastruct.h"

namespace IVOLDUAL {

  // **************************************************
  // Triangulate interval volume hexahedra
  // **************************************************


  /// Triangulate the interval volume.
  void triangulate_interval_volume
  (const IVOLDUAL_DATA & ivoldual_data, 
   DUAL_INTERVAL_VOLUME & dual_interval_volume);

  /// Triangulate the interval volume hexahedra using hexahedra diagonals.
  void triangulate_ivol_hexahedra_diagonal
  (const VERTEX_INDEX_ARRAY & ivolpoly_vert,
   VERTEX_INDEX_ARRAY & tri_vert);


  // **************************************************
  // DEPRECATED. CREATE NON-CONFORMING MESHES.
  // **************************************************

  /// Uniformly triangulate the interval volume hexahedra 
  ///   using an additional vertex inside each hexahedron.
  void triangulate_ivol_hexahedra_tri12_diagonal07
  (const DUALISO_GRID & grid,
   const VERTEX_INDEX_ARRAY & ivolpoly_vert,
   const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   const VERTEX_INDEX first_isov_dual_to_iso_poly,
   VERTEX_INDEX_ARRAY & tri_vert);

  /// Triangulate hexahedron using an additional vertex w0
  /// into 12 hexahedra.
  /// Hex facet triangules should match hex diagonal07 triangulation.
  void triangulate_hexahedron_tri12_diagonal07
  (const VERTEX_INDEX hex_vert[], const VERTEX_INDEX w0,
   std::vector<VERTEX_INDEX> & tet_vert_list);

}

#endif
