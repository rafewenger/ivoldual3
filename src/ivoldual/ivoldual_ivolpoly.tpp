/// @file ivoldual_ivolpoly.tpp
/// @brief Templates for extracting interval volume patches from a polyhedron.
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

#ifndef _IVOLDUAL_IVOLPOLY_
#define _IVOLDUAL_IVOLPOLY_

namespace IVOLDUAL {

  // **************************************************
  // SUBROUTINES FOR EXTRACTING INTERVAL VOLUME PATCHES
  // **************************************************

  /// Compute interval volume table index of polyhedron 
  ///   with primary vertex iv0.
  template <const int NUM_VERTEX_TYPES,
            typename STYPE, typename LOWER_STYPE, typename UPPER_STYPE,
            typename VTYPE, typename INC_TYPE, typename NTYPE, 
            typename CTYPE0, typename CTYPE1,
            typename TABLE_INDEX_TYPE>
  inline void compute_ivoltable_index
  (const STYPE * scalar, 
   const LOWER_STYPE lower_isovalue, const UPPER_STYPE upper_isovalue,
   const VTYPE iv0, const INC_TYPE * increment,
   const NTYPE num_poly_vertices, 
   const CTYPE0 interior_code, const CTYPE1 above_code,
   TABLE_INDEX_TYPE & table_index)
  {
    table_index = 0;
    TABLE_INDEX_TYPE x = 1;
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
      VTYPE iv1 = iv0 + increment[j];
      STYPE s = scalar[iv1];
      if (s > upper_isovalue) {
        // Vertex iv1 is above the interval volume.
        table_index = table_index + above_code*x;
      }
      else if (s >= lower_isovalue) {
        // Vertex iv1 is in the interval volume.
        table_index = table_index + interior_code*x;
      }
      x = x*NUM_VERTEX_TYPES;
    }
  }

  /// Compute interval volume table index of polyhedron 
  ///   with primary vertex iv0.
  template <const int NUM_VERTEX_TYPES,
            typename SGRID_TYPE, typename LOWER_STYPE, typename UPPER_STYPE,
            typename ICUBE_TYPE, typename CTYPE0, typename CTYPE1,
            typename TABLE_INDEX_TYPE>
  inline void compute_ivoltable_index_of_grid_cube
  (const SGRID_TYPE & scalar_grid,
   const LOWER_STYPE lower_isovalue, const UPPER_STYPE upper_isovalue,
   const ICUBE_TYPE icube, const CTYPE0 interior_code, const CTYPE1 above_code,
   TABLE_INDEX_TYPE & table_index)
  {
    compute_ivoltable_index<NUM_VERTEX_TYPES>
      (scalar_grid.ScalarPtrConst(), lower_isovalue, upper_isovalue, icube, 
       scalar_grid.CubeVertexIncrement(), scalar_grid.NumCubeVertices(), 
       interior_code, above_code, table_index);
  }

}

#endif

