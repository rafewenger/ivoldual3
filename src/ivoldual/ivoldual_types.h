/*!
 *  @file ivoldual_types.h
 *  @brief Type definitions for ivoldual.
 *  - Copies types from ijkdual.h.
 *  - Version 0.3.1
 */

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2026 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IVOLDUAL_TYPES_
#define _IVOLDUAL_TYPES_

#include "ijkdual_types.h"

namespace IVOLDUAL {

  // **************************************************
  // SCALAR TYPES
  // **************************************************

  typedef IJKDUAL::SCALAR_TYPE SCALAR_TYPE;
  typedef IJKDUAL::COORD_TYPE COORD_TYPE;
  typedef IJKDUAL::GRID_COORD_TYPE GRID_COORD_TYPE;
  typedef IJKDUAL::VERTEX_INDEX VERTEX_INDEX;
  typedef IJKDUAL::AXIS_SIZE_TYPE AXIS_SIZE_TYPE;
  typedef IJKDUAL::ISO_VERTEX_INDEX ISO_VERTEX_INDEX;
  typedef IJKDUAL::ISO_VERTEX_INDEX IVOL_VERTEX_INDEX;
  typedef IJKDUAL::TABLE_INDEX TABLE_INDEX;
  typedef IJKDUAL::FACET_INDEX FACET_INDEX;
  typedef IJKDUAL::FACET_VERTEX_INDEX FACET_VERTEX_INDEX;
  typedef IJKDUAL::FACET_BITS_TYPE FACET_BITS_TYPE;
  typedef IJKDUAL::DEGREE_TYPE DEGREE_TYPE;

  typedef unsigned char POLY_VERTEX_INDEX;

  /// Boundary bits type.
  /// Number of bits must be at least 2*dimension.
  typedef IJKDUAL::BOUNDARY_BITS_TYPE BOUNDARY_BITS_TYPE;


  // **************************************************
  // ARRAY TYPES
  // **************************************************

  typedef IJKDUAL::COORD_ARRAY COORD_ARRAY;
  typedef IJKDUAL::SCALAR_ARRAY SCALAR_ARRAY;
  typedef IJKDUAL::VERTEX_INDEX_ARRAY VERTEX_INDEX_ARRAY;


  // **************************************************
  // ENUMERATED TYPES
  // **************************************************

  /// Hexahedra triangulation method.
  typedef enum
    { UNDEFINED_HEX_TRI, UNIFORM_HEX_TRI_DIAGONAL07, ONLY_TRI12,
    TRI12_MAX_MIN_ANGLE }
    HEX_TRI_METHOD;

  typedef IJKDUAL::VERTEX_POSITION_METHOD VERTEX_POSITION_METHOD;

}

#endif
