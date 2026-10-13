/// @file isodual_types.h
/// @brief Type definitions for isodual.
/// - Copies types from ijkdual.h.
/// - Version 0.3.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2026 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _ISODUAL_TYPES_
#define _ISODUAL_TYPES_

#include "ijkdual_types.h"

namespace ISODUAL {

  // **************************************************
  // SCALAR TYPES
  // **************************************************

  typedef IJKDUAL::SCALAR_TYPE SCALAR_TYPE;
  typedef IJKDUAL::COORD_TYPE COORD_TYPE;
  typedef IJKDUAL::VERTEX_INDEX VERTEX_INDEX;
  typedef IJKDUAL::AXIS_SIZE_TYPE AXIS_SIZE_TYPE;
  typedef IJKDUAL::ISO_VERTEX_INDEX ISO_VERTEX_INDEX;


  // **************************************************
  // ARRAY TYPES
  // **************************************************

  typedef IJKDUAL::COORD_ARRAY COORD_ARRAY;
  typedef IJKDUAL::SCALAR_ARRAY SCALAR_ARRAY;


  // **************************************************
  // ENUMERATED TYPES
  // **************************************************

  typedef IJKDUAL::VERTEX_POSITION_METHOD VERTEX_POSITION_METHOD;

}

#endif
