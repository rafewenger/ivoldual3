/// @file isodual_datastruct.h
/// @brief Data structure definitions for isodual.
/// - Copies types from ijkdual.h.
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

#ifndef _ISODUAL_DATASTRUCT_
#define _ISODUAL_DATASTRUCT_

#include "isodual_types.h"
#include "ijkdual_datastruct.h"


namespace ISODUAL {

  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef IJKDUAL::DUALISO_GRID DUALISO_GRID;
  typedef IJKDUAL::DUALISO_SCALAR_GRID_BASE DUALISO_SCALAR_GRID_BASE;
  typedef IJKDUAL::DUALISO_SCALAR_GRID DUALISO_SCALAR_GRID;


  // **************************************************
  // DUAL ISOSURFACE VERTICES
  // **************************************************

  typedef IJKDUAL::DUAL_ISOVERT DUAL_ISOVERT;


  // **************************************************
  // DUAL CONTOURING ISOSURFACE CLASS
  // **************************************************

  typedef IJKDUAL::DUAL_ISOSURFACE DUAL_ISOSURFACE;


  // **************************************************
  // DUAL CONTOURING INPUT DATA AND DATA STRUCTURES
  // **************************************************

  typedef IJKDUAL::DUALISO_DATA_FLAGS DUALISO_DATA_FLAGS;
  typedef IJKDUAL::DUALISO_DATA DUALISO_DATA;


  // **************************************************
  // DUALISO TIME
  // **************************************************

  typedef IJKDUAL::DUALISO_TIME DUALISO_TIME;


  // **************************************************
  // DUALISO INFO
  // **************************************************

  typedef IJKDUAL::DUALISO_INFO DUALISO_INFO;


  // **************************************************
  // GRID EDGE
  // **************************************************

  typedef IJKDUAL::GRID_EDGE_ARRAY GRID_EDGE_ARRAY;


  // **************************************************
  // MERGE DATA
  // **************************************************

  typedef IJKDUAL::MERGE_DATA MERGE_DATA;

}

#endif
