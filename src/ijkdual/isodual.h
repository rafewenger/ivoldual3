/// @file isodual.h
/// @brief Generate isosurface in arbitrary dimensions
/// - Version 0.3.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2026 Rephael Wenger

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

/*!
  \mainpage ISODUAL: BASIC DUAL CONTOURING

  ISODUAL is a program for generating isosurfaces 
  using the dual contouring algorithm.  It computes
  isosurfaces from two, three or four dimensional volumetric grid data.
*/

#ifndef _ISODUAL_
#define _ISODUAL_

#include <string>

#include "ijk.tpp"

#include "isodual_types.h"
#include "isodual_datastruct.h"

/// isodual classes and routines.
namespace ISODUAL {

// **************************************************
// DUAL CONTOURING 
// **************************************************

  /// Dual Contouring Algorithm.
  void dual_contouring
    (const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue, 
     DUAL_ISOSURFACE & dual_isosurface, DUALISO_INFO & dualiso_info);


// **************************************************
// DUAL CONTOURING
// **************************************************

  /// Extract isosurface using Manifold Dual Contouring algorithm.
  /// Isosurface must be 2D surface in 3D volume.
  /// Allow multiple isosurface vertices per grid cube.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates.
  /// Isosurface mesh has a manifold representation.
  void dual_contouring_manifold
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const IJKDUAL::ISODUAL_CUBE_TABLE & isodual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   GRID_EDGE_ARRAY & dual_edge,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);

  /// Extract isosurface using Dual Contouring algorithm.
  /// Allow multiple isosurface vertices per grid cube.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates
  /// @param param Algorithm control flags/parameters.
  /// @param dual_edge[] List of grid edges dual to isosurface polytopes.
  /// @param merge_data = Data structure for merging edges.  
  void dual_contouring_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   GRID_EDGE_ARRAY & dual_edge,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);

  /// Extract isosurface using Dual Contouring algorithm.
  /// Allow multiple isosurface vertices per grid cube.
  /// Version which creates isodual table or qdual table.
  void dual_contouring_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   GRID_EDGE_ARRAY & dual_edge,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);

  /// Extract isosurface using Dual Contouring algorithm
  /// Single isosurface vertex per grid cube.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates
  void dual_contouring_single_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const VERTEX_POSITION_METHOD vertex_position_method,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   GRID_EDGE_ARRAY & dual_edge,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info);
}

#endif
