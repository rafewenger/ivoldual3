/// @file ijkdual_types.h
/// @brief Type definitions for ijkdual.
/// - Version 0.3.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2026 Rephael Wenger

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

#ifndef _IJKDUAL_TYPES_
#define _IJKDUAL_TYPES_

#include <string>

#include "ijk.tpp"
#include "ijkscalar_grid.tpp"
#include "ijkmerge.tpp"

namespace IJKDUAL {


// **************************************************
// SCALAR TYPES
// **************************************************

  typedef float SCALAR_TYPE;     ///< Scalar value type.
  typedef float COORD_TYPE;      ///< Isosurface vertex coordinate type.
  typedef int VERTEX_INDEX;      ///< Grid vertex index type.
  typedef float GRADIENT_TYPE;   ///< Gradient coordinate type.
  typedef int GRID_COORD_TYPE;   ///< Grid vertex coordinate type.
  typedef int AXIS_SIZE_TYPE;    ///< Axis size type.
  typedef int ISO_VERTEX_INDEX;  ///< Isosurface vertex index type.
  typedef int MERGE_INDEX;       ///< Merge index type.

  // *** SHOULD CHANGE TO unsigned char ***
  typedef int DIRECTION_TYPE;    ///< Direction type.

  /// Boundary bits type.
  /// Number of bits must be at least 2*dimension.
  typedef unsigned int BOUNDARY_BITS_TYPE;

  /// Edge index type.
  /// Vertex and edge indices must have the same type.
  typedef VERTEX_INDEX EDGE_INDEX;

  /// Facet index type.
  typedef unsigned char FACET_INDEX;

  /// Facet vertex index type.
  typedef unsigned char FACET_VERTEX_INDEX;

  /// Type for integer representing a set of facets.
  /// - k'th bit is 1 if facet k is in the set.
  typedef int FACET_BITS_TYPE;

  /// Vertex degree type.
  typedef unsigned char DEGREE_TYPE;

  /// Isosurface lookup table index type.
  typedef unsigned long TABLE_INDEX;


  // **************************************************
  // ARRAY TYPES
  // **************************************************

  typedef std::vector<COORD_TYPE> COORD_ARRAY;   ///< Grid coordinate array.
  typedef std::vector<VERTEX_INDEX>              /// Vertex index array.
    VERTEX_INDEX_ARRAY; 
  typedef std::vector<SCALAR_TYPE> SCALAR_ARRAY; ///< Scalar array.

  /// Array of facet vertex indices.
  typedef std::vector<FACET_VERTEX_INDEX> FACET_VERTEX_INDEX_ARRAY;

  /// Array of directions.
  typedef std::vector<DIRECTION_TYPE> DIRECTION_ARRAY;


  // **************************************************
  // ENUMERATED TYPES
  // **************************************************

  /// Interpolation type.
  /// - LINEAR_INTERPOLATION: Determine the location of an isosurface vertex
  ///   using linear interpolation on the endpoints of the cube, pyramid
  ///   or simplex containing the isosurface vertex.
  /// - MULTILINEAR_INTERPOLATION: Determine the location of an 
  ///   isosurface vertex using multilinear interpolation on the cube vertices.
  typedef enum { LINEAR_INTERPOLATION, MULTILINEAR_INTERPOLATION }
  INTERPOLATION_TYPE;

  /// Isosurface vertex position method.
  /// CUBE_CENTER: Position isosurface vertices at cube centers.
  /// CENTROID_EDGE_ISO: Position isosurface vertices at the centroid
  ///                    of the edge isosurface intersections.
  /// DIAGONAL_INTERPOLATION: Position some isosurface vertices on diagonals
  ///   using linear interpolation.  Position all other isosurface vertices
  ///   using at the centroid of the edge isosurface intersections.
  typedef enum { CUBE_CENTER, CENTROID_EDGE_ISO, DIAGONAL_INTERPOLATION,
                 IVOL_LIFTED02, QDUAL_INTERPOLATION } 
  VERTEX_POSITION_METHOD;

  /// Quadrilateral triangulation method.
  typedef enum 
  { UNDEFINED_TRI, UNIFORM_TRI, SPLIT_MAX_ANGLE, MAX_MIN_ANGLE,
    ONLY_TRI4, TRI4_BY_DISTANCE, TRI4_MAX_MIN_ANGLE }
  QUAD_TRI_METHOD;

  /// Position method for additional vertex in quadrilateral triangulations.
  typedef enum
  { TRI4_ON_GRID_EDGE, TRI4_CENTROID } TRI4_POSITION_METHOD;

  /// Quad-edge intersection method.
  typedef enum
    { QEI_INTERPOLATE_COORD, QEI_INTERPOLATE_SCALAR, 
      QEI_AVERAGE, QEI_WEIGHTED_AVERAGE }
  QUAD_EDGE_INTERSECTION_METHOD;

  typedef enum { TRIANGLE, QUADRILATERAL, DEGENERATE_POLYGON }
    POLYGON_TYPE;


// **************************************************
// CLASSES
// **************************************************

  typedef IJK::BOX<VERTEX_INDEX> GRID_BOX;  ///< Grid box type.
}

#endif
