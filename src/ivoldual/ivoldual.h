/// @file ivoldual.h
/// @brief Construct interval volume in arbitrary dimensions using dual contouring.
/// - Version 0.3.0

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

/*!
  \mainpage IVOLDUAL: EXTRACT INTERVAL VOLUME

  IVOLDUAL is a program for generating interval volumes
  using the dual contouring.
*/

#ifndef _IVOLDUAL_
#define _IVOLDUAL_

#include <string>

#include "ijk.tpp"

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"
#include "ivoldualtable.h"

/// ivoldual classes and routines.
namespace IVOLDUAL {

// **************************************************
// DUAL CONTOURING INTERVAL VOLUME
// **************************************************

  /// Construct interval volume using dual contouring.
  void dual_contouring_interval_volume
  (const IVOLDUAL_DATA & ivoldual_data, 
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   DUAL_INTERVAL_VOLUME & dual_interval_volume, IVOLDUAL_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// - Returns list of interval volume polytope vertices
  ///   and list of interval volume vertex coordinates.
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IVOLDUAL_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   IVOLDUAL_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// - Returns list of interval volume polytope vertices
  ///   and list of interval volume vertex coordinates.
  /// - Version which creates ivoldual_table and ivolv_list.
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const IVOLDUAL_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   IVOLDUAL_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// - Returns list of interval volume polytope vertices
  ///   and list of interval volume vertex coordinates.
  /// - Version which creates ivoldual_table.
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const IVOLDUAL_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   IVOLDUAL_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// - Returns list of interval volume polytope vertices
  ///   and list of interval volume vertex coordinates.
  /// - Version which creates cube_ivolv_list.
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IVOLDUAL_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   IVOLDUAL_INFO & dualiso_info);

  /// Construct interval volume using dual contouring.
  /// - Returns list of interval volume polytope vertices
  ///   and list of interval volume vertex coordinates.
  void dual_contouring_interval_volume
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IVOLDUAL_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
   std::vector<GRID_CUBE_DATA> & cube_ivolv_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   IVOLDUAL_INFO & dualiso_info);


  // **************************************************
  // ENCODE GRID VERTICES
  // **************************************************

  /// Encode grid vertices as 0, 1, 2 or 3.
  /// - 0: Below lower isovalue.
  /// - 1: Between lower and upper isovalue.
  /// - 2: Between lower and upper isovalue.
  /// - 3: Above upper isovalue.
  /// @param default_interior_code Default code (1 or 2) for vertices
  ///   between lower and upper isovalue.
  void encode_grid_vertices
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   const GRID_VERTEX_ENCODING default_interior_code,
   IVOLDUAL_ENCODED_GRID & encoded_grid,
   IVOLDUAL_INFO & dualiso_info);

  /// Encode grid vertices. Set interior codes based on scalar grid.
  /// - Vertices with scalar value between isovalue0 
  ///     and (isovalue0+isovalue1)/2 receive encoded value 1.
  /// - Vertices with scalar value between (isovalue0+isovalue1)/2 
  ///     and isovalue1 receive encoded value 2.
  /// - Used for case analysis.
  void encode_grid_vertices_set_interior_from_scalar
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
   IVOLDUAL_ENCODED_GRID & encoded_grid,
   IVOLDUAL_INFO & dualiso_info);


  // **************************************************
  // SET IVOLTABLE INFO FOR EACH ACTIVE GRID CUBE
  // **************************************************

  /// Set ivoltable information for each cube in cube_ivolv_list.
  /// @param encoded_grid Encoded grid.  
  ///    Each grid vertex has value 0,1,2 or 3, indicating vertex type.
  /// @param ivoldual_table Interval volume lookup table.
  /// @param[out] cube_ivolv_list Array of active grid cubes.
  void set_cube_ivoltable_info
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   std::vector<GRID_CUBE_DATA> & cube_ivolv_list);


  // **************************************************
  // SET IVOL VERTEX INFORMATION
  // **************************************************

  void set_ivol_vertex_info
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const std::vector<ISO_VERTEX_INDEX> & poly_vert,
   const std::vector<GRID_CUBE_DATA> & cube_list,
   const VERTEX_INDEX index_to_cube_list[],
   const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list);

  void determine_ivol_vertices_missing_incident_hex
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const std::vector<ISO_VERTEX_INDEX> & poly_vert,
   DUAL_IVOLVERT_ARRAY & ivolv_list);

  void determine_ivol_vertices_in_isosurface_boxes
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list);

  void determine_ivol_vertices_in_isosurface_pseudoboxes
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const std::vector<GRID_CUBE_DATA> & cube_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list);

  void determine_ivol_vertices_adjacent_to_upper_and_lower_isosurfaces
  (const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list);

  void determine_thin_regions
  (const DUALISO_GRID & grid,
   const std::vector<GRID_CUBE_DATA> & cube_list,
   const VERTEX_INDEX index_to_cube_list[],
   const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list);

  /// Set flag in_pseudobox to true for all ivol vertices in cube cube_ivolv.
  void set_in_pseudobox
  (const GRID_CUBE_DATA & cube_ivolv,
   DUAL_IVOLVERT_ARRAY & ivolv_list);


  // **************************************************
  // EXTRACT DUAL INTERVAL VOLUME POLYTOPES
  // **************************************************

  /// Extract dual interval volume polytopes.
  /// @param encoded_grid Encoding of each vertex grid.
  /// - 0: Below interval volume.
  /// - 1: In interval volume. Vertex lifted to one +.
  /// - 2: In interval volume. Vertex lifted to two +'s.
  /// - 3: Above interval volume.
  /// @param[out] ivolpoly[] Vertices of interval volume polytopes (hexahedra).
  /// - Poly i has vertices ivolpoly[i*8], ivolpoly[i*8+1], ..., 
  ///     ivolpoly[i*8+7].
  /// - Vertices are given by the index of the cube containing the vertex.
  /// @param[out] poly_vertex[] Index of vertex with respect to the containing
  ///     polytope.  Possible values are 0, 1, ..., 7.
  /// @param[out] ivolpoly_info[] Information about interval volume polytopes.
  ///     ivolpoly_info[i] contains information about polytope i,
  ///     such as whether it is dual to a grid edge or
  ///     a grid vertex, or whether it has reverse orientation.
  void extract_dual_ivolpoly
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   IVOLDUAL_INFO & dualiso_info);

  /// Extract interval volume polytopes dual to grid edges.
  void extract_ivolpoly_dual_to_grid_edges
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   IVOLDUAL_INFO & dualiso_info);

  /// Extract interval volume polytopes dual to grid vertices.
  void extract_ivolpoly_dual_to_grid_vertices
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   IVOLDUAL_INFO & dualiso_info);


  // **************************************************
  // SPLIT DUAL INTERVAL VOLUME VERTICES
  // **************************************************

  /// Split interval volume vertices in each cube.
  /// - Number of vertices in each cube is determined by the cube
  ///   configuration and the interval volume looku table.
  /// @param ivoldual_table Interval volume lookup table.
  void split_dual_ivolvert
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const std::vector<VERTEX_INDEX> & ivolpoly_cube, 
   const std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   std::vector<GRID_CUBE_DATA> & cube_list,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   std::vector<ISO_VERTEX_INDEX> & isopoly,
   int & num_split);


  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Cubes containing vertices have only one ambiguous facet.
  /// @param grid Grid of cubes.
  /// @param ivoldual_table Interval volume lookup table.
  /// @pre Must be a double lookup table containing both the table
  ///      which separates negative cube vertices and the table
  ///      which separates positive cube vertices.
  /// @param index_to_cube_list[] Index to array cube_list[].
  ///      - index_to_cube_list[icube] is the location in cube_list[]
  ///        containing cube icube.
  ///      - index_to_cube_list[icube] is defined only if cube icube
  ///        is active.
  /// @pre index_to_cube_list[] has size at least grid.NumVertices().
  /// @param cube_list[] List of active cubes.
  ///      - Routine may change cube_list[icube].table_index
  ///        for cubes with ambiguous facets.
  /// @param[out] num_split
  void split_non_manifold_ivolv_pairs_ambig
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const VERTEX_INDEX index_to_cube_list[],
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);

  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Version which creates array index_to_cube_list[].
  void split_non_manifold_ivolv_pairs_ambig
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);

  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Cubes containing vertices have only one ambiguous facet.
  /// - Version which allows more splits.
  /// -  Allow split if numv_in_lower_lifted > 1 but numv_in_upper_lifted = 1
  ///    or vice versa.
  void split_non_manifold_ivolv_pairs_ambigB
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const VERTEX_INDEX index_to_cube_list[],
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);

  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Version which creates array index_to_cube_list[].
  /// - Version which allows more splits.
  void split_non_manifold_ivolv_pairs_ambigB
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);

  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Cubes containing vertices have only one ambiguous facet.
  /// - Version which allows more splits.
  /// -  Allow split if numv_in_lower_lifted > 1 but numv_in_upper_lifted = 1
  ///    or vice versa.
  void split_non_manifold_ivolv_pairs_ambigC
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const VERTEX_INDEX index_to_cube_list[],
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);

  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Version which creates array index_to_cube_list[].
  /// - Version which allows more splits.
  void split_non_manifold_ivolv_pairs_ambigC
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);


  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Cubes containing vertices have only one ambiguous facet.
  /// - Version which allows more splits than split_non_manifold_ivolv_pairs_ambigD.
  /// -  Allow split if numv_in_lower_lifted > 1 but numv_in_upper_lifted = 1
  ///    or vice versa.
  /// -  Allow split if new configuration has numv_in_lower_lifted > 1 for both cubes.
  /// -  Allow split if new configuration has numv_in_upper_lifted > 1 for both cubes.
  void split_non_manifold_ivolv_pairs_ambigD
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const VERTEX_INDEX index_to_cube_list[],
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);

  /// Split interval volume vertex pairs which create non-manifold edges.
  /// - Version which creates array index_to_cube_list[].
  /// - Version which allows more splits.
  void split_non_manifold_ivolv_pairs_ambigD
  (const DUALISO_GRID & grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   std::vector<GRID_CUBE_DATA> & cube_list,
   int & num_split);

  // **************************************************
  // POSITION INTERVAL VOLUME VERTICES
  // **************************************************

  /// Position all the dual interval volume vertices
  /// at the centroid of the interpolated isosurface-grid edge
  /// intersection points.
  void position_all_dual_ivol_vertices
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_TYPE * vertex_coord);

  /// Position all the dual interval volume vertices.
  /// - Version with vertex coordinates stored in C++ STL vector vertex_coord.
  void position_all_dual_ivol_vertices
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord);

  /// Position interval volume vertex described by ivolv_info
  /// at the centroid of the interpolated isosurface-grid edge
  /// intersection points.
  void position_dual_ivolv_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const DUAL_IVOLVERT & ivolv_info,
   const CUBE_FACE_INFO & cube,
   COORD_TYPE * vcoord, 
   COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, 
   COORD_TYPE * temp_coord2);

  /// Position interval volume vertex on upper isosurface.
  void position_dual_ivolv_on_upper_isosurface_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue,
   const DUAL_IVOLVERT & ivolv_info,
   const CUBE_FACE_INFO & cube,
   COORD_TYPE * vcoord, 
   COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, 
   COORD_TYPE * temp_coord2);

  /// Position interval volume vertex on lower isosurface.
  void position_dual_ivolv_on_lower_isosurface_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue,
   const DUAL_IVOLVERT & ivolv_info,
   const CUBE_FACE_INFO & cube,
   COORD_TYPE * vcoord, 
   COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, 
   COORD_TYPE * temp_coord2);

  /// Position interval volume vertex described by ivolv_info
  /// in the volume interior.
  void position_dual_ivolv_in_interval_volume_centroid_multi
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const DUAL_IVOLVERT & ivolv_info,
   const CUBE_FACE_INFO & cube,
   COORD_TYPE * vcoord, 
   COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, 
   COORD_TYPE * temp_coord2);

}

#endif
