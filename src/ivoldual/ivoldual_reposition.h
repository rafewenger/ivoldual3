/// @file ivoldual_reposition.h
/// @brief Reposition vertices for mesh quality optimization.
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

#ifndef _IVOLDUAL_REPOSITION_
#define _IVOLDUAL_REPOSITION_

#include "ijk.tpp"

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"
#include "ivoldualtable.h"
#include <unordered_map>

/// ivoldual classes and routines.
namespace IVOLDUAL {

  // **************************************************
  // MESH OPTIMIZATION
  // **************************************************
  
  /// Eliminate non-manifold facets and cubes
  void eliminate_non_manifold_grid
  (IVOLDUAL_DATA & ivoldual_data, 
   const SCALAR_TYPE isovalue0, 
   const SCALAR_TYPE isovalue1, 
   IVOLDUAL_INFO & dualiso_info);

  /// Laplacian Smoothing for small edge length.
  void laplacian_smooth_elength
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord,
   float laplacian_smooth_limit, 
   int iteration);

  /// Laplacian Smoothing for bad Jacobian.
  void laplacian_smooth_jacobian
  (const std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord, 
   float jacobian_limit, 
   int iteration);

  void laplacian_smooth_jacobian
  (const std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord, 
   const std::vector<int> & negative_jabocian_list);

  /// Gradient Smoothing for bad Jacobian.
  void gradient_smooth_jacobian
  (std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   DUAL_IVOLVERT_ARRAY & ivolv_list,
   IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   COORD_ARRAY & vertex_coord, 
   float jacobian_limit, 
   int iteration);

  void expand_flat_hex
  (const std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord, 
   const std::vector<std::vector<int>> & flag_hex);

  void expand_flat_hex_normal_direction
  (const std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord, 
   const std::vector<int> & hex_list);

  void gradient_smooth_jacobian
  (const std::vector<VERTEX_INDEX> & ivolpoly_cube,
   const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   const DUAL_IVOLVERT_ARRAY & ivolv_list,
   COORD_ARRAY & vertex_coord, 
   const std::vector<int> & negative_jabocian_list,
   std::unordered_map<int, COORD_TYPE> & negative_jacobian_value,
   int iter);

  void gradient_move_vertex
	(const std::vector<VERTEX_INDEX> & ivolpoly_cube,
	 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
	 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
	 const DUAL_IVOLVERT_ARRAY & ivolv_list,
	 COORD_ARRAY & vertex_coord, 
	 COORD_TYPE *ver_coord, int ver_index,	 
	 bool flag_onLower,  bool flag_onUpper);

  void move_vertex_along_edge
  (const std::vector<VERTEX_INDEX> & ivolpoly_vert,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   COORD_ARRAY & vertex_coord,
   COORD_TYPE cur_min_jacob,
   int ivert);

  void move_vertex_all_direction
  (const std::vector<VERTEX_INDEX> & ivolpoly_vert,
   IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   COORD_ARRAY & vertex_coord,
   COORD_TYPE cur_min_jacob,
   int ivert, int iter);

  void min_jacob_around_vertex
  (const std::vector<VERTEX_INDEX> & ivolpoly_vert,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   COORD_ARRAY & vertex_coord,
   int & ivert,
   COORD_TYPE & min_jacob_at_cur,
   COORD_TYPE & min_jacob_around_cur);

  void move_vertex_normal_direction
  (const std::vector<VERTEX_INDEX> & ivolpoly_vert,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   COORD_ARRAY & vertex_coord,
   const std::vector<int> & internal_vert,
   COORD_TYPE normal_dir[]);

  void find_optimal_jacobian_point
  (const std::vector<VERTEX_INDEX> & ivolpoly_vert,
   IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
   COORD_ARRAY & vertex_coord, 
   const std::vector<std::vector<int>> & flat_hex, 
   float pre_min_at_facet, float pre_min_around_facet, 
   int ifacet, 
   float & move_dist, std::vector<int> & dir);
   
  void surface_normal_direction
   (COORD_ARRAY & vertex_coord,
    const std::vector<int> & surface_vert,
    COORD_TYPE dir[]);
}

#endif
