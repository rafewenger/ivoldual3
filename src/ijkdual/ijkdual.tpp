/// @file ijkdual.tpp
/// @brief Templates for dual contouring isosurface generation 
///   in arbitrary dimensions.
/// - Version 0.3.1

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

#ifndef _IJKDUAL_TPP_
#define _IJKDUAL_TPP_

#include <string>

#include "ijk.tpp"
#include "ijktime.tpp"
#include "ijkdual_position.tpp"
#include "ijkfacet_intersection_table.tpp"

#include "ijkdual_extract.tpp"

#include "ijkdual_types.h"
#include "ijkdual_datastruct.h"


namespace {

  // *** Local (old) version of split_dual_isovert.
  // *** Global (new) version requires list_loc_of_adjacent_cube
  // *** in GRID_CUBE_ISOVERT.
  /// Split dual isosurface vertices.
  /// @param isopoly_cube[] isopoly_cube[j*num_polyv+k] is the cube containing
  ///    the k'th vertex of isosurface polytope j
  /// @param facet_vertex[] facet_vertex[j*num_polyv+k] is the location
  ///    of the k'th vertex on isosurface polytope j.
  /// @param cube_list[] List of cubes containing isosurface vertices.
  /// @param iso_vlist[] List of isosurface vertices.
  ///    isosurface patch in iso_vlist_cube[i] containing vertex i.
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE,
            typename CINDEX_TYPE, typename FACETV_TYPE,
            typename GRID_CUBE_TYPE, 
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE,  typename NTYPE>
  void split_dual_isovert_local
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE> & isopoly_cube, 
   const std::vector<FACETV_TYPE> & facet_vertex,
   std::vector<GRID_CUBE_TYPE> & cube_list,
   std::vector<ISOV_TYPE> & iso_vlist, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split)
  {
    compute_cube_isotable_info
      (scalar_grid, isodual_table, isovalue, cube_list);

    IJK::split_dual_isovert
      (isodual_table, isopoly_cube, facet_vertex, cube_list, 
       iso_vlist, isopoly, num_split);

  }
  
}

/// ijkdual template functions
namespace IJKDUAL {

  // **************************************************
  // CONSTRUCT DUAL CONTOURING ISOSURFACE MESH
  // **************************************************

  /// Construct isosurface mesh using Dual Contouring algorithm.
  /// Allow multiple isosurface vertices per grid cube.
  /// Returns list of isosurface polytope vertices.
  template <typename GRID_EDGE_TYPE, typename GRID_CUBE_DATA_TYPE,
            typename DUAL_ISOVERT_TYPE>
  void construct_multi_isov_mesh
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_EDGE_TYPE> & dual_edge,
   std::vector<GRID_CUBE_DATA_TYPE> & cube_isov_list,
   std::vector<DUAL_ISOVERT_TYPE> & iso_vlist,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info)
  {
    const int dimension = scalar_grid.Dimension();
    const bool flag_split_non_manifold = param.SplitNonManifoldFlag();
    const bool flag_select_split = param.SelectSplitFlag();
    const bool flag_separate_neg = param.SeparateNegFlag();
    const bool flag_connect_ambiguous = param.ConnectAmbiguousFlag();
    const VERTEX_POSITION_METHOD vpos_method = param.VertexPositionMethod();
    IJK::PROCEDURE_ERROR error("construct_multi_isov_mesh");
    clock_t t0, t1, t2;

    if (scalar_grid.Dimension() != isodual_table.Dimension()) {
      error.AddMessage
        ("Programming error.  Incorrect isodual table dimension.");
      error.AddMessage
        ("  Isodual table dimension does not match scalar grid dimension.");
      error.AddMessage
        ("    Isodual table dimension: ", isodual_table.Dimension(), "");
      error.AddMessage
        ("    Scalar grid dimension: ", scalar_grid.Dimension(), "");
      throw error;
    }

    t0 = clock();

    isopoly_vert.clear();
    dualiso_info.time.Clear();

    std::vector<ISO_VERTEX_INDEX> isopoly;
    std::vector<FACET_VERTEX_INDEX> facet_vertex;
    extract_dual_isopoly
      (scalar_grid, isovalue, isopoly, facet_vertex, dual_edge, dualiso_info);
    t1 = clock();

    std::vector<ISO_VERTEX_INDEX> cube_list;
    std::vector<ISO_VERTEX_INDEX> isopoly_cube;
    merge_identical(isopoly, cube_list, isopoly_cube, merge_data);
    t2 = clock();

    set_grid_cube_indices(cube_list, cube_isov_list);

    VERTEX_INDEX num_split;
    if (flag_split_non_manifold || flag_select_split ||
        flag_connect_ambiguous) {

      int num_non_manifold_split = 0;
      int num_1_2_changed = 0;
      int num_connect_changed = 0;
      int num_ambig_ridge_cubes_changed = 0;
      int num_non_ambig_ridge_cubes_changed = 0;

      IJK::compute_cube_isotable_info
        (scalar_grid, isodual_table, isovalue, cube_isov_list);

      if (flag_split_non_manifold) {
        IJK::FACET_INTERSECTION_TABLE
          <DIRECTION_TYPE,VERTEX_INDEX,BOUNDARY_BITS_TYPE>
          facet_intersection_table(dimension);
        facet_intersection_table.Create();

        // Recompute isotable indices for cubes on grid ridges
        //   to avoid non-manifold vertices.
        IJK::recompute_cube_isotable_info_on_grid_ridges
          (scalar_grid, isodual_table, isovalue, facet_intersection_table,
           flag_separate_neg, cube_isov_list, num_ambig_ridge_cubes_changed,
           num_non_ambig_ridge_cubes_changed);
      }

      if (flag_split_non_manifold) {
        IJK::split_non_manifold_isov_pairs_using_index
          (scalar_grid, isodual_table, cube_isov_list, num_non_manifold_split);
      }

      if (flag_select_split) {
        IJK::select_split_1_2_ambig
          (scalar_grid, isodual_table, isovalue, cube_isov_list, 
           num_1_2_changed);
      }

      if (flag_connect_ambiguous) {
        IJK::select_ambig_to_connect_isosurface
          (scalar_grid, isodual_table, isovalue, cube_isov_list, 
           num_connect_changed);
      }

      IJK::split_dual_isovert
        (isodual_table, isopoly_cube, facet_vertex, cube_isov_list, 
         iso_vlist, isopoly_vert, num_split);

      dualiso_info.multi_isov.num_non_manifold_split = num_non_manifold_split;
      dualiso_info.multi_isov.num_1_2_changed = num_1_2_changed;
      dualiso_info.multi_isov.num_connect_changed = num_connect_changed;
      dualiso_info.multi_isov.num_ambig_ridge_cubes_changed =
        num_ambig_ridge_cubes_changed;
      dualiso_info.multi_isov.num_non_ambig_ridge_cubes_changed =
        num_non_ambig_ridge_cubes_changed;
    }
    else {
      split_dual_isovert
        (scalar_grid, isodual_table, isovalue, 
         isopoly_cube, facet_vertex, cube_isov_list,
         iso_vlist, isopoly_vert, num_split);
    }

    dualiso_info.scalar.num_non_empty_cubes = cube_list.size();
    dualiso_info.multi_isov.num_cubes_multi_isov = num_split;
    dualiso_info.multi_isov.num_cubes_single_isov =
      cube_list.size() - num_split;

    // store times
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
    IJK::clock2seconds(t2-t1, dualiso_info.time.merge);
  }


  // **************************************************
  // DUAL CONTOURING (HYPERCUBES)
  // **************************************************

  /// Extract isosurface using Dual Contouring algorithm.
  /// Allow multiple isosurface vertices per grid cube.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates.
  template <typename GRID_EDGE_TYPE, typename GRID_CUBE_DATA_TYPE,
            typename DUAL_ISOVERT_TYPE>
  void dual_contouring_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_EDGE_TYPE> & dual_edge,
   std::vector<GRID_CUBE_DATA_TYPE> & cube_isov_list,
   std::vector<DUAL_ISOVERT_TYPE> & iso_vlist,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info)
  {
    const VERTEX_POSITION_METHOD vpos_method = param.VertexPositionMethod();
    const COORD_TYPE center_offset = 0.1;
    clock_t t0, t1, t2;

    t0 = clock();

    construct_multi_isov_mesh
      (scalar_grid, isovalue, isodual_table, param, isopoly_vert,
       dual_edge, cube_isov_list, iso_vlist, merge_data, dualiso_info);

    t1 = clock();

    if (vpos_method == CUBE_CENTER) {
      position_all_dual_isovertices_near_cube_center_multi
        (scalar_grid, isodual_table, isovalue, iso_vlist, center_offset, 
         vertex_coord);
    }
    else if (vpos_method == IVOL_LIFTED02) {
      SCALAR_TYPE isovalue0 = isovalue/2;
      SCALAR_TYPE isovalue1 = (isovalue+2)/2;
      position_all_dual_isovertices_ivol_lifted
        (scalar_grid, isodual_table, isovalue0, isovalue1, 
         iso_vlist, vertex_coord);
    }
    else {
      position_all_dual_isovertices_centroid_multi
        (scalar_grid, isodual_table, isovalue, iso_vlist, vertex_coord);
    }

    t2 = clock();

    // store times
    IJK::clock2seconds(t2-t1, dualiso_info.time.position);
    IJK::clock2seconds(t2-t0, dualiso_info.time.total);
  }

  /// Extract isosurface using Dual Contouring algorithm.
  /// Allow multiple isosurface vertices per grid cube.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates
  /// Version without argument cube_isov_list.
  template <typename GRID_EDGE_TYPE, typename DUAL_ISOVERT_TYPE>
  void dual_contouring_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_EDGE_TYPE> & dual_edge,
   std::vector<DUAL_ISOVERT_TYPE> & iso_vlist,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info)
  {
    std::vector<GRID_CUBE_DATA> cube_isov_list;

    dual_contouring_multi_isov
      (scalar_grid, isovalue, isodual_table, param, isopoly_vert,
       dual_edge, cube_isov_list, iso_vlist, vertex_coord,
       merge_data, dualiso_info);
  }

  /// Extract isosurface using Dual Contouring algorithm.
  /// Allow multiple isosurface vertices per grid cube.
  /// Version which creates isodual_table and merge_data.
  /// Version which does not return cube_isov_list or dualiso_info.
  template <typename GRID_EDGE_TYPE, typename DUAL_ISOVERT_TYPE>
  void dual_contouring_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_EDGE_TYPE> & dual_edge,
   std::vector<DUAL_ISOVERT_TYPE> & iso_vlist,
   COORD_ARRAY & vertex_coord)
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
    const bool flag_separate_neg = param.SeparateNegFlag();
    const bool flag_always_separate_opposite(true);
    ISODUAL_CUBE_TABLE_AMBIG 
      isodual_table(dimension, flag_separate_neg, 
                    flag_always_separate_opposite);
    ISO_MERGE_DATA merge_data(dimension, axis_size);
    DUALISO_INFO dualiso_info;
    std::vector<GRID_CUBE_DATA> cube_isov_list;

    dual_contouring_multi_isov
      (scalar_grid, isovalue, isodual_table, param, isopoly_vert,
       dual_edge, cube_isov_list, iso_vlist, vertex_coord,
       merge_data, dualiso_info);
  }

  /// Extract isosurface using Manifold Dual Contouring algorithm.
  /// Isosurface must be 2D surface in 3D volume.
  /// Allow multiple isosurface vertices per grid cube.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates.
  /// Isosurface mesh has a manifold representation.
  template <typename GRID_EDGE_TYPE, typename GRID_CUBE_DATA_TYPE, 
            typename DUAL_ISOVERT_TYPE>
  void dual_contouring_manifold
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_EDGE_TYPE> & dual_edge,
   std::vector<GRID_CUBE_DATA_TYPE> & cube_isov_list,
   std::vector<DUAL_ISOVERT_TYPE> & iso_vlist,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info)
  {
    DUALISO_DATA_FLAGS manifold_param = param;

    manifold_param.allow_multiple_iso_vertices = true;
    manifold_param.flag_split_non_manifold = true;

    dual_contouring_multi_isov
      (scalar_grid, isovalue, isodual_table, manifold_param, 
       isopoly_vert, dual_edge, cube_isov_list, iso_vlist,
       vertex_coord, merge_data, dualiso_info);
  }

  /// Extract isosurface using Manifold Dual Contouring algorithm.
  /// Isosurface must be 2D surface in 3D volume.
  /// Allow multiple isosurface vertices per grid cube.
  /// Returns list of isosurface polytope vertices
  ///   and list of isosurface vertex coordinates.
  /// Isosurface mesh has a manifold representation.
  /// Version without cube_isov_list.
  template <typename GRID_EDGE_TYPE, typename DUAL_ISOVERT_TYPE>
  void dual_contouring_manifold
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue, 
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_FLAGS & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_EDGE_TYPE> & dual_edge,
   std::vector<DUAL_ISOVERT_TYPE> & iso_vlist,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data, 
   DUALISO_INFO & dualiso_info)
  {
    DUALISO_DATA_FLAGS manifold_param = param;

    manifold_param.allow_multiple_iso_vertices = true;
    manifold_param.flag_split_non_manifold = true;

    dual_contouring_multi_isov
      (scalar_grid, isovalue, isodual_table, manifold_param, 
       isopoly_vert, dual_edge, iso_vlist,
       vertex_coord, merge_data, dualiso_info);
  }

}

#endif
