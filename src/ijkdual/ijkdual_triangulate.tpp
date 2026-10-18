/// @file ijkdual_triangulate.tpp
/// @brief Triangulation routines.
/// - Version 0.3.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2016-2026 Rephael Wenger

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


#ifndef _IJKDUAL_TRIANGULATE_TPP_
#define _IJKDUAL_TRIANGULATE_TPP_

#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijkisocoord.tpp"
#include "ijktriangulate_geom.txx"

#include "ijkdual_types.h"


namespace IJKDUAL {

  // **************************************************
  // COMPUTE INTERSECTION/CENTROID
  // **************************************************

  /// Compute intersection of quad and grid edge.
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename GRID_EDGE_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1>
  void compute_intersection_of_quad_and_grid_edge
  (const GRID_TYPE & grid, 
   const ISOV_INDEX_TYPE quad_vert[],
   const GRID_EDGE_TYPE & grid_edge,
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int DIM3(3);
    const int NUM_VERT_PER_QUAD(4);
    COORD_TYPE0 coord0[DIM3], coord1[DIM3];
    COORD_TYPE0 coord01[DIM3], coord23[DIM3];
    const COORD_TYPE0 * isov_coord[NUM_VERT_PER_QUAD];

    
    VTYPE iv0 = grid_edge.Endpoint0();
    DTYPE edge_direction = grid_edge.Direction();
    VTYPE iv1 = grid.NextVertex(iv0, edge_direction);

    grid.ComputeCoord(iv0, coord0);
    grid.ComputeCoord(iv1, coord1);

    const int dir1 = (edge_direction+1)%DIM3;
    const int dir2 = (edge_direction+2)%DIM3;

    for (NTYPE i = 0; i < NUM_VERT_PER_QUAD; i++) {
      const ISOV_INDEX_TYPE isov = quad_vert[i];
      isov_coord[i] = &(coord_array[isov*DIM3]);
    }

    DTYPE dir01, dir12;
    if (IJK::is_dth_coord_in_range
        (dir1, coord0, isov_coord[0], isov_coord[1]) &&
        IJK::is_dth_coord_in_range
        (dir1, coord0, isov_coord[2], isov_coord[3])) {
      dir01 = dir1;
      dir12 = dir2;
    }
    else {
      dir01 = dir2;
      dir12 = dir1;
    }

    IJK::compute_line_segment_hyperplane_intersection
      (DIM3, coord_array, quad_vert[0], quad_vert[1],
       dir01, coord0[dir01], coord01);
    IJK::compute_line_segment_hyperplane_intersection
      (DIM3, coord_array, quad_vert[2], quad_vert[3],
       dir01, coord0[dir01], coord23);
    IJK::compute_line_segment_hyperplane_intersection
      (DIM3, coord01, coord23, dir12, coord0[dir12], intersection_coord);
  }


  /// Compute intersection of quad and grid edge.
  /// - Version with array of quad vertices.
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename GRID_EDGE_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1, typename QUAD_INDEX_TYPE>
  void compute_intersection_of_quad_and_grid_edge
  (const GRID_TYPE & grid, 
   const std::vector<ISOV_INDEX_TYPE> & quad_vert,
   const GRID_EDGE_TYPE & grid_edge,
   std::vector<COORD_TYPE0> & coord_array,
   const QUAD_INDEX_TYPE iquad,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE NUM_VERT_PER_QUAD(4);

    const ISOV_INDEX_TYPE * first_quad_vert = 
      &(quad_vert[iquad*NUM_VERT_PER_QUAD]);

    compute_intersection_of_quad_and_grid_edge
      (grid, first_quad_vert, grid_edge, coord_array, intersection_coord);

  }


  // **************************************************
  // ADD ISOSURFACE VERTICES ON GRID EDGES
  // **************************************************

  /// Add isosurface vertices on grid edges.
  template <typename DATA_TYPE, typename STYPE, typename ISOSURFACE_TYPE>
  void add_isov_on_grid_edges
  (const DATA_TYPE & dualiso_data, const STYPE isovalue,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    const QUAD_EDGE_INTERSECTION_METHOD qei_method =
      dualiso_data.quad_edge_intersection_method;

    if (qei_method == QEI_INTERPOLATE_COORD) {
      add_isov_on_grid_edges_interpolate_coord(dualiso_data, dual_isosurface);
    }
    else {
      add_isov_on_grid_edges_interpolate_scalar
        (dualiso_data, isovalue, dual_isosurface);
    }

  }


  /// Add isosurface vertices on grid edges using interpolation 
  ///   on scalar values.
  template <typename DATA_TYPE, typename STYPE, typename ISOSURFACE_TYPE>
  void add_isov_on_grid_edges_interpolate_scalar
  (const DATA_TYPE & dualiso_data, const STYPE isovalue,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const int numv_per_iso_poly = dual_isosurface.NumVerticesPerIsoPoly();
    const SIZE_TYPE num_poly = 
      dual_isosurface.isopoly_vert.size()/numv_per_iso_poly;
    IJK::PROCEDURE_ERROR error("add_isov_on_grid_edges_interpolate_scalar");

    const SIZE_TYPE first_isov_on_edge = dual_isosurface.NumIsoVert();

    if (!dualiso_data.CheckIsoQuadDualToGridEdges
        ("add_isov_on_grid_edges_interpolate_scalar", error))
      { throw error; }
    if (!dual_isosurface.CheckIsopolyInfo(error)) { throw error; }

    dual_isosurface.first_isov_dual_to_iso_poly = first_isov_on_edge;

    dual_isosurface.vertex_coord.resize
      (dimension*(first_isov_on_edge+num_poly));

    if (dual_isosurface.vertex_coord.size() > 0) {

      COORD_TYPE * first_isov_on_edge_coord =
        IJK::vector2pointerNC(dual_isosurface.vertex_coord)
        + dimension*first_isov_on_edge;

      /* OBSOLETE
      IJK::compute_isov_coord_on_grid_edge_linear
        (dualiso_data.ScalarGrid(), isovalue,
         dual_isosurface.isopoly_info, first_isov_on_edge_coord);
         */

      IJK::compute_all_isov_coord_on_grid_edge_linear
        (dualiso_data.ScalarGrid(), isovalue,
         dual_isosurface.isopoly_info, first_isov_on_edge_coord);
    }

  }


  /// Add isosurface vertices on grid edges using interpolation
  ///   of quadrilateral vertex coordinates.
  template <typename DATA_TYPE, typename ISOSURFACE_TYPE>
  void add_isov_on_grid_edges_interpolate_coord
  (const DATA_TYPE & dualiso_data, ISOSURFACE_TYPE & dual_isosurface)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const int numv_per_iso_poly = dual_isosurface.NumVerticesPerIsoPoly();
    const int num_poly = dual_isosurface.isopoly_vert.size()/numv_per_iso_poly;
    IJK::PROCEDURE_ERROR error("add_isov_on_grid_edges_interpolate_coord");

    const SIZE_TYPE first_isov_on_edge = dual_isosurface.NumIsoVert();

    if (!dualiso_data.CheckIsoQuadDualToGridEdges
        ("add_isov_on_grid_edges_interpolate_coord", error))
      { throw error; }
    if (!dual_isosurface.CheckIsopolyInfo(error)) { throw error; }

    dual_isosurface.first_isov_dual_to_iso_poly = first_isov_on_edge;

    dual_isosurface.vertex_coord.resize
      (dimension*(first_isov_on_edge+num_poly));

    if (dual_isosurface.vertex_coord.size() > 0) {

      COORD_TYPE * first_isov_on_edge_coord =
        IJK::vector2pointerNC(dual_isosurface.vertex_coord)
        + dimension*first_isov_on_edge;

      for (int ipoly = 0; ipoly < num_poly; ipoly++) {

        const ISO_VERTEX_INDEX isov = first_isov_on_edge + ipoly;
        COORD_TYPE * isov_coord =
          &(dual_isosurface.vertex_coord[isov*dimension]);

        compute_intersection_of_quad_and_grid_edge
          (dualiso_data.ScalarGrid(), dual_isosurface.isopoly_vert,
           dual_isosurface.isopoly_info[ipoly], dual_isosurface.vertex_coord, 
           ipoly, isov_coord);
      }
    }
  }


  // ***********************************************************************
  // ADD NEW ISOSURFACE VERTICES AT CENTROIDS OF ISOSURFACE POLY VERTICES
  // ***********************************************************************


  /// Add new isosurface vertices at centroid of isosurface poly vertices.
  template <typename DATA_TYPE, typename ISOSURFACE_TYPE>
  void add_isov_at_poly_centroids
  (const DATA_TYPE & dualiso_data, ISOSURFACE_TYPE & dual_isosurface)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const int numv_per_iso_poly = dual_isosurface.NumVerticesPerIsoPoly();
    const int num_poly = dual_isosurface.isopoly_vert.size()/numv_per_iso_poly;
    IJK::PROCEDURE_ERROR error("add_isov_at_poly_centroids");

    const SIZE_TYPE first_new_isov = dual_isosurface.NumIsoVert();
    const SIZE_TYPE num_new_isov = num_poly;

    dual_isosurface.first_isov_dual_to_iso_poly = first_new_isov;

    dual_isosurface.vertex_coord.resize
      (dimension*(first_new_isov+num_new_isov));

    if (dual_isosurface.vertex_coord.size() > 0) {

      for (int ipoly = 0; ipoly < num_poly; ipoly++) {

        const ISO_VERTEX_INDEX isov = first_new_isov + ipoly;
        COORD_TYPE * isov_coord =
          &(dual_isosurface.vertex_coord[isov*dimension]);

        const ISO_VERTEX_INDEX * first_isopoly_vert = 
          &(dual_isosurface.isopoly_vert[ipoly*numv_per_iso_poly]);
        
        IJK::compute_quad_centroid
          (dimension, first_isopoly_vert, dual_isosurface.vertex_coord, 
           isov_coord);
      }
    }

  }


  /// Add new vertex coordinates at centroid of pentagon vertices.
  template <typename VTYPE, typename NTYPE, typename CTYPE>
  void add_vertex_coord_at_pentagon_centroids
  (const VTYPE pentagon_vert[],
   const NTYPE num_pentagon,
   std::vector<CTYPE> & vertex_coord)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int DIM3(3);
    const int NUM_VERT_PER_PENTAGON(5);
    const SIZE_TYPE vertex_coord_size = vertex_coord.size();

    vertex_coord.resize(vertex_coord_size + DIM3*num_pentagon);

    for (int ipoly = 0; ipoly < num_pentagon; ipoly++) {

      COORD_TYPE * new_coord =
        &(vertex_coord[vertex_coord_size + ipoly*DIM3]);

      const ISO_VERTEX_INDEX * first_pentagon_vert = 
        &(pentagon_vert[ipoly*NUM_VERT_PER_PENTAGON]);
        
      IJK::compute_pentagon_centroid
        (DIM3, first_pentagon_vert, vertex_coord, new_coord);
    }

  }

  /// Add new vertex coordinates at centroid of pentagon vertices.
  /// - Version using C++ STL vector for array pentagon_vert[].
  template <typename VTYPE, typename CTYPE>
  void add_vertex_coord_at_pentagon_centroids
  (const std::vector<VTYPE> & pentagon_vert,
   std::vector<CTYPE> & vertex_coord)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const int num_pentagon = pentagon_vert.size()/NUM_VERT_PER_PENTAGON;

    add_vertex_coord_at_pentagon_centroids
      (IJK::vector2pointer(pentagon_vert), num_pentagon, vertex_coord);
  }


  /// Add new vertex coordinates at centroid of hexahedra vertices.
  template <typename DTYPE, typename VTYPE, 
            typename NTYPE, typename CTYPE>
  void add_vertex_coord_at_hexahedra_centroids
  (const DTYPE dimension,
   const VTYPE hexahedra_vert[],
   const NTYPE num_hexahedra,
   std::vector<CTYPE> & vertex_coord)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_HEXAHEDRA(8);
    const SIZE_TYPE vertex_coord_size = vertex_coord.size();

    vertex_coord.resize(vertex_coord_size + dimension*num_hexahedra);

    for (int ipoly = 0; ipoly < num_hexahedra; ipoly++) {

      COORD_TYPE * new_coord =
        &(vertex_coord[vertex_coord_size + ipoly*dimension]);

      const ISO_VERTEX_INDEX * first_hexahedron_vert = 
        &(hexahedra_vert[ipoly*NUM_VERT_PER_HEXAHEDRA]);
        
      IJK::compute_hexahedron_centroid
        (dimension, first_hexahedron_vert, vertex_coord, new_coord);
    }

  }


  /// Add new vertex coordinates at centroid of hexahedra vertices.
  template <typename DATA_TYPE, typename ISOSURFACE_TYPE>
  void add_vertex_coord_at_hexahedra_centroids
  (const DATA_TYPE & dualiso_data, ISOSURFACE_TYPE & dual_isosurface)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const int NUM_VERT_PER_HEXAHEDRA(8);
    const SIZE_TYPE first_new_isov = dual_isosurface.NumIsoVert();
    const SIZE_TYPE num_hexahedra = 
      dual_isosurface.isopoly_vert.size()/NUM_VERT_PER_HEXAHEDRA;

    dual_isosurface.first_isov_dual_to_iso_poly = first_new_isov;

    const ISO_VERTEX_INDEX * first_hex_vert = 
      &(dual_isosurface.isopoly_vert.front());

    add_vertex_coord_at_hexahedra_centroids
      (dimension, first_hex_vert, num_hexahedra, 
       dual_isosurface.vertex_coord);
  }


  // ***************************************************************
  // CONVERT QUADRILATERALS TO TRIANGLES
  // ***************************************************************

  /// Convert quadrilaterals (embedded in 3D) to triangles.
  /// Assumes input isosurface is quadrilaterals embedded in 3D.
  template <typename DATA_TYPE, typename ISOSURFACE_TYPE>
  void convert_quad_to_tri
  (const DATA_TYPE & dualiso_data, ISOSURFACE_TYPE & dual_isosurface)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_QUAD = 4;
    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const int numv_per_iso_poly = dual_isosurface.NumVerticesPerIsoPoly();
    const COORD_TYPE max_small_magnitude = dualiso_data.MaxSmallMagnitude();
    const QUAD_TRI_METHOD quad_tri_method = 
      dualiso_data.QuadTriangulationMethod();
    VERTEX_INDEX_ARRAY quad_vert(dual_isosurface.isopoly_vert);
    IJK::PROCEDURE_ERROR error("convert_quad_to_tri");

    if (numv_per_iso_poly != NUM_VERT_PER_QUAD) {
      error.AddMessage("Number of vertices per isosurface polygon is ",
                       numv_per_iso_poly, ".");
      error.AddMessage("Number of vertices per isosurface polygon must be ",
                       NUM_VERT_PER_QUAD, ".");
      throw error;
    }

    IJK::reorder_quad_vertices(quad_vert);

    if (quad_tri_method == MAX_MIN_ANGLE) {

      IJK::triangulate_quad_max_min_angle
        (dimension, dual_isosurface.vertex_coord, quad_vert, 
         max_small_magnitude, dual_isosurface.tri_vert);
    }
    else if (quad_tri_method == SPLIT_MAX_ANGLE) {
      IJK::triangulate_quad_split_max_angle
        (dimension, dual_isosurface.vertex_coord, quad_vert, 
         max_small_magnitude, dual_isosurface.tri_vert);
    }
    else if (dualiso_data.flag_tri4_quad) {
      const VERTEX_INDEX num_iso_vert = dual_isosurface.NumIsoVert();
      const VERTEX_INDEX num_iso_poly = dual_isosurface.NumIsoPoly();

      if (num_iso_vert != num_iso_poly+dual_isosurface.first_isov_dual_to_iso_poly) {
        error.AddMessage
          ("Programming error.  Incorrect number of isosurface vertices.");
        error.AddMessage
          ("  Check that add_isov_on_grid_edges is called before convert_quad_to_tri.");
        throw error;
      }
    
      if (dual_isosurface.vertex_coord.size() > 0) {

        const ISO_VERTEX_INDEX first_isov_dual_to_iso_poly =
          dual_isosurface.first_isov_dual_to_iso_poly;
        const COORD_TYPE min_distance_use_tri4 = 
          dualiso_data.MinDistanceUseTri4();
        const COORD_TYPE min_distance_allow_tri4 = 
          dualiso_data.MinDistanceAllowTri4();

        if (quad_tri_method == ONLY_TRI4) {
          IJK::triangulate_quad_with_vertices
            (quad_vert, first_isov_dual_to_iso_poly, dual_isosurface.tri_vert);
        }
        else if (quad_tri_method == TRI4_BY_DISTANCE) {
          triangulate_quad_tri4_by_distance
            (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord, 
             quad_vert, dual_isosurface.isopoly_info,
             first_isov_dual_to_iso_poly, min_distance_use_tri4, 
             max_small_magnitude, dual_isosurface.tri_vert);
        }
        else {


          if (dualiso_data.tri4_position_method == TRI4_CENTROID) {
            triangulate_quad_tri4_max_min_angle
              (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord, 
               quad_vert, first_isov_dual_to_iso_poly, max_small_magnitude, 
               dual_isosurface.tri_vert);
          }
          else {
            triangulate_dual_quad_tri4_max_min_angle
              (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord, 
               quad_vert, dual_isosurface.isopoly_info,
               first_isov_dual_to_iso_poly,
               min_distance_allow_tri4, max_small_magnitude, 
               dual_isosurface.tri_vert);
          }
        }
      }
    }
    else {
      IJK::triangulate_quad(quad_vert, dual_isosurface.tri_vert);
    }
  }

}

#endif
