/// @file ijkdual_position.txx
/// @brief Position dual isosurface vertices.
/// - Version 0.3.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2026 Rephael Wenger

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

#ifndef IJKDUAL_POSITION_TXX_
#define IJKDUAL_POSITION_TXX_

#include <vector>

#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijkisopoly.txx"
#include "ijkscalar_grid.tpp"


namespace IJKDUAL {

  // **************************************************
  // SINGLE ISOSURFACE VERTEX IN A GRID CUBE
  // **************************************************

  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE,
            typename CTYPE>
  /// Position dual isosurface vertices in cube centers
  void position_all_dual_isovertices_cube_center
  (const GRID_TYPE & grid,
   const std::vector<ISOV_INDEX_TYPE> & vlist, CTYPE * coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = grid.Dimension();

    for (VTYPE i = 0; i < vlist.size(); i++) {
      VTYPE iv = vlist[i];

      grid.ComputeCubeCenterCoord(iv, coord+i*dimension);
    }
  }

  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE,
            typename CTYPE>
  /// Position dual isosurface vertices in cube centers.
  /// C++ STL vector format for array coord[].
  void position_all_dual_isovertices_cube_center
  (const GRID_TYPE & grid,
   const std::vector<ISOV_INDEX_TYPE> & vlist, 
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();

    coord.resize(vlist.size()*dimension);
    position_all_dual_isovertices_cube_center(grid, vlist, &(coord.front()));
  }


  /// Position dual isosurface vertices in centroid 
  ///   of isosurface-edge intersections
  template <typename GRID_TYPE, typename STYPE, typename ISOV_INDEX_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ISOV_INDEX_TYPE> & vlist, CTYPE * coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> vcoord(dimension);
    IJK::ARRAY<CTYPE> coord0(dimension);
    IJK::ARRAY<CTYPE> coord1(dimension);
    IJK::ARRAY<CTYPE> coord2(dimension);

    for (VTYPE i = 0; i < vlist.size(); i++) {
      VTYPE iv = vlist[i];

      NTYPE num_intersected_edges = 0;
      IJK::set_coord(dimension, 0.0, vcoord.Ptr());

      for (DTYPE edge_dir = 0; edge_dir < dimension; edge_dir++)
        for (NTYPE k = 0; k < scalar_grid.NumFacetVertices(); k++) {
          VTYPE iend0 = scalar_grid.FacetVertex(iv, edge_dir, k);
          VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);

          STYPE s0 = scalar_grid.Scalar(iend0);
          bool is_end0_positive = true;
          if (s0 < isovalue)
            { is_end0_positive = false; };

          STYPE s1 = scalar_grid.Scalar(iend1);
          bool is_end1_positive = true;
          if (s1 < isovalue)
            { is_end1_positive = false; };

          if (is_end0_positive != is_end1_positive) {

            scalar_grid.ComputeCoord(iend0, coord0.Ptr());
            scalar_grid.ComputeCoord(iend1, coord1.Ptr());

            IJK::linear_interpolate_coord
              (dimension, s0, coord0.Ptr(), s1, coord1.Ptr(), 
               isovalue, coord2.Ptr());

            IJK::add_coord(dimension, vcoord.Ptr(), coord2.Ptr(), vcoord.Ptr());

            num_intersected_edges++;
          }

        }

      if (num_intersected_edges > 0) {
        IJK::multiply_coord
          (dimension, 1.0/num_intersected_edges, vcoord.Ptr(), 
           coord+i*dimension);
      }
      else {
        scalar_grid.ComputeCubeCenterCoord(iv, coord+i*dimension);
      }

    }
  }


  /// Position dual isosurface vertices in centroid 
  ///   of isosurface-edge intersections
  template <typename GRID_TYPE, typename STYPE, typename ISOV_INDEX_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ISOV_INDEX_TYPE> & vlist, 
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(vlist.size()*dimension);
    position_all_dual_isovertices_centroid
      (scalar_grid, isovalue, vlist, &(coord.front()));
  }


  // **************************************************
  // ALLOW MULTIPLE ISOSURFACE VERTICES IN A GRID CUBE
  // **************************************************

  /// Position dual isosurface vertex at centroid.
  /// @pre scalar_grid.Dimension() == cube.Dimension().
  /// @param[out] vcoord[] Vertex coordinates.
  /// @param temp_coord0[] Temporary coordinate array.
  /// @param temp_coord1[] Temporary coordinate array.
  /// @param temp_coord2[] Temporary coordinate array.
  template <typename GRID_TYPE, 
            typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE, 
            typename CUBE_TYPE,
            typename CTYPE, typename CTYPE2>
  void position_dual_isov_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const DUAL_ISOV_TYPE & isov_info,
   const CUBE_TYPE & cube,
   CTYPE * vcoord, 
   CTYPE2 * temp_coord0, CTYPE2 * temp_coord1, CTYPE2 * temp_coord2)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = scalar_grid.Dimension();
    const VTYPE icube = isov_info.cube_index;
    const NTYPE ipatch = isov_info.patch_index;
    const TABLE_INDEX it = isov_info.table_index;

    NTYPE num_intersected_edges = 0;
    IJK::set_coord(dimension, 0.0, vcoord);

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(it, ie)) {
        if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
          NTYPE k0 = cube.EdgeEndpoint(ie, 0);
          NTYPE k1 = cube.EdgeEndpoint(ie, 1);
          VTYPE iend0 = scalar_grid.CubeVertex(icube, k0);
          VTYPE iend1 = scalar_grid.CubeVertex(icube, k1);
          STYPE s0 = scalar_grid.Scalar(iend0);
          STYPE s1 = scalar_grid.Scalar(iend1);

          scalar_grid.ComputeCoord(iend0, temp_coord0);
          scalar_grid.ComputeCoord(iend1, temp_coord1);

          if ((s0 < isovalue && s1 < isovalue) ||
              (s0 > isovalue && s1 > isovalue)) {
            // Use edge midpoint.
            IJK::linear_interpolate_coord
              (dimension, 0.5, temp_coord0, temp_coord1, temp_coord2);
          }
          else {
            IJK::linear_interpolate_coord
              (dimension, s0, temp_coord0, s1, temp_coord1, isovalue, 
               temp_coord2);
          }

          IJK::add_coord(dimension, vcoord, temp_coord2, vcoord);

          num_intersected_edges++;
        }
      }
    }

    if (num_intersected_edges > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /// Position dual isosurface vertices using centroids.
  template <typename GRID_TYPE, 
            typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & iso_vlist,
   CTYPE * coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> coord0(dimension);
    IJK::ARRAY<CTYPE> coord1(dimension);
    IJK::ARRAY<CTYPE> coord2(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    for (NTYPE i = 0; i < iso_vlist.size(); i++) {
      position_dual_isov_centroid_multi
        (scalar_grid, isodual_table, isovalue, iso_vlist[i], cube,
         coord+i*dimension, coord0.Ptr(), coord1.Ptr(), coord2.Ptr());
    }
  }

  /// Position dual isosurface vertices using centroids.
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & iso_vlist, 
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(iso_vlist.size()*dimension);
    position_all_dual_isovertices_centroid_multi
      (scalar_grid, isodual_table,isovalue, iso_vlist,
       &(coord.front()));
  }

  /// Position dual isosurface vertices near cube centers.
  /// More than one vertex can be in a cube.
  /// If cube contains multiple isosurface then vertices are positioned
  ///   near but not on cube center.
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE0, typename CTYPE1>
  void position_all_dual_isovertices_near_cube_center_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & iso_vlist,
   const CTYPE0 offset,
   CTYPE1 * coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = scalar_grid.Dimension();
    CTYPE1 vcoord[dimension];
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    IJK::UNIT_CUBE<int,int,int> unit_cube(dimension);
    bool intersects_facet[2*dimension];

    for (VTYPE i = 0; i < iso_vlist.size(); i++) {
      VTYPE icube = iso_vlist[i].cube_index;
      NTYPE ipatch = iso_vlist[i].patch_index;
      TABLE_INDEX it = iso_vlist[i].table_index;

      if (isodual_table.NumIsoVertices(it) == 1) {
        scalar_grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
      }
      else {

        scalar_grid.ComputeCubeCenterCoord(icube, vcoord);

        // Set intersects facet to false.
        for (DTYPE d = 0; d < dimension; d++) {
          intersects_facet[2*d] = false;
          intersects_facet[2*d+1] = false;
        }

        for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
          if (isodual_table.IsBipolar(it, ie)) {
            if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
              NTYPE k0 = cube.EdgeEndpoint(ie, 0);
              NTYPE k1 = cube.EdgeEndpoint(ie, 1);

              for (DTYPE d = 0; d < dimension; d++) {
                NTYPE c = unit_cube.VertexCoord(k0, d);
                if (c == unit_cube.VertexCoord(k1, d)) {
                  if (c > 0) 
                    { intersects_facet[2*d+1] = true; }
                  else
                    { intersects_facet[2*d] = true; }
                }
              }
            }
          }
        }

        for (DTYPE d = 0; d < dimension; d++) {

          if (intersects_facet[2*d] != intersects_facet[2*d+1]) {
            if (intersects_facet[2*d]) 
              { vcoord[d] -= offset; }
            else
              { vcoord[d] += offset; }
          }
        }
        IJK::copy_coord(dimension, vcoord, coord+i*dimension);
      }
    }
  }

  /// Position dual isosurface vertices near cube centers.
  /// - More than one vertex can be in a cube.
  /// - C++ STL vector format for array coord[].
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE, 
            typename STYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE0, typename CTYPE1>
  void position_all_dual_isovertices_near_cube_center_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & iso_vlist,
   const CTYPE0 offset,
   std::vector<CTYPE1> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(iso_vlist.size()*dimension);
    position_all_dual_isovertices_near_cube_center_multi
      (scalar_grid, isodual_table, isovalue, iso_vlist, offset, 
       &(coord.front()));
  }

  // **************************************************
  // POSITION LIFTED INTERVAL VOLUME VERTICES
  // **************************************************

  /// Position interval volume vertices which have been lifted
  ///   to one higher dimension.
  /// @param isovalue0 Lower isovalue.
  /// @param isovalue1 Upper isovalue.
  /// @pre isovalue0 < isovalue1.
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISOVAL0_TYPE, typename ISOVAL1_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_ivol_lifted
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISOVAL0_TYPE isovalue0,
   const ISOVAL1_TYPE isovalue1,
   const std::vector<DUAL_ISOV_TYPE> & iso_vlist,
   CTYPE * coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::AXIS_SIZE_TYPE ATYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename GRID_TYPE::SCALAR_TYPE STYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    if (dimension < 1) { return; }

    const VTYPE numv_in_grid_facet_maxd = 
      scalar_grid.AxisIncrement(dimension-1);
    IJK::ARRAY<CTYPE> coord0(dimension);
    IJK::ARRAY<CTYPE> coord1(dimension);
    IJK::ARRAY<CTYPE> coord2(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    for (NTYPE isov = 0; isov < iso_vlist.size(); isov++) {

      const VTYPE icube = iso_vlist[isov].cube_index;

      if (icube < numv_in_grid_facet_maxd) {
        // Apply isovalue1 to (*,*,*,0) grid cubes.
        position_dual_isov_centroid_multi
          (scalar_grid, isodual_table, isovalue1, iso_vlist[isov], cube,
           coord+isov*dimension, coord0.Ptr(), coord1.Ptr(), coord2.Ptr());
      }
      else {
        // Apply isovalue0 to (*,*,*,1) grid cubes.
        position_dual_isov_centroid_multi
          (scalar_grid, isodual_table, isovalue0, iso_vlist[isov], cube,
           coord+isov*dimension, coord0.Ptr(), coord1.Ptr(), coord2.Ptr());

      }
    }
  }


  /// Position interval volume vertices which have been lifted
  ///   to one higher dimension.
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE, 
            typename STYPE0, typename STYPE1,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_ivol_lifted
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE0 isovalue0,
   const STYPE1 isovalue1,
   const std::vector<DUAL_ISOV_TYPE> & iso_vlist,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(iso_vlist.size()*dimension);
    position_all_dual_isovertices_ivol_lifted
      (scalar_grid, isodual_table, isovalue0, isovalue1, iso_vlist,
       &(coord.front()));
  }

}

#endif
