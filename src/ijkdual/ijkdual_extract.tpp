/// @file ijkdual_extract.tpp
/// @brief Template functions for extracting dual isosurface mesh
/// - Version 0.3.1

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


#ifndef _IJKDUAL_EXTRACT_TPP_
#define _IJKDUAL_EXTRACT_TPP_

#include "ijkgrid_macros.h"
#include "ijkdual_types.h"

#include "ijkisopoly.tpp"
#include "ijkmerge.tpp"
#include "ijktime.tpp"


namespace IJKDUAL {

  // ***************************************************
  // EXTRACT ROUTINES
  // ***************************************************

  /*!
   *  @brief Extract isosurface polytopes.
   *  - Returns list representing isosurface polytopes
   *  @param scalar_grid = scalar grid data
   *  @param isovalue = isosurface scalar value
   *  @param iso_poly[] = vector of isosurface polygope vertices
   *    iso_simplices[numv_per_poly*ip+k] = 
   *      cube containing k'th vertex of polytope ip.
   */
  template <typename GTYPE, typename STYPE, typename INFO_TYPE>
  void extract_dual_isopoly
  (const GTYPE & scalar_grid, const STYPE isovalue, 
   std::vector<ISO_VERTEX_INDEX> & iso_poly,
   INFO_TYPE & dualiso_info)
  {
    dualiso_info.time.extract = 0;

    clock_t t0 = clock();

    // initialize output
    iso_poly.clear();

    if (scalar_grid.NumCubeVertices() < 1) { return; }

    IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

      extract_dual_isopoly_around_bipolar_edge
        (scalar_grid, isovalue, iend0, edge_dir, iso_poly);
    }

    clock_t t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }

  /// Extract isosurface polytopes
  /// Returns list representing isosurface polytopes
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polygope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  /// @param dual_edge[] = Array of dual edges sorted by edge direction,
  ///   by grid lines and by order along grid lines.
  ///   If two collinear dual edges e1 and e2 share an endpoint,
  ///     with e1 appearing to the left/below e2, 
  ///     then e1 will be followed by e2 in array dual_edge[].
  template <typename GTYPE, typename STYPE, typename ETYPE, typename INFO_TYPE>
  void extract_dual_isopoly
  (const GTYPE & scalar_grid, const STYPE isovalue, 
   std::vector<ISO_VERTEX_INDEX> & iso_poly,
   std::vector<ETYPE> & dual_edge,
   INFO_TYPE & dualiso_info)
  {
    dualiso_info.time.extract = 0;

    clock_t t0 = clock();

    // initialize output
    iso_poly.clear();

    if (scalar_grid.NumCubeVertices() < 1) { return; }

    IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {
      extract_dual_isopoly_around_bipolar_edge_I
        (scalar_grid, isovalue, iend0, edge_dir, iso_poly, dual_edge);
    }

    clock_t t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }

  /// Extract isosurface polytopes
  /// Returns list representing isosurface polytopes
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polygope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  /// @param facet_vertex[i] = Edge of cube containing iso_poly[i].
  template <typename GTYPE, typename STYPE, typename INFO_TYPE>
  void extract_dual_isopoly
  (const GTYPE & scalar_grid, const STYPE isovalue, 
   std::vector<ISO_VERTEX_INDEX> & iso_poly,
   std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   INFO_TYPE & dualiso_info)
  {
    dualiso_info.time.extract = 0;

    clock_t t0 = clock();

    // initialize output
    iso_poly.clear();
    facet_vertex.clear();

    if (scalar_grid.NumCubeVertices() < 1) { return; }

    IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {

      extract_dual_isopoly_around_bipolar_edge
        (scalar_grid, isovalue, iend0, edge_dir, iso_poly, facet_vertex);
    }

    clock_t t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }


  /// Extract isosurface polytopes
  /// Returns list representing isosurface polytopes
  /// @param scalar_grid = scalar grid data
  /// @param isovalue = isosurface scalar value
  /// @param iso_poly[] = vector of isosurface polygope vertices
  ///   iso_simplices[numv_per_poly*ip+k] = 
  ///     cube containing k'th vertex of polytope ip.
  /// @param facet_vertex[i] = Edge of cube containing iso_poly[i].
  /// @param dual_edge[] = Array of dual edges sorted by edge direction,
  ///   by grid lines and by order along grid lines.
  ///   If two collinear dual edges e1 and e2 share an endpoint,
  ///     with e1 appearing to the left/below e2, 
  ///     then e1 will be followed by e2 in array dual_edge[].
  template <typename GTYPE, typename STYPE, typename ETYPE, typename INFO_TYPE>
  void extract_dual_isopoly
  (const GTYPE & scalar_grid,
   const STYPE isovalue, std::vector<ISO_VERTEX_INDEX> & iso_poly,
   std::vector<FACET_VERTEX_INDEX> & facet_vertex,
   std::vector<ETYPE> & dual_edge,
   INFO_TYPE & dualiso_info)
  {
    dualiso_info.time.extract = 0;

    clock_t t0 = clock();

    // initialize output
    iso_poly.clear();
    facet_vertex.clear();

    if (scalar_grid.NumCubeVertices() < 1) { return; }

    IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, scalar_grid, VERTEX_INDEX) {
      extract_dual_isopoly_around_bipolar_edge_I
        (scalar_grid, isovalue, iend0, edge_dir, iso_poly, facet_vertex,
         dual_edge);
    }

    clock_t t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }


};

#endif
