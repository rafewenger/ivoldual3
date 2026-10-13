/// \file ijktriangulate.txx
/// ijk templates for triangulating polyhedra
/// Version 0.2.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010-2017 Rephael Wenger

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

#ifndef _IJKTRIANGULATE_
#define _IJKTRIANGULATE_

#include "ijkcube.tpp"
#include "ijkmesh.tpp"
#include "ijkmesh_datastruct.tpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


namespace IJK {

  // **************************************************
  //! @name TRIANGULATE POLYGONS
  // **************************************************

  ///@{

  /// Triangulate a polygon by adding diagonals from vertex poly_vert[0].
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - poly_vert[0] is a vertex of all new triangles.
  /// - Add new triangles to vector tri_vert.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1>
  void triangulate_polygon
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, 
   std::vector<VTYPE1> & tri_vert)
  {
    VTYPE0 v0 = poly_vert[0];
    for (NTYPE i = 1; i+1 < num_poly_vert; i++) {
      add_triangle_vertices(v0, poly_vert[i], poly_vert[i+1], tri_vert);
    }
  }


  /// Triangulate a polygon by adding diagonals from vertex poly_vert[0].
  /// - Reverse triangle orientations.
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - poly_vert[0] is a vertex of all new triangles.
  /// - Add new triangles to vector tri_vert.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1>
  void triangulate_polygon_reverse_orient
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, 
   std::vector<VTYPE1> & tri_vert)
  {
    VTYPE0 v0 = poly_vert[0];
    for (NTYPE i = 1; i+1 < num_poly_vert; i++) {
      add_triangle_vertices(v0, poly_vert[i+1], poly_vert[i], tri_vert);
    }
  }


  /// Triangulate a polygon by adding diagonals from vertex poly_vert[0].
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - Add new triangles to vector tri_vert.
  /// @param flag_reverse_orient If true, reverse triangle orientations.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1>
  void triangulate_polygon_from_v0
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const bool flag_reverse_orient, std::vector<VTYPE1> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_polygon_reverse_orient(num_poly_vert, poly_vert, tri_vert);
    }
    else {
      triangulate_polygon(num_poly_vert, poly_vert, tri_vert);
    }
  }


  /// Triangulate a polygon using diagonals from vertex poly_vert[index_v0].
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - poly_vert[index_v0] is a vertex of all new triangles.
  /// - Add new triangles to vector tri_vert.
  template <typename NTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_polygon
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const VTYPE1 index_v0,
   std::vector<VTYPE2> & tri_vert)
  {
    VTYPE0 v0 = poly_vert[index_v0];
    NTYPE i1 = (index_v0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != index_v0) {
      add_triangle_vertices(v0, poly_vert[i1], poly_vert[i2], tri_vert);
      i1 = i2;
      i2 = (i1+1)%num_poly_vert;
    }
  }


  /// Triangulate a polygon using diagonals from vertex poly_vert[index_v0].
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - poly_vert[index_v0] is a vertex of all new triangles.
  /// - Add new triangles to vector tri_vert.
  /// - Reverse triangle orientations.
  template <typename NTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_polygon_reverse_orient
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const VTYPE1 index_v0,
   std::vector<VTYPE2> & tri_vert)
  {
    VTYPE0 v0 = poly_vert[index_v0];
    NTYPE i1 = (index_v0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != index_v0) {
      add_triangle_vertices(v0, poly_vert[i2], poly_vert[i1], tri_vert);
      i1 = i2;
      i2 = (i1+1)%num_poly_vert;
    }
  }


  /// Triangulate a polygon by adding a triangle between each polygon edge
  ///   and vertex w0.
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - Add new triangles to vector tri_vert.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2>
  void triangulate_polygon_with_vertex
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, const VTYPE1 w0,
   std::vector<VTYPE2> & tri_vert)
  {
    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
      NTYPE i1 = (i0+1)%num_poly_vert;
      add_triangle_vertices(w0, poly_vert[i0], poly_vert[i1], tri_vert);
    }
  }


  /// Triangulate a polygon by adding a triangle between each polygon edge
  ///   and vertex w0.
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - Add new triangles to vector tri_vert.
  /// - Reverse triangle orientation.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2>
  void triangulate_polygon_with_vertex_reverse_orient
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, const VTYPE1 w0,
   std::vector<VTYPE2> & tri_vert)
  {
    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
      NTYPE i1 = (i0+1)%num_poly_vert;
      add_triangle_vertices(w0, poly_vert[i1], poly_vert[i0], tri_vert);
    }
  }


  /// Triangulate a polygon by adding a triangle between each polygon edge
  ///   and vertex w0.
  /// - Polygon vertices are listed in clockwise or counter-clockwise order
  ///   around the polygon.
  /// - Add new triangles to vector tri_vert.
  /// @param flag_reverse_orient If true, reverse triangle orientations.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2>
  void triangulate_polygon_with_vertex
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, const VTYPE1 w0,
   const bool flag_reverse_orient, std::vector<VTYPE2> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_polygon_with_vertex_reverse_orient
        (num_poly_vert, poly_vert, w0, tri_vert);
    }
    else {
      triangulate_polygon_with_vertex
        (num_poly_vert, poly_vert, w0, tri_vert);
    }
  }


  /// Triangulate list of polygons.
  template <typename NTYPE0, typename NTYPE1, 
            typename VTYPE0, typename VTYPE1,
            typename ITYPE>
  void triangulate_polygon_list
  (const NTYPE0 * num_poly_vert, const VTYPE0 * poly_vert,
   const ITYPE * first_poly_vert, const NTYPE1 num_poly,
   std::vector<VTYPE1> & tri_vert)
  {
    for (NTYPE1 ipoly = 0; ipoly < num_poly; ipoly++) {
      triangulate_polygon
        (num_poly_vert[ipoly], poly_vert+first_poly_vert[ipoly], tri_vert);
    }
  }

  /// Triangulate list of polygons.
  /// - C++ STL vector format for num_poly_vert[], poly_vert[], 
  ///   first_poly_vert[].
  template <typename NTYPE, typename VTYPE0, typename VTYPE1,
            typename ITYPE>
  void triangulate_polygon_list
  (const std::vector<NTYPE> & num_poly_vert, 
   const std::vector<VTYPE0> & poly_vert,
   const std::vector<ITYPE> & first_poly_vert,
   std::vector<VTYPE1> & tri_vert)
  {
    triangulate_polygon_list
      (IJK::vector2pointer(num_poly_vert), IJK::vector2pointer(poly_vert),
       IJK::vector2pointer(first_poly_vert), num_poly_vert.size(), tri_vert);
  }

  ///@}


  // **************************************************
  //! @name TRIANGULATE QUADRILATERALS
  // **************************************************

  ///@{

  /// Triangulate a set of quadrilaterals.
  /// - Quadrilateral vertices are listed in clockwise or counter-clockwise 
  ///   order around the polygon.
  /// - Add new triangles to vector tri_vert.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1>
  void triangulate_quad
  (const VTYPE0 * quad_vert, const NTYPE num_quad,
   std::vector<VTYPE1> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {

      NTYPE k = iquad*NUM_VERT_PER_QUAD;
      triangulate_polygon(NUM_VERT_PER_QUAD, quad_vert+k, tri_vert);
    }
  }

  /// Triangulate a set of quadrilaterals.
  /// - C++ STL vector format for quad_vert.
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_quad
  (const std::vector<VTYPE0> quad_vert, std::vector<VTYPE1> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;

    SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    triangulate_quad(IJK::vector2pointer(quad_vert), num_quad, tri_vert);
  }


  /// Triangulate a set of quadrilaterals by adding a triangle between
  ///   each quad edge and a new vertex for each quad.
  /// - Quadrilateral vertices are listed in clockwise or counter-clockwise
  ///   order around the polygon.
  /// - Quad i is triangulated using vertex (first_vertex+i).
  /// - Add new triangles to vector tri_vert.
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2>
  void triangulate_quad_with_vertices
  (const VTYPE0 * quad_vert, const NTYPE num_quad, 
   const VTYPE1 first_vertex, std::vector<VTYPE2> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {

      NTYPE k = iquad*NUM_VERT_PER_QUAD;
      triangulate_polygon_with_vertex
        (NUM_VERT_PER_QUAD, quad_vert+k, first_vertex+iquad, tri_vert);
    }
  }

  /// Triangulate a set of quadrilaterals by adding a triangle between
  ///   each quad edge and a new vertex for each quad.
  /// - C++ STL vector format for quad_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_with_vertices
  (const std::vector<VTYPE0> & quad_vert,
   const VTYPE1 first_vertex, std::vector<VTYPE2> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;

    SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    triangulate_quad_with_vertices
      (IJK::vector2pointer(quad_vert), num_quad, first_vertex, tri_vert);
  }


  /// Triangulate a set of quadrilaterals with diagonal (0,2) or (1,3).
  /// - Quadrilateral vertices are listed in clockwise or counter-clockwise
  ///   order around the polygon.
  /// - Add new triangles to vector tri_vert.
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_quad_using_diagonal
  (const VTYPE0 * quad_vert, const bool flag_diag02,
   std::vector<VTYPE1> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);

    if (flag_diag02) {
      triangulate_polygon(NUM_VERT_PER_QUAD, quad_vert, 0, tri_vert);
    }
    else {
      triangulate_polygon(NUM_VERT_PER_QUAD, quad_vert, 1, tri_vert);
    }
  }

  /// Triangulate a set of quadrilaterals with diagonal (0,2) or (1,3).
  /// - Quadrilateral vertices are listed in clockwise or counter-clockwise
  ///   order around the polygon.
  /// - Add new triangles to vector tri_vert.
  /// - Reverse triangle_orientation
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_quad_using_diagonal_reverse_orient
  (const VTYPE0 * quad_vert, const bool flag_diag02,
   std::vector<VTYPE1> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);

    if (flag_diag02) {
      triangulate_polygon_reverse_orient
        (NUM_VERT_PER_QUAD, quad_vert, 0, tri_vert);
    }
    else {
      triangulate_polygon_reverse_orient
        (NUM_VERT_PER_QUAD, quad_vert, 1, tri_vert);
    }
  }


  /// Triangulate a set of quadrilaterals with diagonal (0,2) or (1,3).
  /// - Quadrilateral vertices are listed in clockwise or counter-clockwise
  ///   order around the polygon.
  /// - Add new triangles to vector tri_vert.
  /// @param flag_reverse_orient If true, reverse triangle orientation.
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_quad_using_diagonal
  (const VTYPE0 * quad_vert, const bool flag_diag02,
   const bool flag_reverse_orient, std::vector<VTYPE1> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_quad_using_diagonal_reverse_orient
        (quad_vert, flag_diag02, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal(quad_vert, flag_diag02, tri_vert);
    }
  }


  ///@}


  // **************************************************
  //! @name TRIANGULATE PENTAGONS
  // **************************************************

  ///@{

  /// Triangulate a pentagon by adding triangles
  ///   (v0,v1,v2), (v0, v2, v3) and (v0, v3, v4).
  /// - Pentagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the pentagon.
  /// - Add new triangles to vector tri_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPEB>
  void triangulate_pentagon
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v0, v1, v2, tri_vert);
    add_triangle_vertices(v0, v2, v3, tri_vert);
    add_triangle_vertices(v0, v3, v4, tri_vert);
  }

  /// Triangulate a pentagon by adding triangles
  ///   (v0,v1,v2), (v0, v2, v3) and (v0, v3, v4).
  /// - Version allowing reverse orientation.
  /// @param flag_reverse_orient If true, reverse orientation.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPEB>
  void triangulate_pentagon
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const bool flag_reverse_orient,
   std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_pentagon(v0, v4, v3, v2, v1, tri_vert);
    }
    else {
      triangulate_pentagon(v0, v1, v2, v3, v4, tri_vert);
    }
  }


  /// Triangulate a pentagon by adding triangles incident on pentagon_vert[ivX].
  /// - Pentagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the pentagon.
  /// - Add new triangles to vector tri_vert.
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_pentagon
  (const VTYPE pentagon_vert[], const ITYPE ivX,
   std::vector<VTYPEB> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_PENTAGON(5);
    const ITYPE i1 = (ivX+1)%NUM_VERT_PER_PENTAGON;
    const ITYPE i2 = (ivX+2)%NUM_VERT_PER_PENTAGON;
    const ITYPE i3 = (ivX+3)%NUM_VERT_PER_PENTAGON;
    const ITYPE i4 = (ivX+4)%NUM_VERT_PER_PENTAGON;

    triangulate_pentagon
      (pentagon_vert[ivX], pentagon_vert[i1], pentagon_vert[i2],
       pentagon_vert[i3], pentagon_vert[i4], tri_vert);
  }

  /// Triangulate a pentagon.
  /// - Triangles are specified by ears.
  /// @param ear0 Index of first ear in triangulation. In range [0..5].
  /// @param ear1 Index of second ear in triangulation. In range [0..5].
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, 
            typename EAR_TYPE0, typename EAR_TYPE1,
            typename VTYPEB>
  void triangulate_pentagon_by_ears
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const EAR_TYPE0 ear0, const EAR_TYPE1 ear1,
   std::vector<VTYPEB> & tri_vert)
  {
    switch(ear0) {

    case 0:
      if (ear1 == 1 || ear1 == 3) 
        { triangulate_pentagon(v4, v0, v1, v2, v3, tri_vert); }
      else
        { triangulate_pentagon(v1, v2, v3, v4, v0, tri_vert); }
      break;

    case 1:
      if (ear1 == 2 || ear1 == 4) 
        { triangulate_pentagon(v0, v1, v2, v3, v4, tri_vert); }
      else
        { triangulate_pentagon(v2, v3, v4, v0, v1, tri_vert); }
      break;

    case 2:
      if (ear1 == 0 || ear1 == 3) 
        { triangulate_pentagon(v1, v2, v3, v4, v0, tri_vert); }
      else
        { triangulate_pentagon(v3, v4, v0, v1, v2, tri_vert); }
      break;

    case 3:
      if (ear1 == 0 || ear1 == 2) 
        { triangulate_pentagon(v4, v0, v1, v2, v3, tri_vert); }
      else
        { triangulate_pentagon(v2, v3, v4, v0, v1, tri_vert); }
      break;

    case 4:
    default:
      if (ear1 == 0 || ear1 == 2) 
        { triangulate_pentagon(v3, v4, v0, v1, v2, tri_vert); }
      else
        { triangulate_pentagon(v0, v1, v2, v3, v4, tri_vert); }
      break;
    }
  }


  /// Triangulate a pentagon by adding triangles
  ///   (w0,v0,v1), (w0,v1,v2), (w0,v2,v3), (w0,v3,v4) and (w0, v4, v0).
  /// - Pentagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the pentagon.
  /// - Add new triangles to vector tri_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename WTYPE, typename VTYPEB>
  void triangulate_pentagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const WTYPE w0,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v0, v1, tri_vert);
    add_triangle_vertices(w0, v1, v2, tri_vert);
    add_triangle_vertices(w0, v2, v3, tri_vert);
    add_triangle_vertices(w0, v3, v4, tri_vert);
    add_triangle_vertices(w0, v4, v0, tri_vert);
  }


  /// Triangulate a pentagon by adding triangles
  ///   (w0,v1,v0), (w0,v2,v1), (w0,v3,v2), (w0,v4,v3) and (w0, v0, v4).
  /// - Pentagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the pentagon.
  /// - Add new triangles to vector tri_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename WTYPE, typename VTYPEB>
  void triangulate_pentagon_with_vertex_reverse_orient
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const WTYPE w0,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v1, v0, tri_vert);
    add_triangle_vertices(w0, v2, v1, tri_vert);
    add_triangle_vertices(w0, v3, v2, tri_vert);
    add_triangle_vertices(w0, v4, v3, tri_vert);
    add_triangle_vertices(w0, v0, v4, tri_vert);
  }


  /// Triangulate a pentagon by adding a triangle between each pentagon edge
  ///   and vertex w0.
  /// - Pentagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the pentagon.
  /// - Add new triangles to vector tri_vert.
  /// @param flag_reverse_orient If true, reverse triangle orientations.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename WTYPE, typename VTYPEB>
  void triangulate_pentagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const WTYPE w0,
   const bool flag_reverse_orient, std::vector<VTYPEB> & tri_vert)
  {

    if (flag_reverse_orient) {
      triangulate_pentagon_with_vertex_reverse_orient
        (v0, v1, v2, v3, v4, w0, tri_vert);
    }
    else {
      triangulate_pentagon_with_vertex
        (v0, v1, v2, v3, v4, w0, tri_vert);
    }
  }

  ///@}



  // **************************************************
  //! @name TRIANGULATE HEXAGONS
  // **************************************************

  ///@{

  /// Triangulate a hexagon by adding triangle (iv0,iv2,iv4).
  /// - Hexagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the hexagon.
  /// - Add new triangles to vector tri_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon_using_triangle_024
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v0, v1, v2, tri_vert);
    add_triangle_vertices(v2, v3, v4, tri_vert);
    add_triangle_vertices(v4, v5, v0, tri_vert);
    add_triangle_vertices(v0, v2, v4, tri_vert);
  }

  /// Triangulate a hexagon by adding triangle (iv4,iv2,iv0).
  /// - Triangles have reverse orientation of hexagon.
  /// - Hexagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the hexagon.
  /// - Add new triangles to vector tri_vert.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon_using_triangle_024_reverse_orient
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v2, v1, v0, tri_vert);
    add_triangle_vertices(v4, v3, v2, tri_vert);
    add_triangle_vertices(v0, v5, v4, tri_vert);
    add_triangle_vertices(v4, v2, v0, tri_vert);
  }

  /// Triangulate a hexagon by adding triangle (iv0,iv2,iv4).
  /// - Version with flag_reverse_orient.
  /// @param flag_reverse_orient 
  ///   If true, reverse orientation in creating triangles from the hexagon.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon_using_triangle_024
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const bool flag_reverse_orient,
   std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_hexagon_using_triangle_024_reverse_orient
        (v0, v1, v2, v3, v4, v5, tri_vert);
    }
    else {
      triangulate_hexagon_using_triangle_024
        (v0, v1, v2, v3, v4, v5, tri_vert);
    }
  }

  /// Triangulate a hexagon.
  /// - Triangles are specified by ears.
  /// - Hexagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the hexagon.
  /// - Add new triangles to vector tri_vert.
  /// @param ear0 Index of first ear in triangulation. In range [0..5].
  /// @param ear1 Index of second ear in triangulation. In range [0..5].
  /// @param ear2 Index of third ear in triangulation. In range [0..5].
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename EAR_TYPE0, typename EAR_TYPE1, typename EAR_TYPE2,
            typename VTYPEB>
  void triangulate_hexagon_by_ears
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const EAR_TYPE0 ear0, const EAR_TYPE1 ear1, const EAR_TYPE2 ear2,
   std::vector<VTYPEB> & tri_vert)
  {
    EAR_TYPE1 ear1B;
    EAR_TYPE2 ear2B;

    if (ear1 > ear0) { ear1B = ear1-1; }
    else { ear1B = ear1; }
    if (ear2 > ear0) { ear2B = ear2-1; }
    else { ear2B = ear2; }

    switch(ear0) {

    case 0:
      add_triangle_vertices(v5, v0, v1, tri_vert);
      triangulate_pentagon_by_ears
        (v1, v2, v3, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 1:
      add_triangle_vertices(v0, v1, v2, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v2, v3, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 2:
      add_triangle_vertices(v1, v2, v3, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v3, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 3:
      add_triangle_vertices(v2, v3, v4, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v2, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 4:
      add_triangle_vertices(v3, v4, v5, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v2, v3, v5, ear1B, ear2B, tri_vert);
      break;

    case 5:
    default:
      add_triangle_vertices(v4, v5, v0, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v2, v3, v4, ear1B, ear2B, tri_vert);
      break;
    }
  }

  /// Triangulate a hexagon.
  /// - Version with array of hexagon vertices.
  /// - Hexagon vertices are listed in clockwise or counter-clockwise 
  ///   order around the hexagon.
  template <typename VTYPE,
            typename EAR_TYPE0, typename EAR_TYPE1, typename EAR_TYPE2,
            typename VTYPEB>
  void triangulate_hexagon_by_ears
  (const VTYPE hex_vert[],
   const EAR_TYPE0 ear0, const EAR_TYPE1 ear1, const EAR_TYPE2 ear2,
   std::vector<VTYPEB> & tri_vert)
  {
    triangulate_hexagon_by_ears
      (hex_vert[0], hex_vert[1], hex_vert[2], hex_vert[3], hex_vert[4],
       hex_vert[5], ear0, ear1, ear2, tri_vert);
  }

  ///@}


  // **************************************************
  // @name CLASS HEX_TRIANGULATION_INFO
  // **************************************************

  ///@{

  /// Information on hex triangulation.
  template <typename ANCHOR_TYPE, typename FLAGS_TYPE>
  class HEX_TRIANGULATION_INFO {
  public:

    /// Hex vertex (0,1,...,7) which is the triangulation anchor.
    /// All triangulation tetrahedra are incident on the triangulation_anchor.
    ANCHOR_TYPE triangulation_anchor;

    /// Flags indicating triangulation of hex facet diagonals.
    FLAGS_TYPE triangulation_flags;

    HEX_TRIANGULATION_INFO() 
    {
      triangulation_flags = 0;
    }

    template <typename FTYPE>
    bool TriangulationFlag(const FTYPE jfacet)
    {
      const FLAGS_TYPE mask = (FLAGS_TYPE(1) << jfacet);
      return(!((mask & triangulation_flags) == 0));
    }

    template <typename FTYPE>
    void FlipTriangulationFlag(const FTYPE jfacet)
    {
      const FLAGS_TYPE mask = (FLAGS_TYPE(1) << jfacet);
      triangulation_flags = (triangulation_flags ^ mask);
    }
  };

  ///@}



  // **************************************************
  //! @name TRIANGULATE SQUARE PYRAMIDS
  // **************************************************

  ///@{

  /// Triangulate square pyramid adding diagonal (v0,v2).
  /// - Vertices of pyramid base are listed in clockwise/counter-clockwise
  ///   order around the square base.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_square_pyramid_diagonal02
  (const VTYPE0 apex, const VTYPE1 v0, const VTYPE1 v1,
   const VTYPE1 v2, const VTYPE1 v3,
   std::vector<VTYPE2> & tet_vert_list)
  {
    add_tetrahedron_vertices(apex, v0, v2, v3, tet_vert_list);
    add_tetrahedron_vertices(apex, v2, v0, v1, tet_vert_list);
  }


  /// Triangulate square pyramid adding diagonal (v1,v3).
  /// - Vertices of pyramid base are listed in clockwise/counter-clockwise
  ///   order around the square base.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_square_pyramid_diagonal13
  (const VTYPE0 apex, const VTYPE1 v0, const VTYPE1 v1,
   const VTYPE1 v2, const VTYPE1 v3,
   std::vector<VTYPE2> & tet_vert_list)
  {
    add_tetrahedron_vertices(apex, v1, v3, v0, tet_vert_list);
    add_tetrahedron_vertices(apex, v3, v1, v2, tet_vert_list);
  }


  /// Triangulate square pyramid.  Use facet_triangulation_bits
  ///   to determine triangulation.
  /// - Vertices of pyramid base are listed in clockwise/counter-clockwise
  ///   order around the square base.
  /// - If ((facet_triangulation_bits & mask) == 0), then
  ///   triangulate facet with diagonal (v0,v2).
  /// - Otherwise, triangulate facet with diagonal (v1,v2).
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename BITS_TYPE>
  void triangulate_square_pyramid_diagonal_bits
  (const VTYPE0 apex, const VTYPE1 v0, const VTYPE1 v1,
   const VTYPE1 v2, const VTYPE1 v3,
   const BITS_TYPE facet_triangulation_bits,
   const BITS_TYPE mask,
   std::vector<VTYPE2> & tet_vert_list)
  {
    if ((facet_triangulation_bits & mask) == 0) {
      triangulate_square_pyramid_diagonal02
        (apex, v0, v1, v2, v3, tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal13
        (apex, v0, v1, v2, v3, tet_vert_list);
    }
  }

  ///@}


  // **************************************************
  //! @name TRIANGULATE HEXAHEDRA
  // **************************************************

  ///@{

  /// Triangulate a hexahedron with tetrahedra all incident
  ///   on the diagonal (hex_vert[0], hex_vert[7])
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_hexahedron_diagonal07
  (const VTYPE0 hex_vert[], std::vector<VTYPE1> & tet_vert_list)
  {
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[3], hex_vert[2], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[2], hex_vert[6], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[6], hex_vert[4], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[4], hex_vert[5], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[5], hex_vert[1], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[1], hex_vert[3], tet_vert_list);
  }


  /// Triangulate a hexahedron using a triangulation whose tetrahedra
  ///   all contain vertex v0.  Triangulation of the facets not
  ///   containing v0 is determined by the facet_triangulation_bits.
  /// @param hex_vert[] Array of 8 vertices of the hexahedron.
  /// @param v0 Vertex index between 0 and 7, inclusive.
  ///    All tetrahedra are incident on v0.
  /// @param facet_triangulation_bits Bits indicating the triangulation
  ///   of facets not incident on v0.
  ///   - If bit i is 0, then the triangulation diagonal of facet i
  ///     is either incident on vertex 0 or on vertex 7.
  ///   - If bit i is 1, then the triangulation diagonal of facet i
  ///     is not incident on vertex 0 or on vertex 7.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename BITS_TYPE>
  void triangulate_hexahedron_anchored
  (const VTYPE0 hex_vert[], const VTYPE1 v0,
   const BITS_TYPE facet_triangulation_bits,
   std::vector<VTYPE2> & tet_vert_list)
  {
    const int NUM_CUBE_FACETS(6);
    const BITS_TYPE mask[NUM_CUBE_FACETS] = 
      { BITS_TYPE(1), BITS_TYPE(2), BITS_TYPE(4),
        BITS_TYPE(8), BITS_TYPE(16), BITS_TYPE(32) };

    if (((v0 & mask[0]) == 0)) {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[1], hex_vert[3], hex_vert[7], hex_vert[5],
         facet_triangulation_bits, mask[3], tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[0], hex_vert[4], hex_vert[6], hex_vert[2],
         facet_triangulation_bits, mask[0], tet_vert_list);
    }

    if (((v0 & mask[1]) == 0)) {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[2], hex_vert[6], hex_vert[7], hex_vert[3],
         facet_triangulation_bits, mask[4], tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[0], hex_vert[1], hex_vert[5], hex_vert[4],
         facet_triangulation_bits, mask[1], tet_vert_list);
    }

    if (((v0 & mask[2]) == 0)) {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[4], hex_vert[5], hex_vert[7], hex_vert[6],
         facet_triangulation_bits, mask[5], tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[0], hex_vert[2], hex_vert[3], hex_vert[1],
         facet_triangulation_bits, mask[2], tet_vert_list);
    }

  }

  /// Triangulate all hexahedron in a hex mesh using triangulation anchors
  ///   and triangulation flags.<br>
  /// Return list of tetrahedra vertices.
  /// @param hex_mesh Hexahedral mesh. 
  ///    Contains member fields poly_data[].triangulation_anchor and
  ///    poly_data[].triangulation_flags.
  template <typename HEX_MESH_TYPE, typename VTYPE>
  void triangulate_anchored_hex_mesh
  (const HEX_MESH_TYPE & hex_mesh,
   std::vector<VTYPE> & tet_vert_list)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE ipoly = 0; ipoly < hex_mesh.NumPoly(); ipoly++) {
      triangulate_hexahedron_anchored
        (hex_mesh.VertexList(ipoly), 
         hex_mesh.poly_data[ipoly].triangulation_anchor,
         hex_mesh.poly_data[ipoly].triangulation_flags,
         tet_vert_list);
    }
  }


  /// Set the triangulation anchors for a hexahedral mesh.
  /// - No facet contains two distinct triangulation anchors 
  ///        (associated with different hexahedra.)
  template <typename VERTEX_POLY_INCIDENCE_TYPE,
            typename HEX_MESH_TYPE>
  void set_hex_triangulation_anchors
  (const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
   HEX_MESH_TYPE & hex_mesh)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename HEX_MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    std::vector<bool> flag_anchor_set(hex_mesh.NumPoly(), false);

    for (NTYPE ipoly = 0; ipoly < hex_mesh.NumPoly(); ipoly++) {
      if (!flag_anchor_set[ipoly]) {
        hex_mesh.poly_data[ipoly].triangulation_anchor = 0;
        flag_anchor_set[ipoly] = true;
        const VTYPE iv0 = hex_mesh.Vertex(ipoly, 0);

        for (NTYPE j = 0; j < vertex_poly_incidence.NumIncidentPoly(iv0); j++) {
          const NTYPE jpoly = vertex_poly_incidence.IncidentPoly(iv0, j);

          if (flag_anchor_set[jpoly]) { continue; }

          NTYPE kloc;
          if (does_list_contain
              (hex_mesh.VertexList(jpoly), hex_mesh.NumPolyVert(jpoly), 
               iv0, kloc)) {
            hex_mesh.poly_data[jpoly].triangulation_anchor = kloc;
            flag_anchor_set[jpoly] = true;
          }
        }
      }
    }

  }


  // Set triangulation flags for hexahedra in hex_mesh
  // @pre Triangulation anchors are already set.
  // @pre No facet contains two distinct triangulation anchors 
  //        (associated with different hexahedra.)
  template <typename HEX_MESH_TYPE>
  void set_hex_triangulation_flags(HEX_MESH_TYPE & hex_mesh)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename HEX_MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE DIM3(3);
    const IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube_face_info(DIM3);
    const NTYPE NUM_VERT_PER_TETRAHEDRON(4);
    const NTYPE NUM_CUBE_VERT(8);
    const NTYPE NUM_CUBE_FACETS(6);
    const char DEFAULT_TRIANGULATION_FLAG[NUM_CUBE_VERT]
      = { char(0), char(06), char(05), char(030),
          char(03), char(050), char(060), char(0) };
    const bool VERTEX_FACET_FLAGS[NUM_CUBE_VERT][NUM_CUBE_FACETS]
      = { { true, true, true, false, false, false},
          { false, true, true, true, false, false},
          { true, false, true, false, true, false},
          { false, false, true, true, true, false},
          { true, true, false, false, false, true},
          { false, true, false, true, false, true},
          { true, false, false, false, true, true},
          { false, false, false, true, true, true} };

    for (NTYPE ipoly = 0; ipoly < hex_mesh.NumPoly(); ipoly++) {
      const VTYPE iv = hex_mesh.poly_data[ipoly].triangulation_anchor;
      hex_mesh.poly_data[ipoly].triangulation_flags =
        DEFAULT_TRIANGULATION_FLAG[iv];
    }
  
    IJK::HEXMESH_ADJACENCY<NTYPE,NTYPE,NTYPE> hexmesh_adjacency;
    VTYPE diagonal0_endpoint[2];
    VTYPE diagonal1_endpoint[2];

    hexmesh_adjacency.SetFromMeshOfCubes(hex_mesh, cube_face_info);

    for (NTYPE ipoly0 = 0; ipoly0 < hex_mesh.NumPoly(); ipoly0++) {

      for (NTYPE jfacet0 = 0; jfacet0 < NUM_CUBE_FACETS; jfacet0++) {
        if (hexmesh_adjacency.IsAdjacent(ipoly0, jfacet0)) {
          const NTYPE ipoly1 = hexmesh_adjacency.AdjacentPoly(ipoly0, jfacet0);
          const NTYPE jfacet1 = 
            hexmesh_adjacency.AdjacentPolyFacet(ipoly0, jfacet0);
          const bool tri0_flag = 
            hex_mesh.poly_data[ipoly0].TriangulationFlag(jfacet0);
          const bool tri1_flag = 
            hex_mesh.poly_data[ipoly1].TriangulationFlag(jfacet1);

          get_hexahedron_facet_diagonal
            (hex_mesh.VertexList(ipoly0), jfacet0, tri0_flag, 
             diagonal0_endpoint);
          get_hexahedron_facet_diagonal
            (hex_mesh.VertexList(ipoly1), jfacet1, tri1_flag, 
             diagonal1_endpoint);

          if (diagonal0_endpoint[0] != diagonal1_endpoint[0]) {
            const VTYPE iv1 = hex_mesh.poly_data[ipoly1].triangulation_anchor;
            if (VERTEX_FACET_FLAGS[iv1][jfacet1]) {
              hex_mesh.poly_data[ipoly0].FlipTriangulationFlag(jfacet0);
            }
            else {
              hex_mesh.poly_data[ipoly1].FlipTriangulationFlag(jfacet1);
            }
          }
        }
      }
    }

  }

  /// Triangulate a hexahedral mesh.<br>
  /// Return list of tetrahedra vertices.
  /// - Version which takes a hex mesh with member fields storing
  ///   triangulation anchors and triangulation flags.
  template <typename HEX_MESH_WITH_TRI_DATA, typename VTYPE>
  void triangulate_hex_meshX
  (HEX_MESH_WITH_TRI_DATA & hex_mesh, std::vector<VTYPE> & tet_vert_list)
  {
    typedef typename HEX_MESH_WITH_TRI_DATA::NUMBER_TYPE NTYPE;

    IJK::VERTEX_POLY_INCIDENCE<NTYPE,NTYPE> vertex_poly_incidence;

    vertex_poly_incidence.Set(hex_mesh);

    set_hex_triangulation_anchors(vertex_poly_incidence, hex_mesh);
    set_hex_triangulation_flags(hex_mesh);

    triangulate_anchored_hex_mesh(hex_mesh, tet_vert_list);
  }


  /// Triangulate a hexahedral mesh.<br>
  /// Return list of tetrahedra vertices.
  template <typename HEX_MESH_TYPE, typename VTYPE>
  void triangulate_hex_mesh
  (const HEX_MESH_TYPE & hex_mesh, std::vector<VTYPE> & tet_vert_list)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NTYPE;

    IJK::VERTEX_POLY_INCIDENCE<NTYPE,NTYPE> vertex_poly_incidence;
    IJK::POLYMESH_DATA<VTYPE,NTYPE,
      IJK::HEX_TRIANGULATION_INFO<char,char>> hex_data;

    hex_data.Set(hex_mesh);
    triangulate_hex_meshX(hex_data, tet_vert_list);
  }

  ///@}

}

#endif
