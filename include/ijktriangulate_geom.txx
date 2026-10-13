/// \file ijktriangulate_geom.txx
/// ijk templates for triangulating polytopes.
/// - Includes code based on geomtry.
/// - Version 0.2.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2018 Rephael Wenger

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

#ifndef _IJKTRIANGULATE_GEOM_
#define _IJKTRIANGULATE_GEOM_

#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijktriangulate.txx"


namespace IJK {

  // **************************************************
  //! @name COMPUTE COS MIN ANGLE
  // **************************************************

  ///@{

  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a quadrilateral.
   *  - Quadrilateral is split by diagonal (coord0[],coord2[]) 
   *    into two triangles.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename COS_TYPE, typename MTYPE>
  void compute_cos_min_quad_tri02_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA, cosB;
    bool flagA, flagB;

    compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coord2, max_small_magnitude, cosA, flagA);
    compute_cos_min_triangle_angle
      (dimension, coord0, coord3, coord2, max_small_magnitude, cosB, flagB);

    if (flagA || flagB) {
      cos_min_angle = 0;
      flag_zero = true;
    }
    else {
      if (cosA > cosB) { cos_min_angle = cosA; }
      else { cos_min_angle = cosB; }
      flag_zero = false;
    }
  }


  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a quadrilateral.
   *  - Quadrilateral is split by diagonal (coord1[],coord3[]) 
   *    into two triangles.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle have two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename MTYPE>
  void compute_cos_min_quad_tri13_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, CTYPE4 & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_quad_tri02_angle
      (dimension, coord1, coord2, coord3, coord0, max_small_magnitude,
       cos_min_angle, flag_zero);
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
   *  - Quadrilateral is either split by diagonal (coord0[], coord2[])
   *    or by diagonal (coord1[],coord3[]) into two triangles.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_diag02 If true, triangulation with min angle
   *       has diagonal (coord0[], coord2[]).  Otherwise, triangulation 
   *       with min angle has diagonal (coord1[], coord2[]).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    with min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename COS_TYPE, typename MTYPE>
  void compute_cos_min_quad_tri02_tri13_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   bool & flag_diag02, bool & flag_zero)
  {
    COS_TYPE cos_min_diag02, cos_min_diag13;
    bool flag_zero_diag02, flag_zero_diag13;

    compute_cos_min_quad_tri02_angle
      (dimension, coord0, coord1, coord2, coord3, max_small_magnitude, 
       cos_min_diag02, flag_zero_diag02);
    compute_cos_min_quad_tri13_angle
      (dimension, coord0, coord1, coord2, coord3, max_small_magnitude, 
       cos_min_diag13, flag_zero_diag13);

    int index_selected;
    select_min
      (cos_min_diag02, flag_zero_diag02, 
       cos_min_diag13, flag_zero_diag13, 
       cos_min_angle, index_selected, flag_zero);

    if (index_selected == 0) { flag_diag02 = true; }
    else { flag_diag02 = false; }
  }


  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a quadrilateral.
   *  - Quadrilateral is split into four triangles formed from coordX[]
   *    and each of the four quadrilateral edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the four triangles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coordX[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPEX, typename COS_TYPE,
            typename MTYPE>
  void compute_cos_min_quad_tri4_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPEX * coordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_QUAD = 4;
    COS_TYPE cos[NUM_VERT_PER_QUAD];
    bool flag_tri_zero[NUM_VERT_PER_QUAD];

    compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coordX, max_small_magnitude, 
       cos[0], flag_tri_zero[0]);
    compute_cos_min_triangle_angle
      (dimension, coord1, coord2, coordX, max_small_magnitude, 
       cos[1], flag_tri_zero[1]);
    compute_cos_min_triangle_angle
      (dimension, coord2, coord3, coordX, max_small_magnitude, 
       cos[2], flag_tri_zero[2]);
    compute_cos_min_triangle_angle
      (dimension, coord3, coord0, coordX, max_small_magnitude, 
       cos[3], flag_tri_zero[3]);

    bool is_cos_min_angle_set = false;
    for (int i = 0; i < NUM_VERT_PER_QUAD; i++) {

      if (flag_tri_zero[i]) {
        flag_zero = true;
        cos_min_angle = 0;
        return;
      }

      if (!is_cos_min_angle_set || cos[i] > cos_min_angle) {
        cos_min_angle = cos[i];
        is_cos_min_angle_set = true;
      }
    }

    flag_zero = false;
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals or
   *    the triangulation into four triangles.
   *  - Quadrilateral is either split by diagonal (coord0[], coord2[])
   *    or by diagonal (coord1[],coord3[]) into two triangles or
   *    split into 4 triangles all incident on coordX[].
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coordX[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_tri4 If true, split into 4 triangles.
   *  @param[out] flag_diag02 If true, triangulation with min angle
   *       has diagonal (coord0[], coord2[]).  Otherwise, triangulation 
   *       with min angle has diagonal (coord1[], coord2[]).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    with min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPEX, typename COS_TYPE,
            typename MTYPE>
  void compute_cos_min_quad_tri02_tri13_tri4_angle
  (const DTYPE dimension, const CTYPE0 * vcoord0, const CTYPE1 * vcoord1, 
   const CTYPE2 * vcoord2, const CTYPE3 * vcoord3, const CTYPEX * vcoordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   bool & flag_tri4, bool & flag_diag02, bool & flag_zero)
  {
    COS_TYPE cos_min_diag, cos_min_tri4;
    bool flag_zero_diag, flag_zero_tri4;

    compute_cos_min_quad_tri02_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
       cos_min_diag, flag_diag02, flag_zero_diag);
    compute_cos_min_quad_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoordX,
       max_small_magnitude, cos_min_tri4, flag_zero_tri4);

    int index_selected;
    select_min
      (cos_min_tri4, flag_zero_tri4, 
       cos_min_diag, flag_zero_diag, 
       cos_min_angle, index_selected, flag_zero);

    if (index_selected == 0) { flag_tri4 = true; }
    else { flag_tri4 = false; }
  }


  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a polygon from polygon vertex iv0.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param iv0 Triangulation is from iv0.
   *  @pre iv0 is in range [0,num_poly_vert-1].
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename IVTYPE,typename COS_TYPE>
  void compute_cos_min_triangulation_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   const IVTYPE index_v0, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    NTYPE i1 = (index_v0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != index_v0) {
      
      const CTYPE * v0_coord = vcoord+poly_vert[index_v0]*dimension;
      const CTYPE * v1_coord = vcoord+poly_vert[i1]*dimension;
      const CTYPE * v2_coord = vcoord+poly_vert[i2]*dimension;

      compute_cos_min_triangle_angle
        (dimension, v0_coord, v1_coord, v2_coord, max_small_magnitude, 
         cosA, flagA);

      if (flagA) { 
        cos_min_angle = 0;
        flag_zero = true;
        return;
      }
      else {
        if (cosA > cos_min_angle) 
          { cos_min_angle = cosA; }
      }
      i1 = i2;
      i2 = (i1+1)%num_poly_vert;
    }

    return;
  }


  ///  Compute the cosine of the smallest triangle angle
  ///   in the triangulation of a polygon from polygon vertex iv0.
  /// - C++ STL vector format for vcoord.
  template <typename DTYPE, typename NTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename IVTYPE,typename COS_TYPE>
  void compute_cos_min_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vcoord, 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   const IVTYPE index_v0, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), num_poly_vert,
       poly_vert, max_small_magnitude, index_v0, cos_min_angle, flag_zero);
  }

  /*!
   *  Compute the cosine of the smallest triangle angle in the 
   *    triangulation of a polygon from an additional vertex at coordX[].
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param pcoord[] Triangulation is from pcoord[]
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  /* NOT TESTED....
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE0, typename CTYPE1, typename VTYPE,
            typename MTYPE, typename IVTYPE,typename COS_TYPE,
            typename ITYPE>
  void compute_cos_min_triangulation_angle_add_vertex
  (const DTYPE dimension, const CTYPE0 vcoord[], const NTYPE num_poly_vert, 
   const VTYPE poly_vert[], const CTYPE1 coordX[], 
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero, 
   ITYPE & index_of_triangle_with_min_cos)
  {
    CTYPE0 * vcoord0[], vcoord1[];
    COS_TYPE cos_min_angle_i;
    bool flag_zero_i;

    cos_min_angle = -1;
    flag_zero = false;
    index_of_triangle_with_min_cos = 0;

    if (num_poly_vert < 3) { 
      flag_zero = true;
      return; 
    }

    vcoord0 = vcoord + poly_vert[0]*dimension;
    vcoord1 = vcoord + poly_vert[1]*dimension;
    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, coordX, max_small_magnitude, 
       cos_min_angle, flag_zero);

    for (NTYPE i0 = 1; i0 < num_poly_vert; i0++) {
      const NTYPE i1 = (i0+1)%num_poly_vert;
      vcoord0 = vcoord + poly_vert[i0]*dimension;
      vcoord1 = vcoord + poly_vert[i1]*dimension;

      compute_cos_min_triangle_angle
        (dimension, vcoord0, vcoord1, coordX, max_small_magnitude, 
         cos_min_angle_i, flag_zero_i);
      if (flag_zero_i) { flag_zero = flag_zero_i; }
      if (cos_min_angle_i > cos_min_angle) {
        cos_min_angle = cos_min_angle_i;
        index_of_triangle_with_min_cos = i0;
      }
    }
  }
  */


  /*!
   *  Compute triangulation vertex which maximizes min triangulation angle.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] List of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] iv_max_min Index in poly_vert[] of vertex 
   *    which maximizes min triangulation angle.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if all triangulations have some triangle
   *    with two or three edges with length less than or equal 
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename VTYPE, typename MTYPE, 
            typename IVTYPE, typename COS_TYPE>
  void compute_vertex_to_max_min_triangulation_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    // Initialize
    cos_min_angle = 1;
    flag_zero = true;
    iv_max_min = 0;

    for (VTYPE iv = 0; iv < num_poly_vert; iv++) {
      compute_cos_min_triangulation_angle
        (dimension, vcoord, num_poly_vert, poly_vert,
         max_small_magnitude, iv, cosA, flagA);

      if (!flagA && cosA < cos_min_angle) {
        flag_zero = false;
        iv_max_min = iv;
        cos_min_angle = cosA;
      }
    }

    return;
  }


  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a pentagon into five triangles.
   *  - Quadrilateral is split into five triangles formed from coord4[]
   *    and each of the five pentagon edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the five triagles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vcoord0[] Coordinates of pentagon vertex 0.
   *  @param vcoord1[] Coordinates of pentagon vertex 1.
   *  @param vcoord2[] Coordinates of pentagon vertex 2.
   *  @param vcoord3[] Coordinates of pentagon vertex 3.
   *  @param vcoord4[] Coordinates of pentagon vertex 4.
   *  @param vcoordX[] Coordinates of vertex at center of star.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename CTYPEX,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_pentagon_tri5_angle
  (const DTYPE dimension, const CTYPE0 * vcoord0, const CTYPE1 * vcoord1, 
   const CTYPE2 * vcoord2, const CTYPE3 * vcoord3, const CTYPE4 * vcoord4,
   const CTYPEX * vcoordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_PENTAGON = 5;
    COS_TYPE cos[NUM_VERT_PER_PENTAGON];
    bool flag_tri_zero[NUM_VERT_PER_PENTAGON];

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoordX, max_small_magnitude,
       cos[0], flag_tri_zero[0]);
    compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord2, vcoordX, max_small_magnitude,
       cos[1], flag_tri_zero[1]);
    compute_cos_min_triangle_angle
      (dimension, vcoord2, vcoord3, vcoordX, max_small_magnitude,
       cos[2], flag_tri_zero[2]);
    compute_cos_min_triangle_angle
      (dimension, vcoord3, vcoord4, vcoordX, max_small_magnitude,
       cos[3], flag_tri_zero[3]);
    compute_cos_min_triangle_angle
      (dimension, vcoord4, vcoord0, vcoordX, max_small_magnitude,
       cos[4], flag_tri_zero[4]);

    bool is_cos_min_angle_set = false;
    for (int i = 0; i < NUM_VERT_PER_PENTAGON; i++) {

      if (flag_tri_zero[i]) {
        flag_zero = true;
        cos_min_angle = 0;
        return;
      }

      if (!is_cos_min_angle_set || cos[i] > cos_min_angle) {
        cos_min_angle = cos[i];
        is_cos_min_angle_set = true;
      }
    }

    flag_zero = false;
  }



  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a pentagon into five triangles.
   *  - Version with pentagon vertices listed in array pentagon_vert[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vertex_coord[] Vertex coordinates.
   *    vertex_coord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param pentagon_vert[] List of pentagon vertices in clockwise
   *    or counter-clockwise order around the pentagon.
   *  @param ivX Triangulate into five triangles by connecting each
   *     pentagon edge to ivX.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE0, typename VTYPE1,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_pentagon_tri5_angle
  (const DTYPE dimension, const CTYPE vertex_coord[], 
   const VTYPE0 pentagon_vert[], const VTYPE1 ivX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const CTYPE * vcoord0 = vertex_coord + pentagon_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord + pentagon_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord + pentagon_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord + pentagon_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord + pentagon_vert[4]*dimension;
    const CTYPE * vcoordX = vertex_coord + ivX*dimension;

    compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoordX,
       max_small_magnitude, cos_min_angle, flag_zero);
  }

  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a pentagon into five triangles.
   *  - Version with pentagon vertices listed in array pentagon_vert[].
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE0, typename VTYPE1,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_pentagon_tri5_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const VTYPE0 pentagon_vert[], const VTYPE1 ivX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_pentagon_tri5_angle
      (dimension, IJK::vector2pointer(vertex_coord), pentagon_vert, ivX,
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /// Compute min angle in triangulation of a pentagon
  ///   with all triangles incident on vertex_coord0.
  template <typename CTYPE, typename MTYPE, typename COS_TYPE>
  void compute_cos_min_pentagon_triangulation_angle
  (const CTYPE vertex_coord0,
    const CTYPE vertex_coord1,
    const CTYPE vertex_coord2,
    const CTYPE vertex_coord3,
    const CTYPE vertex_coord4,
    const MTYPE max_small_magnitude,
    COS_TYPE & cos_min_angle, bool & flag_zero)
   {
     const int DIM3(3);
     COS_TYPE cos012, cos023, cos034;
     bool flag_zero_012, flag_zero_023, flag_zero_034;

     // Initialize.
     cos_min_angle = 1;
     flag_zero = true;

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord1, vertex_coord2, max_small_magnitude,
        cos012, flag_zero_012);
     if (flag_zero_012) { return; }

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord2, vertex_coord3, max_small_magnitude,
        cos023, flag_zero_023);
     if (flag_zero_023) { return; }

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord3, vertex_coord4, max_small_magnitude,
        cos034, flag_zero_034);
     if (flag_zero_034) { return; }

     cos_min_angle = cos012;
     if (cos_min_angle < cos023) { cos_min_angle = cos023; }
     if (cos_min_angle < cos034) { cos_min_angle = cos034; }

     flag_zero = false;
   }

   ///@}

   // **************************************************
   //! @name COMPUTE DISTANCE
   // **************************************************

   ///@{

   /// Compute the distance from pcoord[] to the closest grid facet
   ///   incident on edge.
   /// @tparam GRID_TYPE Must have function SpacingPtrConst().
   /// @param grid Regular grid.
   template <typename GRID_TYPE, typename CTYPE, 
             typename VTYPE, typename DIR_TYPE, typename DIST_TYPE>
   void compute_unscaled_distance_to_closest_facet
   (const GRID_TYPE & grid, const CTYPE pcoord[],
    const VTYPE endpoint0, const DIR_TYPE edge_dir, DIST_TYPE & distance)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

     const DTYPE dimension = grid.Dimension();
     IJK::ARRAY<CTYPE> end0_coord(dimension);

     grid.ComputeCoord(endpoint0, end0_coord.Ptr());

     // Initialize
     distance = 0;

     bool is_distance_set = false;
     for (DTYPE d = 0; d < dimension; d++) {
       if (d != edge_dir) {
         DIST_TYPE distance2 = end0_coord[d]-pcoord[d]/grid.Spacing(d);
         if (distance2 < 0) { distance2 = -distance2; }
         if (!is_distance_set || distance2 < distance) {
           distance = distance2; 
           is_distance_set = true;
         }
       }
     }
   }

   /// Compute the distance from pcoord[] to the closest and farthest 
   ///   grid facets incident on edge.
   /// @tparam GRID_TYPE Must have function SpacingPtrConst().
   /// @param grid Regular grid.
   template <typename GRID_TYPE, typename CTYPE, 
             typename VTYPE, typename DIR_TYPE, typename DIST_TYPE>
   void compute_unscaled_distances_to_facets
   (const GRID_TYPE & grid, const CTYPE pcoord[],
    const VTYPE endpoint0, const DIR_TYPE edge_dir, 
    DIST_TYPE & closest_distance, DIST_TYPE & farthest_distance)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

     const DTYPE dimension = grid.Dimension();
     IJK::ARRAY<CTYPE> end0_coord(dimension);

     grid.ComputeCoord(endpoint0, end0_coord.Ptr());

     // Initialize
     closest_distance = 0;
     farthest_distance = 1;

     bool is_distance_set = false;
     for (DTYPE d = 0; d < dimension; d++) {
       if (d != edge_dir) {
         DIST_TYPE distance2 = end0_coord[d]-pcoord[d]/grid.Spacing(d);
         if (distance2 < 0) { distance2 = -distance2; }
         if (!is_distance_set) {
           closest_distance = distance2;
           farthest_distance = distance2;
         }
         else {
           if (distance2 < closest_distance) 
             { closest_distance = distance2; }
           if (distance2 > farthest_distance) 
             { farthest_distance = distance2; }
         }
         is_distance_set = true;
       }
     }
   }

   /// Compute the minimum distance from the four quad vertices 
   ///   to the closest grid facet incident on edge (endpoint0, endpoint1).
   /// @tparam GRID_TYPE Must have function SpacingPtrConst().
   /// @param grid Regular grid.
   template <typename GRID_TYPE, typename CTYPE, 
             typename VTYPE0, typename VTYPE1, typename DIR_TYPE,
             typename DIST_TYPE>
   void compute_unscaled_min_distance_to_closest_facet
   (const GRID_TYPE & grid, const CTYPE vertex_coord[],
    const VTYPE0 quad_vert[],
    const VTYPE1 endpoint0, const DIR_TYPE edge_dir, DIST_TYPE & min_distance)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
     typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

     const DTYPE dimension = grid.Dimension();

     const NTYPE NUM_VERT_PER_QUAD = 4;

     const CTYPE * pcoord0 = vertex_coord+dimension*quad_vert[0];
     compute_unscaled_distance_to_closest_facet
       (grid, pcoord0, endpoint0, edge_dir, min_distance);

     for (NTYPE i = 1; i < NUM_VERT_PER_QUAD; i++) {
       const CTYPE * pcoord = vertex_coord+dimension*quad_vert[i];

       DIST_TYPE distance;
       compute_unscaled_distance_to_closest_facet
         (grid, pcoord, endpoint0, edge_dir, distance);

       if (distance < min_distance) 
         { min_distance = distance; }
     }
   }


   /// Compute the minimum distance from the four quad vertices 
   ///   to the closest and farthest grid facets incident on an edge.
   /// @tparam GRID_TYPE Must have function SpacingPtrConst().
   /// @param grid Regular grid.
   /// @param[out] min_distance Minimum distance from any quad vert
   ///    to the facets incident on the edge.
   /// @param[out] min_max_distance Minimum of the maximum distance 
   ///    from a quad vert to the facets incident on the edge.
   template <typename GRID_TYPE, typename CTYPE, 
             typename VTYPE0, typename VTYPE1, typename DIR_TYPE,
             typename DIST_TYPE>
   void compute_unscaled_distances_from_quad_vert_to_facets
   (const GRID_TYPE & grid, const CTYPE vertex_coord[],
    const VTYPE0 quad_vert[],
    const VTYPE1 endpoint0, const DIR_TYPE edge_dir, 
    DIST_TYPE & min_distance,
    DIST_TYPE & min_max_distance)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
     typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

     const DTYPE dimension = grid.Dimension();

     const NTYPE NUM_VERT_PER_QUAD = 4;

     const CTYPE * pcoord0 = vertex_coord+dimension*quad_vert[0];
     compute_unscaled_distances_to_facets
       (grid, pcoord0, endpoint0, edge_dir, min_distance, min_max_distance);

     for (NTYPE i = 1; i < NUM_VERT_PER_QUAD; i++) {
       const CTYPE * pcoord = vertex_coord+dimension*quad_vert[i];

       DIST_TYPE closest_distance, farthest_distance;
       compute_unscaled_distances_to_facets
         (grid, pcoord, endpoint0, edge_dir, 
          closest_distance, farthest_distance);

       if (closest_distance < min_distance)
         { min_distance = closest_distance; }
       if (farthest_distance < min_max_distance)
         { min_max_distance = farthest_distance; }
     }
   }


   /// Return true if some quad vertex is close to one of the facets
   ///   incident on edge (endpoint0, endpoint1).
   /// @tparam GRID_TYPE Must have function SpacingPtrConst()
   ///   and EdgeDirecton().
   /// @param grid Regular grid.
   template <typename GRID_TYPE, typename CTYPE, 
             typename VTYPE0, typename VTYPE1, 
             typename DIR_TYPE,
             typename DIST_TYPE>
   bool is_some_quad_vert_close_to_facet
   (const GRID_TYPE & grid, const CTYPE vertex_coord[],
    const VTYPE0 quad_vert[], const VTYPE1 endpoint0, const DIR_TYPE edge_dir,
    const DIST_TYPE min_distance)
   {
     CTYPE distance;
     compute_unscaled_min_distance_to_closest_facet
       (grid, vertex_coord, quad_vert, endpoint0, edge_dir, distance);

     if (distance < min_distance) { return(true); }
     else { return(false); }
   }

   ///@}


   // **************************************************
   /// @name TRIANGULATE POLYGON
   // **************************************************

   ///@{

   /// Triangulate a single polygon by adding diagonals from the polygon
   ///   vertex with the maximum angle
   /// - Split maximum polygon angle.
   /// - Polygon vertices are listed in clockwise or counter-clockwise order
   ///   around the polygon.
   /// - Add new triangles to vector tri_vert.
   /// - If some polygon edge has length (near) zero, the triangulation
   ///   will be from an arbitrary vertex.
   template <typename DTYPE, typename CTYPE, typename NTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE>
   void triangulate_polygon_split_max_angle
   (const DTYPE dimension,
    const CTYPE * vertex_coord,
    const NTYPE num_poly_vert,
    const VTYPE0 * poly_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     CTYPE min_cos;
     bool flag_zero;

     if (num_poly_vert < 3) { return; };

     NTYPE k_min_cos = 0;        // Index of the vertex with min_cos
     NTYPE iv0 = poly_vert[num_poly_vert-1];
     NTYPE iv1 = poly_vert[0];
     NTYPE iv2 = poly_vert[1];

     compute_cos_triangle_angle_coord_list
       (dimension, vertex_coord, iv0, iv1, iv2, max_small_magnitude,
        min_cos, flag_zero);

     for (NTYPE i = 1; i < num_poly_vert; i++) {
       iv0 = poly_vert[i-1];
       iv1 = poly_vert[i];
       iv2 = poly_vert[((i+1)%num_poly_vert)];
       CTYPE cos_v0;

       compute_cos_triangle_angle_coord_list
         (dimension, vertex_coord, iv0, iv1, iv2, max_small_magnitude,
          cos_v0, flag_zero);

       if (cos_v0 < min_cos) {
         k_min_cos = i;
         min_cos = cos_v0;
       }
     }

     triangulate_polygon(num_poly_vert, poly_vert, k_min_cos, tri_vert);
   }


   /// Triangulate list of polygons.
   template <typename DTYPE, typename CTYPE,
             typename NTYPE0, typename NTYPE1, 
             typename VTYPE0, typename VTYPE1,
             typename MTYPE, typename ITYPE>
   void triangulate_polygon_list_split_max_angle
   (const DTYPE dimension, const CTYPE * vertex_coord,
    const NTYPE0 * num_poly_vert, const VTYPE0 * poly_vert,
    const ITYPE * first_poly_vert, const NTYPE1 num_poly,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     for (NTYPE1 ipoly = 0; ipoly < num_poly; ipoly++) {
       triangulate_polygon_split_max_angle
         (dimension, vertex_coord, num_poly_vert[ipoly], 
          poly_vert+first_poly_vert[ipoly], max_small_magnitude, tri_vert);
     }
   }

   /// Triangulate list of polygons.
   /// - C++ STL vector format for vert_coord[], num_poly_vert[],
   ///   poly_vert[], and first_poly_vert[].
   template <typename DTYPE, typename CTYPE, typename NTYPE,
             typename VTYPE0, typename VTYPE1,
             typename MTYPE, typename ITYPE>
   void triangulate_polygon_list_split_max_angle
   (const DTYPE dimension, const std::vector<CTYPE> & vert_coord,
    const std::vector<NTYPE> & num_poly_vert, 
    const std::vector<VTYPE0> & poly_vert,
    const std::vector<ITYPE> & first_poly_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     triangulate_polygon_list_split_max_angle
       (dimension, IJK::vector2pointer(vert_coord), 
        IJK::vector2pointer(num_poly_vert), IJK::vector2pointer(poly_vert),
        IJK::vector2pointer(first_poly_vert), num_poly_vert.size(),
        max_small_magnitude, tri_vert);
   }

   ///@}


   // **************************************************
   /// @name TRIANGULATE QUADRILATERAL
   // **************************************************

   ///@{

   /// Triangulate a list of quadrilaterals.
   template <typename DTYPE, typename CTYPE, typename NTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE>
   void triangulate_quad_split_max_angle
   (const DTYPE dimension,
    const CTYPE * vert_coord,
    const VTYPE0 * quad_vert,
    const NTYPE num_quad,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     const NTYPE NUM_VERT_PER_QUAD = 4;

     for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
       NTYPE k = iquad*NUM_VERT_PER_QUAD;
       triangulate_polygon_split_max_angle
         (dimension, vert_coord, NUM_VERT_PER_QUAD, quad_vert+k, 
          max_small_magnitude, tri_vert);
     }
   }

   /// Triangulate a list of quadrilaterals.
   /// - C++ STL vector format for vert_coord and quad_vert.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE>
   void triangulate_quad_split_max_angle
   (const DTYPE dimension,
    const std::vector<CTYPE> & vert_coord,
    const std::vector<VTYPE0> & quad_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;

     SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
     triangulate_quad_split_max_angle
       (dimension, IJK::vector2pointer(vert_coord),
        IJK::vector2pointer(quad_vert), num_quad, max_small_magnitude,
        tri_vert);

   }

   /// Triangulate a single quad.
   /// - Maximize minimum triangle angle.
   /// - Quad vertices are listed in clockwise or counter-clockwise order
   ///   around the quadrilateral.
   /// - Add new triangles to vector tri_vert.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE>
   void triangulate_quad_max_min_angle
   (const DTYPE dimension,
    const CTYPE * vert_coord,
    const VTYPE0 * quad_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const CTYPE * vcoord0 = vert_coord+quad_vert[0]*dimension;
     const CTYPE * vcoord1 = vert_coord+quad_vert[1]*dimension;
     const CTYPE * vcoord2 = vert_coord+quad_vert[2]*dimension;
     const CTYPE * vcoord3 = vert_coord+quad_vert[3]*dimension;
     CTYPE cos_min_angle;
     bool flag_diag02;
     bool flag_zero;

     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
        cos_min_angle, flag_diag02, flag_zero);

     triangulate_quad_using_diagonal(quad_vert, flag_diag02, tri_vert);
   }


   /// Triangulate a list of quadrilaterals.
   template <typename DTYPE, typename CTYPE, typename NTYPE,
             typename VTYPE0, typename VTYPE1,
             typename MTYPE>
   void triangulate_quad_max_min_angle
   (const DTYPE dimension,
    const CTYPE * vert_coord,
    const VTYPE0 * quad_vert,
    const NTYPE num_quad,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     const NTYPE NUM_VERT_PER_QUAD = 4;

     for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
       NTYPE k = iquad*NUM_VERT_PER_QUAD;
       triangulate_quad_max_min_angle
         (dimension, vert_coord, quad_vert+k, max_small_magnitude, tri_vert);
     }
   }

   /// Triangulate a list of quadrilaterals.
   /// - C++ STL vector format for vert_coord and quad_vert.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE>
   void triangulate_quad_max_min_angle
   (const DTYPE dimension,
    const std::vector<CTYPE> & vert_coord,
    const std::vector<VTYPE0> & quad_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

     triangulate_quad_max_min_angle
       (dimension, IJK::vector2pointer(vert_coord),
        IJK::vector2pointer(quad_vert), num_quad, max_small_magnitude,
        tri_vert);
   }


   /// Triangulate polygon.  
   /// - All triangulation triangles are incident on one vertex.
   /// - Maximize minimum triangle angle.
   /// - Add new triangles to vector tri_vert.
   /// @param[out] index_tri_vert Index of triangulation vertex.
   ///    All triangulation triangles are incident on poly_vert[index_tri_vert].
   template <typename DTYPE, typename NTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE,
             typename COS_TYPE, typename IVTYPE>
   void triangulate_polygon_max_min_angle
   (const DTYPE dimension, const CTYPE * vert_coord, 
    const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert,
    IVTYPE & index_tri_vert,
    COS_TYPE & cos_min_angle)
   {
     bool flagA;

     compute_vertex_to_max_min_triangulation_angle
       (dimension, vert_coord, num_poly_vert, poly_vert,
        max_small_magnitude, index_tri_vert, cos_min_angle, flagA);

     if (flagA) {
       index_tri_vert = 0;
       cos_min_angle = 0;
     }

     triangulate_polygon(num_poly_vert, poly_vert, index_tri_vert, tri_vert);
   }

   /// Triangulate polygon.
   /// - All triangulation triangles are incident on one vertex.
   /// - Maximize minimum triangle angle.
   /// - Version which does not return index_tri_vert or cos_min_angle.
   template <typename DTYPE, typename NTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE>
   void triangulate_polygon_max_min_angle
   (const DTYPE dimension, const CTYPE * vert_coord, 
    const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     NTYPE index_tri_vert;
     CTYPE cos_min_angle;

     triangulate_polygon_max_min_angle
       (dimension, vert_coord, num_poly_vert, poly_vert,
        max_small_magnitude, tri_vert, index_tri_vert, cos_min_angle);
   }

   /// Triangulate polygon.
   /// - All triangulation triangles are incident on one vertex.
   /// - Maximize minimum triangle angle.
   /// - Version which does not return index_tri_vert or cos_min_angle.
   /// - C++ STL vector format for vert_coord and poly_vert.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename MTYPE>
   void triangulate_polygon_max_min_angle
   (const DTYPE dimension, 
    const std::vector<CTYPE> & vert_coord, 
    const std::vector<VTYPE0> & poly_vert,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE1> & tri_vert)
   {
     triangulate_polygon_max_min_angle
       (dimension, IJK::vector2pointer(vert_coord),
        poly_vert.size(), IJK::vector2pointer(poly_vert),
        max_small_magnitude, tri_vert);
   }

   /// Return true if vertex iw separates iv0 and iv1 from iv2 and iv3
   ///   in direction d.
   template <typename CTYPE, typename VTYPE0, typename VTYPE1, 
             typename DIR_TYPE>
   bool does_vertex_separate
   (const CTYPE vertex_coord[],
    const VTYPE0 iw, const VTYPE1 iv0, const VTYPE1 iv1,
    const VTYPE1 iv2, const VTYPE1 iv3, const DIR_TYPE d)
   {
     const DIR_TYPE DIM3(3);
     const CTYPE wc = vertex_coord[DIM3*iw+d];
     const CTYPE vc0 = vertex_coord[DIM3*iv0+d];
     const CTYPE vc1 = vertex_coord[DIM3*iv1+d];
     const CTYPE vc2 = vertex_coord[DIM3*iv2+d];
     const CTYPE vc3 = vertex_coord[DIM3*iv3+d];

     if (vc0 <= wc && vc1 <= wc && vc2 >= wc && vc3 >= wc)
       { return(true); }

     if (vc0 >= wc && vc1 >= wc && vc2 <= wc && vc3 <= wc)
       { return(true); }

     return(false); 
   }

   /// Return true if quad vertices surround the dual edge.
   template <typename CTYPE, typename VTYPE0, typename VTYPE1,
             typename GRID_EDGE_TYPE>
   bool does_quad_surround_dual_edge
   (const CTYPE vertex_coord[], 
    const VTYPE0 quad_vert[],
    const GRID_EDGE_TYPE & dual_edge, 
    const VTYPE1 vertex_on_dual_edge)
   {
     typedef typename GRID_EDGE_TYPE::DIRECTION_TYPE DIR_TYPE;

     const DIR_TYPE DIM3(3);
     const DIR_TYPE edge_dir = dual_edge.direction;

     const DIR_TYPE d1 = (edge_dir+1)%DIM3;
     const DIR_TYPE d2 = (edge_dir+2)%DIM3;

     if (does_vertex_separate
         (vertex_coord, vertex_on_dual_edge, 
          quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3], d1)) {

       if (does_vertex_separate
           (vertex_coord, vertex_on_dual_edge, 
            quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[0], d2)) {

         return(true);
       }
     }

     if (does_vertex_separate
         (vertex_coord, vertex_on_dual_edge, 
          quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3], d2)) {

       if (does_vertex_separate
           (vertex_coord, vertex_on_dual_edge, 
            quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[0], d1)) {

         return(true);
       }
     }

     return(false);
   }


   /// Triangulate a set of quadrilaterals.
   /// - Add 4 triangles when all quad vertices are min distance from facets
   ///   incident on dual quad edge.
   ///   <br> Otherwise, use diagonal which minimizes quad angle.
   /// - Split quad i into four triangles by adding a triangle between
   ///   each quad edge and vertex (first_vertex+i).
   /// - Quadrilateral vertices are listed in clockwise or counter-clockwise 
   ///   order around the polygon.
   /// - Add new triangles to vector tri_vert.
   template <typename GRID_TYPE, typename CTYPE, typename NTYPE, 
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename GRID_EDGE_TYPE, typename MTYPE, typename DIST_TYPE>
   void triangulate_quad_tri4_by_distance
   (const GRID_TYPE & grid, const CTYPE vertex_coord[],
    const VTYPE0 quad_vert[], const NTYPE num_quad,
    const GRID_EDGE_TYPE dual_edge[], const VTYPE1 first_vertex, 
    const DIST_TYPE min_distance,  const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

     const DTYPE dimension = grid.Dimension();
     const NTYPE NUM_VERT_PER_QUAD = 4;

     for (NTYPE iquad = 0; iquad < num_quad; iquad++) {

       const NTYPE k = iquad*NUM_VERT_PER_QUAD;
       VTYPE0 iv0 = dual_edge[iquad].Endpoint0();
       DTYPE edge_direction = dual_edge[iquad].Direction();

       if (is_some_quad_vert_close_to_facet
           (grid, vertex_coord, quad_vert+k, iv0, 
            edge_direction, min_distance)) {

         // Split by some diagonal into two triangles.
         triangulate_quad_max_min_angle
           (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
            tri_vert);
       }
       else {
         // Split into four triangles.
         triangulate_polygon_with_vertex
           (NUM_VERT_PER_QUAD, quad_vert+k, first_vertex+iquad, tri_vert);
       }

     }
   }

   /// Triangulate a set of quadrilaterals.
   /// - Add 4 triangles when all quad vertices are min distance from facets
   ///   incident on dual quad edge.
   ///   <br> Otherwise, use diagonal which minimizes quad angle.
   /// - Split quad i into four triangles by adding a triangle between
   ///   each quad edge and vertex (first_vertex+i).
   /// - Quadrilateral vertices are listed in clockwise or counter-clockwise
   ///   order around the polygon.
   /// - Add new triangles to vector tri_vert.
   /// - C++ STL vector format for vert_coord and quad_vert and dual_edge[].
   template <typename GRID_TYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename GRID_EDGE_TYPE, typename MTYPE, typename DIST_TYPE>
   void triangulate_quad_tri4_by_distance
   (const GRID_TYPE & grid, const std::vector<CTYPE> & vertex_coord,
    const std::vector<VTYPE0> & quad_vert, 
    const std::vector<GRID_EDGE_TYPE> & dual_edge, 
    const VTYPE1 first_vertex, 
    const DIST_TYPE min_distance,  const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

     triangulate_quad_tri4_by_distance
       (grid, IJK::vector2pointer(vertex_coord),
        IJK::vector2pointer(quad_vert), num_quad,
        IJK::vector2pointer(dual_edge), first_vertex,
        min_distance, max_small_magnitude, tri_vert);
   }


   /// Triangulate a single quad into two or four triangles.
   /// - Minimize maximum triangle angle.
   /// - Quad vertices are listed in clockwise or counter-clockwise order
   ///   around the quadrilateral.
   /// - Add new triangles to vector tri_vert.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_quad_tri4_max_min_angle
   (const DTYPE dimension,
    const CTYPE vert_coord[],
    const VTYPE0 quad_vert[],
    const VTYPE1 ivert_on_dual_edge,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const CTYPE * vcoord0 = vert_coord+quad_vert[0]*dimension;
     const CTYPE * vcoord1 = vert_coord+quad_vert[1]*dimension;
     const CTYPE * vcoord2 = vert_coord+quad_vert[2]*dimension;
     const CTYPE * vcoord3 = vert_coord+quad_vert[3]*dimension;
     const CTYPE * vcoordX = vert_coord+ivert_on_dual_edge*dimension;
     CTYPE cos_min_angle;
     bool flag_tri4, flag_diag02, flag_zero;

     compute_cos_min_quad_tri02_tri13_tri4_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoordX,
        max_small_magnitude, cos_min_angle, flag_tri4, flag_diag02, flag_zero);


     if (flag_tri4) {
       triangulate_polygon_with_vertex
         (NUM_VERT_PER_QUAD, quad_vert, ivert_on_dual_edge, tri_vert);
     }
     else {
       triangulate_quad_using_diagonal(quad_vert, flag_diag02, tri_vert);
     }
   }


   /// Triangulate a list of isosurface quadrilaterals.
   /// - Each quadrilateral is dual to a grid edge.
   /// @min_dist_allow_tri4 Allow triangulation into 4 triangles when all
   ///   quad vertices are min_distance_allow_tri4 distance from facets
   ///   incident on dual quad edges.
   ///   If some quad vertex is closer than min_distance_allow_tri4 to
   ///   a facet incident on the dual quad edge, use a diagonal triangulation
   ///   into two triangles.
   template <typename GRID_TYPE, typename CTYPE, typename NTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename GRID_EDGE_TYPE,
             typename DIST_TYPE, typename MTYPE>
   void triangulate_dual_quad_tri4_max_min_angle
   (const GRID_TYPE & grid, 
    const CTYPE vertex_coord[],
    const VTYPE0 quad_vert[],
    const NTYPE num_quad,
    const GRID_EDGE_TYPE dual_edge[],
    const VTYPE1 first_vertex_on_dual_edge, 
    const DIST_TYPE min_distance_allow_tri4,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

     const DTYPE dimension = grid.Dimension();
     const DTYPE DIM3(3);
     const NTYPE NUM_VERT_PER_QUAD(4);
     CTYPE min_distance, min_max_distance;
     IJK::PROCEDURE_ERROR error("triangulate_dual_quad_tri4_max_min_angle");

     if (num_quad > 0 && dual_edge == NULL) {
       error.AddMessage
         ("Programming error.  Array dual_edge[] is empty.");
       throw error;
     }

     for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
       const VTYPE2 iw = first_vertex_on_dual_edge+iquad;
       const NTYPE k = iquad*NUM_VERT_PER_QUAD;

       if (dimension == DIM3 &&
           !does_quad_surround_dual_edge
           (vertex_coord, quad_vert+k, dual_edge[iquad], iw)) {

         // Split by some diagonal into two triangles.
         triangulate_quad_max_min_angle
           (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
            tri_vert);
       }
       else {
         VTYPE0 iv0 = dual_edge[iquad].Endpoint0();
         DTYPE edge_direction = dual_edge[iquad].Direction();

         compute_unscaled_distances_from_quad_vert_to_facets
           (grid, vertex_coord, quad_vert+k, iv0, edge_direction, 
            min_distance, min_max_distance);

         if (min_max_distance > min_distance_allow_tri4) {
           // Split into two or four triangles, whichever maximizes min angle.
           triangulate_quad_tri4_max_min_angle
             (dimension, vertex_coord, quad_vert+k, iw, max_small_magnitude, 
              tri_vert);
         }
         else {
           // Split by some diagonal into two triangles.
           triangulate_quad_max_min_angle
             (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
              tri_vert);
         }
       }
     }
   }

   /// Triangulate a list of isosurface quadrilaterals.
   /// - Each quadrilateral is dual to a grid edge.
   /// - C++ STL vector format for vert_coord and quad_vert.
   template <typename GRID_TYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2, 
             typename GRID_EDGE_TYPE,
             typename DIST_TYPE, typename MTYPE>
   void triangulate_dual_quad_tri4_max_min_angle
   (const GRID_TYPE & grid, 
    const std::vector<CTYPE> & vert_coord,
    const std::vector<VTYPE0> & quad_vert,
    const std::vector<GRID_EDGE_TYPE> & dual_edge,
    const VTYPE1 first_vertex_on_dual_edge, 
    const DIST_TYPE min_distance_allow_tri4,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

     if (num_quad != dual_edge.size()) {
       IJK::PROCEDURE_ERROR error("triangulate_dual_quad_tri4_max_min_angle");
       error.AddMessage
         ("Programming error. Number of dual edges does not equal number of isosurface quadrilaterals.");
       error.AddMessage("  Number of dual edges: ", dual_edge.size(), ".");
       error.AddMessage("  Number of isosurface quadrilaterals: ", 
                        num_quad, ".");
       throw error;
     }

     triangulate_dual_quad_tri4_max_min_angle
       (grid, IJK::vector2pointer(vert_coord),
        IJK::vector2pointer(quad_vert), num_quad, 
        IJK::vector2pointer(dual_edge),
        first_vertex_on_dual_edge,
        min_distance_allow_tri4, max_small_magnitude,
        tri_vert);
   }


   /// Triangulate a list of isosurface quadrilaterals.
   template <typename GRID_TYPE, typename CTYPE, typename NTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_quad_tri4_max_min_angle
   (const GRID_TYPE & grid, 
    const CTYPE vertex_coord[],
    const VTYPE0 quad_vert[],
    const NTYPE num_quad,
    const VTYPE1 first_vertex_on_dual_edge, 
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

     const DTYPE dimension = grid.Dimension();
     const DTYPE DIM3(3);
     const NTYPE NUM_VERT_PER_QUAD(4);
     CTYPE min_distance, min_max_distance;
     IJK::PROCEDURE_ERROR error("triangulate_quad_tri4_max_min_angle");

     if (num_quad > 0 && quad_vert == NULL) {
       error.AddMessage
         ("Programming error.  Array quad_vert[] is empty.");
       throw error;
     }

     for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
       const VTYPE2 iw = first_vertex_on_dual_edge+iquad;
       const NTYPE k = iquad*NUM_VERT_PER_QUAD;

       // Split into two or four triangles, whichever maximizes min angle.
       triangulate_quad_tri4_max_min_angle
         (dimension, vertex_coord, quad_vert+k, iw, max_small_magnitude, 
          tri_vert);
     }
   }

   /// Triangulate a list of isosurface quadrilaterals.
   /// - C++ STL vector format for vert_coord and quad_vert.
   template <typename GRID_TYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2, 
             typename MTYPE>
   void triangulate_quad_tri4_max_min_angle
   (const GRID_TYPE & grid, 
    const std::vector<CTYPE> & vert_coord,
    const std::vector<VTYPE0> & quad_vert,
    const VTYPE1 first_additional_vertex,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

     triangulate_quad_tri4_max_min_angle
       (grid, IJK::vector2pointer(vert_coord),
        IJK::vector2pointer(quad_vert), num_quad, 
        first_additional_vertex,
        max_small_magnitude,
        tri_vert);
   }

   /// Triangulate a single quad into four triangles or two triangles
   ///   using diag 02.
   /// - Minimize maximum triangle angle.
   /// - Quad vertices are listed in clockwise or counter-clockwise order
   ///   around the quadrilateral.
   /// - Add new triangles to vector tri_vert.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_quad_tri4_or_tri02_max_min_angle
   (const DTYPE dimension,
    const CTYPE * vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivert_on_dual_edge,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const CTYPE * vcoord0 = vert_coord+quad_vert[0]*dimension;
     const CTYPE * vcoord1 = vert_coord+quad_vert[1]*dimension;
     const CTYPE * vcoord2 = vert_coord+quad_vert[2]*dimension;
     const CTYPE * vcoord3 = vert_coord+quad_vert[3]*dimension;
     const CTYPE * vcoord4 = vert_coord+ivert_on_dual_edge*dimension;
     CTYPE cos_min_diag02, cos_min_tri4;
     bool flag_zero_diag02, flag_zero_tri4;

     compute_cos_min_quad_tri02_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
        cos_min_diag02, flag_zero_diag02);
     compute_cos_min_quad_tri4_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4,
        max_small_magnitude, cos_min_tri4, flag_zero_tri4);

     bool flag_tri4;

     if (flag_zero_tri4) { 
       flag_tri4 = false;
     }
     else if (flag_zero_diag02) {
       flag_tri4 = true;
     }
     else {
       if (cos_min_diag02 < cos_min_tri4) 
         { flag_tri4 = false; }
       else
         { flag_tri4 = true; }
     }

     if (flag_tri4) {
       triangulate_polygon_with_vertex
         (NUM_VERT_PER_QUAD, quad_vert, ivert_on_dual_edge, tri_vert);
     }
     else {
       triangulate_polygon(NUM_VERT_PER_QUAD, quad_vert, 0, tri_vert);
     }
   }


   /// Triangulate a single quad into four triangles or two triangles
   ///   using diag 02.
   /// - Minimize maximum triangle angle.
   /// - C++ STL vector format for vert_coord.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_quad_tri4_or_tri02_max_min_angle
   (const DTYPE dimension,
    const std::vector<CTYPE> & vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivert_on_dual_edge,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     triangulate_quad_tri4_or_tri02_max_min_angle
       (dimension, IJK::vector2pointer(vert_coord), quad_vert,
        ivert_on_dual_edge, max_small_magnitude, tri_vert);
   }


   /// Triangulate a single quad into four triangles or two triangles
   ///   using diag 13.
   /// - Minimize maximum triangle angle.
   /// - Quad vertices are listed in clockwise or counter-clockwise order
   ///   around the quadrilateral.
   /// - Add new triangles to vector tri_vert.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_quad_tri4_or_tri13_max_min_angle
   (const DTYPE dimension,
    const CTYPE * vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivert_on_dual_edge,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

     const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
     const CTYPE * vcoord0 = vert_coord+quad_vert[0]*dimension;
     const CTYPE * vcoord1 = vert_coord+quad_vert[1]*dimension;
     const CTYPE * vcoord2 = vert_coord+quad_vert[2]*dimension;
     const CTYPE * vcoord3 = vert_coord+quad_vert[3]*dimension;
     const CTYPE * vcoord4 = vert_coord+ivert_on_dual_edge*dimension;
     CTYPE cos_min_diag13, cos_min_tri4;
     bool flag_zero_diag13, flag_zero_tri4;

     compute_cos_min_quad_tri13_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
        cos_min_diag13, flag_zero_diag13);
     compute_cos_min_quad_tri4_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4,
        max_small_magnitude, cos_min_tri4, flag_zero_tri4);

     bool flag_tri4;


     if (flag_zero_tri4) { 
       flag_tri4 = false;
     }
     else if (flag_zero_diag13) {
       flag_tri4 = true;
     }
     else {
       if (cos_min_diag13 < cos_min_tri4) 
         { flag_tri4 = false; }
       else
         { flag_tri4 = true; }
     }

     if (flag_tri4) {
       triangulate_polygon_with_vertex
         (NUM_VERT_PER_QUAD, quad_vert, ivert_on_dual_edge, tri_vert);
     }
     else {
       triangulate_polygon(NUM_VERT_PER_QUAD, quad_vert, 1, tri_vert);
     }
   }


   /// Triangulate a single quad into four triangles or two triangles
   ///   using diag 13.
   /// - Minimize maximum triangle angle.
   /// - C++ STL vector format for vert_coord.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_quad_tri4_or_tri13_max_min_angle
   (const DTYPE dimension,
    const std::vector<CTYPE> & vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivert_on_dual_edge,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     triangulate_quad_tri4_or_tri13_max_min_angle
       (dimension, IJK::vector2pointer(vert_coord), quad_vert,
        ivert_on_dual_edge, max_small_magnitude, tri_vert);
   }

   ///@}


   // **************************************************
   /// @name TRIANGULATE PENTAGON
   // **************************************************

   ///@{

   /// Compute vertex to maximize min triangulation angle of a pentagon.
   /// - Faster than calling compute_vertex_to_max_min_triangulation_angle.
   /// @param vertex_coord0 Coordinates of vertex 0.
   /// @param vertex_coord1 Coordinates of vertex 1.
   /// @param vertex_coord2 Coordinates of vertex 2.
   /// @param vertex_coord3 Coordinates of vertex 3.
   /// @param vertex_coord4 Coordinates of vertex 4.
   /// @param max_small_magnitude Vectors with magnitude less than or
   ///    equal to max_small_magnitude are set to 0.
   /// @pre max_small_magnitude >= 0.
   /// @param[out] iv_max_min Index in pentagon_vert[] of vertex 
   ///    which maximizes min triangulation angle.
   ///  @param[out] cos_min_angle Cosine of min triange angle.
   ///  @param[out] flag_zero True, if all triangulations have some triangle
   ///    with two or three edges with length less than or equal 
   ///    to max_small_magnitude.
   template <typename CTYPE, typename MTYPE, typename IVTYPE, typename COS_TYPE>
   void compute_vertex_to_max_min_pentagon_triangulation_angle
   (const CTYPE vertex_coord0[], const CTYPE vertex_coord1[], 
    const CTYPE vertex_coord2[], const CTYPE vertex_coord3[], 
    const CTYPE vertex_coord4[],
    const MTYPE max_small_magnitude,
    IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
   {
     const int DIM3(3);
     const int NUM_VERT_PER_PENTAGON(5);
     const CTYPE * vcoord[NUM_VERT_PER_PENTAGON] =
       { vertex_coord0, vertex_coord1, vertex_coord2, 
         vertex_coord3, vertex_coord4 };
     CTYPE cos_ear_triangle[NUM_VERT_PER_PENTAGON];
     bool flag_zero_ear[NUM_VERT_PER_PENTAGON];
     CTYPE cos_center_triangle[NUM_VERT_PER_PENTAGON];
     bool flag_zero_center_triangle[NUM_VERT_PER_PENTAGON];

     // Initialize
     iv_max_min = 0;
     cos_min_angle = 1;
     flag_zero = true;

     for (int i0 = 0; i0 < NUM_VERT_PER_PENTAGON; i0++) {
       const int i1 = (i0+1)%NUM_VERT_PER_PENTAGON;
       const int i2 = (i0+2)%NUM_VERT_PER_PENTAGON;
       const int i4 = (i0+4)%NUM_VERT_PER_PENTAGON;

       compute_cos_min_triangle_angle
         (DIM3, vcoord[i0], vcoord[i1], vcoord[i2], max_small_magnitude, 
          cos_ear_triangle[i1], flag_zero_ear[i1]);

       compute_cos_min_triangle_angle
         (DIM3, vcoord[i0], vcoord[i2], vcoord[i4], max_small_magnitude, 
          cos_center_triangle[i2], flag_zero_center_triangle[i2]);
     }

     bool flag_found = false;
     for (int i0 = 0; i0 < NUM_VERT_PER_PENTAGON; i0++) {
       const int i1 = (i0+1)%NUM_VERT_PER_PENTAGON;
       const int i2 = (i0+2)%NUM_VERT_PER_PENTAGON;

       if (flag_zero_ear[i0] || flag_zero_ear[i2] ||
           flag_zero_center_triangle[i1]) { continue; }

       COS_TYPE cos_min_i1 = cos_center_triangle[i1];
       if (cos_ear_triangle[i0] > cos_min_i1)
         { cos_min_i1 = cos_ear_triangle[i0]; }
       if (cos_ear_triangle[i2] > cos_min_i1)
         { cos_min_i1 = cos_ear_triangle[i2]; }

       if (flag_found) {
         if (cos_min_angle > cos_min_i1) {
           iv_max_min = i1;
           cos_min_angle = cos_min_i1;
         }
       }
       else {
         iv_max_min = i1;
         cos_min_angle = cos_min_i1;
         flag_found = true;
       }
     }

     flag_zero = !flag_found;
   }

   /// Compute vertex to maximize min triangulation angle of a pentagon.
   /// - Faster than calling compute_vertex_to_max_min_triangulation_angle.
   /// @param vertex_coord[] List of vertex coordinates.
   ///   vertex_coord[i*3+j] is j'th coordinate of i'th vertex.
   /// @param pentagon_vert[] List of pentagon vertices in clockwise
   ///   or counter-clockwise order around the pentagon.
   /// @param max_small_magnitude Vectors with magnitude less than or
   ///    equal to max_small_magnitude are set to 0.
   /// @pre max_small_magnitude >= 0.
   /// @param[out] iv_max_min Index in pentagon_vert[] of vertex 
   ///    which maximizes min triangulation angle.
   ///  @param[out] cos_min_angle Cosine of min triange angle.
   ///  @param[out] flag_zero True, if all triangulations have some triangle
   ///    with two or three edges with length less than or equal 
   ///    to max_small_magnitude.
   template <typename CTYPE, typename VTYPE, typename MTYPE,
             typename IVTYPE, typename COS_TYPE>
   void compute_vertex_to_max_min_pentagon_triangulation_angle
   (const CTYPE vertex_coord[],
    const VTYPE pentagon_vert[],
    const MTYPE max_small_magnitude,
    IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
   {
     const int DIM3(3);
     const int NUM_VERT_PER_PENTAGON(5);
     const CTYPE * vcoord[NUM_VERT_PER_PENTAGON] =
       { vertex_coord+pentagon_vert[0]*DIM3,
         vertex_coord+pentagon_vert[1]*DIM3,
         vertex_coord+pentagon_vert[2]*DIM3,
         vertex_coord+pentagon_vert[3]*DIM3,
         vertex_coord+pentagon_vert[4]*DIM3 };

     compute_vertex_to_max_min_pentagon_triangulation_angle
       (vcoord[0], vcoord[1], vcoord[2], vcoord[3], vcoord[4],
        max_small_magnitude, iv_max_min, cos_min_angle, flag_zero);
   }

   /// Compute vertex to maximize min triangulation angle of a pentagon.
   /// Version using C++ STL vector for array vertex_coord[].
   /// - Faster than calling compute_vertex_to_max_min_triangulation_angle.
   /// @param vertex_coord[] List of vertex coordinates.
   ///   vertex_coord[i*3+j] is j'th coordinate of i'th vertex.
   /// @param pentagon_vert[] List of pentagon vertices in clockwise
   ///   or counter-clockwise order around the pentagon.
   /// @param max_small_magnitude Vectors with magnitude less than or
   ///    equal to max_small_magnitude are set to 0.
   /// @pre max_small_magnitude >= 0.
   /// @param[out] iv_max_min Index in pentagon_vert[] of vertex 
   ///    which maximizes min triangulation angle.
   ///  @param[out] cos_min_angle Cosine of min triange angle.
   ///  @param[out] flag_zero True, if all triangulations have some triangle
   ///    with two or three edges with length less than or equal 
   ///    to max_small_magnitude.
   template <typename CTYPE, typename VTYPE, typename MTYPE,
             typename IVTYPE, typename COS_TYPE>
   void compute_vertex_to_max_min_pentagon_triangulation_angle
   (const std::vector<CTYPE> & vertex_coord,
    const VTYPE pentagon_vert[],
    const MTYPE max_small_magnitude,
    IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
    {
      compute_vertex_to_max_min_pentagon_triangulation_angle
        (IJK::vector2pointer(vertex_coord), pentagon_vert,
         max_small_magnitude, iv_max_min, cos_min_angle, flag_zero);
    }


   /// Triangulate a pentagon.
   /// - Maximize minimum triangle angle.
   /// - Pentagon vertices are listed in clockwise or counter-clockwise order
   ///   around the quadrilateral.
   /// - Add new triangles to vector tri_vert.
   template <typename CTYPE, typename VTYPE, typename MTYPE,
             typename TRI_VTYPE>
   void triangulate_pentagon_max_min_angle
   (const CTYPE vertex_coord[],
    const VTYPE pentagon_vert[],
    const MTYPE max_small_magnitude,
    std::vector<TRI_VTYPE> & tri_vert)
   {
     const int NUM_VERT_PER_PENTAGON(5);
     int iv;
     CTYPE cos_min_angle;
     bool flag_zero;

     compute_vertex_to_max_min_pentagon_triangulation_angle
       (vertex_coord, pentagon_vert, max_small_magnitude, iv,
        cos_min_angle, flag_zero);

     if (flag_zero) {
       triangulate_pentagon(pentagon_vert, 0, tri_vert);
     }
     else {
       triangulate_pentagon(pentagon_vert, iv, tri_vert);
     }
   }

   /// Triangulate a pentagon into three or five triangles.
   /// - Maximize minimum triangle angle.
   /// - Pentagon vertices are listed in clockwise or counter-clockwise order
   ///   around the quadrilateral.
   /// - Add new triangles to vector tri_vert.
   /// @param ivX Triangulate into five triangles by connecting each
   ///     pentagon edge to ivX.
   template <typename CTYPE, typename VTYPE0, typename VTYPE1,
             typename MTYPE, typename TRI_VTYPE>
   void triangulate_pentagon_tri5_max_min_angle
   (const std::vector<CTYPE> & vertex_coord,
    const VTYPE0 pentagon_vert[],
    const VTYPE1 ivX,
    const MTYPE max_small_magnitude,
    std::vector<TRI_VTYPE> & tri_vert)
   {
     const int DIM3(3);
     const int NUM_VERT_PER_PENTAGON(5);
     int iv;
     CTYPE cos_min_angle, cos_min_tri5_angle;
     bool flag_zero, flag_tri5_zero;

     compute_vertex_to_max_min_pentagon_triangulation_angle
       (vertex_coord, pentagon_vert, max_small_magnitude, iv,
        cos_min_angle, flag_zero);

     compute_cos_min_pentagon_tri5_angle
       (DIM3, vertex_coord, pentagon_vert, ivX, max_small_magnitude,
        cos_min_tri5_angle, flag_tri5_zero);

     if (flag_zero && flag_tri5_zero) {
       triangulate_pentagon(pentagon_vert, 0, tri_vert);
     }
     else if (flag_tri5_zero) {
       triangulate_pentagon(pentagon_vert, iv, tri_vert);
     }
     else if (flag_zero) {
       triangulate_polygon_with_vertex
         (NUM_VERT_PER_PENTAGON, pentagon_vert, ivX, tri_vert);
     }
     else {
       // Triangulations into three and five vertices are both valid.
       if (cos_min_angle < cos_min_tri5_angle) {
         triangulate_pentagon(pentagon_vert, iv, tri_vert);
       }
       else {
         triangulate_polygon_with_vertex
           (NUM_VERT_PER_PENTAGON, pentagon_vert, ivX, tri_vert);
       }
     }
   }


   template <typename CTYPE, typename VTYPE, typename MTYPE,
             typename TRI_VTYPE>
   void triangulate_pentagon_max_min_angle
   (const std::vector<CTYPE> & vertex_coord,
    const VTYPE pentagon_vert[],
    const MTYPE max_small_magnitude,
    std::vector<TRI_VTYPE> & tri_vert)
   {
     triangulate_pentagon_max_min_angle
       (IJK::vector2pointer(vertex_coord), pentagon_vert, 
        max_small_magnitude, tri_vert);
   }


   template <typename CTYPE, typename VTYPE, typename NTYPE,
             typename MTYPE, typename TRI_VTYPE>
   void triangulate_pentagon_list_max_min_angle
   (const std::vector<CTYPE> & vertex_coord,
    const VTYPE pentagon_vert[],
    const NTYPE num_pentagon,
    const MTYPE max_small_magnitude,
    std::vector<TRI_VTYPE> & tri_vert)
   {
     const int NUM_VERT_PER_PENTAGON(5);

     for (int i = 0; i < num_pentagon; i++) {
       triangulate_pentagon_max_min_angle
         (vertex_coord, pentagon_vert+i*NUM_VERT_PER_PENTAGON, 
          max_small_magnitude, tri_vert);
     }
   }


   template <typename CTYPE, typename VTYPE,
             typename MTYPE, typename TRI_VTYPE>
   void triangulate_pentagon_list_max_min_angle
   (const std::vector<CTYPE> & vertex_coord,
    const std::vector<VTYPE> & pentagon_vert,
    const MTYPE max_small_magnitude,
    std::vector<TRI_VTYPE> & tri_vert)
   {
     const int NUM_VERT_PER_PENTAGON(5);
     const int num_pentagon = pentagon_vert.size()/NUM_VERT_PER_PENTAGON;

     triangulate_pentagon_list_max_min_angle
       (vertex_coord, IJK::vector2pointer(pentagon_vert), num_pentagon,
        max_small_magnitude, tri_vert);
   }


   /// Triangulate list of pentagons into three or five triangles each.
   /// @param first_new_vert First new vertex.
   ///   - Pentagon i is triangulated into five triangles by connecting
   ///     first_new_vert+i to each edge of pentagon i.
   template <typename CTYPE, typename VTYPE0, typename VTYPE1, 
             typename NTYPE, typename MTYPE, typename TRI_VTYPE>
   void triangulate_pentagon_list_tri5_max_min_angle
   (const std::vector<CTYPE> & vertex_coord,
    const VTYPE0 pentagon_vert[],
    const NTYPE num_pentagon,
    const VTYPE1 first_new_vert,
    const MTYPE max_small_magnitude,
    std::vector<TRI_VTYPE> & tri_vert)
   {
     const int NUM_VERT_PER_PENTAGON(5);

     for (int i = 0; i < num_pentagon; i++) {
       triangulate_pentagon_tri5_max_min_angle
         (vertex_coord, pentagon_vert+i*NUM_VERT_PER_PENTAGON,
          first_new_vert+i, max_small_magnitude, tri_vert);
     }
   }


   /// Triangulate list of pentagons into three or five triangles each.
   /// - Version using C++ STL vector to store array pentagon_vert[].
   template <typename CTYPE, typename VTYPE0, typename VTYPE1,
             typename MTYPE, typename TRI_VTYPE>
   void triangulate_pentagon_list_tri5_max_min_angle
   (const std::vector<CTYPE> & vertex_coord,
    const std::vector<VTYPE0> & pentagon_vert,
    const VTYPE1 first_new_vert,
    const MTYPE max_small_magnitude,
    std::vector<TRI_VTYPE> & tri_vert)
   {
     const int NUM_VERT_PER_PENTAGON(5);
     const int num_pentagon = pentagon_vert.size()/NUM_VERT_PER_PENTAGON;

     triangulate_pentagon_list_tri5_max_min_angle
       (vertex_coord, IJK::vector2pointer(pentagon_vert), num_pentagon,
        first_new_vert, max_small_magnitude, tri_vert);
   }

   ///@}


   // ***************************************************************
   /// @name TRIANGULATE HEXAGON
   // ***************************************************************

   /*!
    *  Compute the cosine of the smallest triangle angle
    *    in the triangulation of a hexahedron.
    *  - Hexadhedron is split into six triangles formed from coordX[]
    *    and each of the six hexahedron edges.
    *  - Note: Cosine of the smallest angle is the largest cosine 
    *    of the four triangles.
    *  @param dimension Coordinate dimension (= number of coordinates.)
    *  @param coord0[] Input coordinates.
    *  @param coord1[] Input coordinates.
    *  @param coord2[] Input coordinates.
    *  @param coord3[] Input coordinates.
    *  @param coord4[] Input coordinates.
    *  @param coord5[] Input coordinates.
    *  @param coordX[] Input coordinates.
    *  @param max_small_magnitude Vectors with magnitude less than or
    *           equal to max_small_magnitude are set to 0.
    *  @pre max_small_magnitude >= 0.
    *  @param[out] cos_min_angle Cosine of min triange angle.
    *  @param[out] flag_zero True, if some triangle has two or three edges
    *    with length less than or equal to max_small_magnitude.
    */
   template <typename DTYPE, 
             typename CTYPE0, typename CTYPE1, typename CTYPE2, 
             typename CTYPE3, typename CTYPE4, typename CTYPE5,
             typename CTYPEX, typename COS_TYPE,
             typename MTYPE>
   void compute_cos_min_hex_tri6_angle
   (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
    const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPE4 * coord4,
    const CTYPE5 * coord5, const CTYPEX * coordX,
    const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
   {
     const int NUM_VERT_PER_HEX = 6;
     COS_TYPE cos[NUM_VERT_PER_HEX];
     bool flag_tri_zero[NUM_VERT_PER_HEX];

     compute_cos_min_triangle_angle
       (dimension, coord0, coord1, coordX, max_small_magnitude, 
        cos[0], flag_tri_zero[0]);
     compute_cos_min_triangle_angle
       (dimension, coord1, coord2, coordX, max_small_magnitude, 
        cos[1], flag_tri_zero[1]);
     compute_cos_min_triangle_angle
       (dimension, coord2, coord3, coordX, max_small_magnitude, 
        cos[2], flag_tri_zero[2]);
     compute_cos_min_triangle_angle
       (dimension, coord3, coord4, coordX, max_small_magnitude, 
        cos[3], flag_tri_zero[3]);
     compute_cos_min_triangle_angle
       (dimension, coord4, coord5, coordX, max_small_magnitude, 
        cos[4], flag_tri_zero[4]);
       (dimension, coord4, coord5, coordX, max_small_magnitude, 
        cos[5], flag_tri_zero[5]);

     bool is_cos_min_angle_set = false;
     for (int i = 0; i < NUM_VERT_PER_HEX; i++) {

       if (flag_tri_zero[i]) {
         flag_zero = true;
         cos_min_angle = 0;
         return;
       }

       if (!is_cos_min_angle_set || cos[i] > cos_min_angle) {
         cos_min_angle = cos[i];
         is_cos_min_angle_set = true;
       }
     }

     flag_zero = false;
   }


   ///@{

   /// Compute triangulation which maximizes min triangulation angle 
   ///    of a hexagon and includes triangle 012.
   /// @param[out] ear[] Array of ears defining triangulation.
   /// @pre Array ear[] is pre-allocated to size at least 3.
   template <typename CTYPE, typename MTYPE,
             typename EAR_TYPE, typename COS_TYPE>
   void compute_hex_triangulation_with_max_min_angle_tri012
   (const CTYPE vertex_coord0,
    const CTYPE vertex_coord1,
    const CTYPE vertex_coord2,
    const CTYPE vertex_coord3,
    const CTYPE vertex_coord4,
    const CTYPE vertex_coord5,
    const MTYPE max_small_magnitude,
    EAR_TYPE ear[], COS_TYPE & cos_min_angle, bool & flag_zero)
   {
     const int DIM3(3);
     const int NUM_VERT_PER_HEXAGON(6);
     COS_TYPE cos_tri_012, cos_pentagon;
     bool flag_zero_tri_012, flag_zero_pentagon;
     int iv_max_min;

     // Initialize
     cos_min_angle = 1;
     flag_zero = true;
     ear[0] = 1;
     ear[1] = 2;
     ear[2] = 3;

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord1, vertex_coord2, max_small_magnitude,
        cos_tri_012, flag_zero_tri_012);

     if (!flag_zero_tri_012) {
       compute_vertex_to_max_min_pentagon_triangulation_angle
         (vertex_coord0, vertex_coord2, vertex_coord3, vertex_coord4, 
          vertex_coord5, max_small_magnitude, 
          iv_max_min, cos_pentagon, flag_zero_pentagon);

       if (!flag_zero_pentagon) {
         cos_min_angle = std::max(cos_tri_012, cos_pentagon);
         flag_zero = false;

         if (iv_max_min != 4) {
           ear[1] = (iv_max_min+2)%NUM_VERT_PER_HEXAGON;
           ear[2] = (iv_max_min+3)%NUM_VERT_PER_HEXAGON;
         }
         else {
           ear[1] = 0;
           ear[2] = 2;
         }
       }
     }
   }

   /// Compute triangulation which maximizes min triangulation angle 
   ///    of a hexagon and includes triangles 012 and 345.
   /// @param[out] ear[] Array of ears defining triangulation.
   /// @pre Array ear[] is pre-allocated to size at least 3.
   template <typename CTYPE, typename MTYPE, typename COS_TYPE>
   void compute_hex_triangulation_with_max_min_angle_tri012_tri345
   (const CTYPE vertex_coord0,
    const CTYPE vertex_coord1,
    const CTYPE vertex_coord2,
    const CTYPE vertex_coord3,
    const CTYPE vertex_coord4,
    const CTYPE vertex_coord5,
    const MTYPE max_small_magnitude,
    bool & flag_diag03,
    COS_TYPE & cos_min_angle, 
    bool & flag_zero)
   {
     const int DIM3(3);
     COS_TYPE cos012, cos345, cos_diag03, cos_diag25, cos_quad;
     bool flag_zero_012, flag_zero_345, flag_zero_diag03, flag_zero_diag25;

     // Initialize.
     cos_min_angle = 1;
     flag_zero = true;
     flag_diag03 = true;

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord1, vertex_coord2, max_small_magnitude,
        cos012, flag_zero_012);
     if (flag_zero_012) { return; }

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord3, vertex_coord4, vertex_coord5, max_small_magnitude,
        cos345, flag_zero_345);
     if (flag_zero_345) { return; }

     IJK::compute_cos_min_quad_tri02_angle
       (DIM3, vertex_coord0, vertex_coord2, vertex_coord3, vertex_coord5,
        max_small_magnitude, cos_diag03, flag_zero_diag03);
     IJK::compute_cos_min_quad_tri13_angle
       (DIM3, vertex_coord0, vertex_coord2, vertex_coord3, vertex_coord5,
        max_small_magnitude, cos_diag25, flag_zero_diag25);
     if (flag_zero_diag03 && flag_zero_diag25) { return; }

     if (flag_zero_diag03) { flag_diag03 = false;  }
     else if (flag_zero_diag25) { flag_diag03 = true; }
     else if (cos_diag25 < cos_diag03) { flag_diag03 = false; }
     else { flag_diag03 = true; }

     if (flag_diag03) { cos_quad = cos_diag03; }
     else { cos_quad = cos_diag25; }

     flag_zero = false;
   }

   /// Compute min angle in triangulation of hexagon
   ///   with all triangles incident on vertex_coord0.
   template <typename CTYPE, typename MTYPE, typename COS_TYPE>
   void compute_cos_min_hex_triangulation_angle
   (const CTYPE vertex_coord0,
    const CTYPE vertex_coord1,
    const CTYPE vertex_coord2,
    const CTYPE vertex_coord3,
    const CTYPE vertex_coord4,
    const CTYPE vertex_coord5,
    const MTYPE max_small_magnitude,
    COS_TYPE & cos_min_angle, bool & flag_zero)
   {
     const int DIM3(3);
     COS_TYPE cos012, cos023, cos034, cos045;
     bool flag_zero_012, flag_zero_023, flag_zero_034, flag_zero_045;

     // Initialize.
     cos_min_angle = 1;
     flag_zero = true;

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord1, vertex_coord2, max_small_magnitude,
        cos012, flag_zero_012);
     if (flag_zero_012) { return; }

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord2, vertex_coord3, max_small_magnitude,
        cos023, flag_zero_023);
     if (flag_zero_023) { return; }

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord3, vertex_coord4, max_small_magnitude,
        cos034, flag_zero_034);
     if (flag_zero_034) { return; }

     IJK::compute_cos_min_triangle_angle
       (DIM3, vertex_coord0, vertex_coord4, vertex_coord5, max_small_magnitude,
        cos045, flag_zero_045);
     if (flag_zero_045) { return; }

     cos_min_angle = cos012;
     if (cos_min_angle < cos023) { cos_min_angle = cos023; }
     if (cos_min_angle < cos034) { cos_min_angle = cos034; }
     if (cos_min_angle < cos045) { cos_min_angle = cos045; }

     flag_zero = false;
   }


   /// Compute triangulation which maximizes min triangulation angle 
   ///    of a hexagon.
   /// @param vertex_coord[] List of vertex coordinates.
   ///   vertex_coord[i*3+j] is j'th coordinate of i'th vertex.
   /// @param hex_vert[] List of hexagon vertices in clockwise
   ///   or counter-clockwise order around the hexagon.
   /// @param max_small_magnitude Vectors with magnitude less than or
   ///    equal to max_small_magnitude are set to 0.
   /// @pre max_small_magnitude >= 0.
   /// @param[out] ear[] Array of ears defining triangulation.
   /// @pre Array ear[] is pre-allocated to size at least 3.
   /// @param[out] cos_min_angle Cosine of min triange angle.
   /// @param[out] flag_zero True, if all triangulations have some triangle
   ///   with two or three edges with length less than or equal 
   ///   to max_small_magnitude.
   template <typename CTYPE, typename VTYPE, typename MTYPE,
             typename EAR_TYPE, typename COS_TYPE>
   void compute_hex_triangulation_with_max_min_angle
   (const CTYPE vertex_coord[],
    const VTYPE hex_vert[],
    const MTYPE max_small_magnitude,
    EAR_TYPE ear[], COS_TYPE & cos_min_angle, bool & flag_zero)
   {
     const int DIM3(3);
     const int NUM_VERT_PER_HEXAGON(6);
     const int NUM_EAR(3);
     const CTYPE * vcoord[NUM_VERT_PER_HEXAGON] =
       { vertex_coord+hex_vert[0]*DIM3,
         vertex_coord+hex_vert[1]*DIM3,
         vertex_coord+hex_vert[2]*DIM3,
         vertex_coord+hex_vert[3]*DIM3,
         vertex_coord+hex_vert[4]*DIM3,
         vertex_coord+hex_vert[5]*DIM3 };
     COS_TYPE cos_hex_tri_012, cos_hex_tri_501, cos_hex_tri_123_450;
     COS_TYPE cos_hex_tri_v3, cos_hex_tri_v4;
     bool flag_zero_012, flag_zero_501, flag_zero_123_450;
     bool flag_zero_v3, flag_zero_v4;
     EAR_TYPE ear_tri012[NUM_EAR], ear_tri501[NUM_EAR];
     bool flag_diag14;

     typedef enum { HEX_TRI_012, HEX_TRI_501, HEX_TRI_V3, HEX_TRI_V4, 
                    HEX_TRI_123_450 } HEX_TRI_TYPE;
     HEX_TRI_TYPE hex_tri_type;

     // Initialize
     cos_min_angle = 1;
     flag_zero = true;

     compute_hex_triangulation_with_max_min_angle_tri012
       (vcoord[0], vcoord[1], vcoord[2], vcoord[3], vcoord[4], vcoord[5],
        ear_tri012, cos_hex_tri_012, flag_zero_012);

     compute_hex_triangulation_with_max_min_angle_tri012
       (vcoord[5], vcoord[0], vcoord[1], vcoord[2], vcoord[3], vcoord[4],
        ear_tri501, cos_hex_tri_501, flag_zero_501);

     compute_cos_min_hex_triangulation_angle
       (vcoord[3], vcoord[4], vcoord[5], vcoord[0], vcoord[1], vcoord[2],
        max_small_magnitude, cos_hex_tri_v3, flag_zero_v3);

     compute_cos_min_hex_triangulation_angle
       (vcoord[4], vcoord[5], vcoord[0], vcoord[1], vcoord[2], vcoord[3],
        max_small_magnitude, cos_hex_tri_v4, flag_zero_v4);

     compute_hex_triangulation_with_max_min_angle_tri012_tri345
       (vcoord[1], vcoord[2], vcoord[3], vcoord[4], vcoord[5], vcoord[0],
        max_small_magnitude, 
        flag_diag14, cos_hex_tri_123_450, flag_zero_123_450);


     if (!flag_zero_012) {
       cos_min_angle = cos_hex_tri_012;
       flag_zero = false;
       hex_tri_type = HEX_TRI_012;
     }

     if (!flag_zero_501) {
       if (flag_zero) {
         cos_min_angle = cos_hex_tri_501;
         flag_zero = false;
         hex_tri_type = HEX_TRI_501;
       }
       else if (cos_hex_tri_501 < cos_min_angle) {
         cos_min_angle = cos_hex_tri_501;
         hex_tri_type = HEX_TRI_501;
       }
     }

     if (!flag_zero_v3) {
       if (flag_zero) {
         cos_min_angle = cos_hex_tri_v3;
         flag_zero = false;
         hex_tri_type = HEX_TRI_V3;
       }
       else if (cos_hex_tri_v3 < cos_min_angle) {
         cos_min_angle = cos_hex_tri_v3;
         hex_tri_type = HEX_TRI_V3;
       }
     }

     if (!flag_zero_v4) {
       if (flag_zero) {
         cos_min_angle = cos_hex_tri_v4;
         flag_zero = false;
         hex_tri_type = HEX_TRI_V4;
       }
       else if (cos_hex_tri_v4 < cos_min_angle) {
         cos_min_angle = cos_hex_tri_v4;
         hex_tri_type = HEX_TRI_V4;
       }
     }

     if (!flag_zero_123_450) {
       if (flag_zero) {
         cos_min_angle = cos_hex_tri_123_450;
         flag_zero = false;
         hex_tri_type = HEX_TRI_123_450;
       }
       else if (cos_hex_tri_v4 < cos_min_angle) {
         cos_min_angle = cos_hex_tri_123_450;
         hex_tri_type = HEX_TRI_123_450;
       }
     }

     if (flag_zero) { 
       ear[0] = 0;
       ear[1] = 1;
       ear[2] = 2;
       return; 
     }

     if (hex_tri_type == HEX_TRI_012) {
       ear[0] = ear_tri012[0];
       ear[1] = ear_tri012[1];
       ear[2] = ear_tri012[2];
     }
     else if (hex_tri_type == HEX_TRI_501) {
       ear[0] = (ear_tri501[0]+5)%NUM_VERT_PER_HEXAGON;
       ear[1] = (ear_tri501[1]+5)%NUM_VERT_PER_HEXAGON;
       ear[2] = (ear_tri501[2]+5)%NUM_VERT_PER_HEXAGON;
     }
     else if (hex_tri_type == HEX_TRI_V3) {
       ear[0] = 4;
       ear[1] = 5;
       ear[2] = 0;
     }
     else if (hex_tri_type == HEX_TRI_V4) {
       ear[0] = 5;
       ear[1] = 0;
       ear[2] = 1;
     }
     else {
       ear[0] = 2;
       ear[1] = 5;

       if (flag_diag14) { ear[2] = 3; }
       else { ear[2] = 1; }
     }

   }

   ///@}


   // **************************************************
   /// @name TRIANGULATE TWO QUADRILATERALS
   // **************************************************

   ///@{

   /// Triangulate two quadrilaterals sharing an edge.
   /// - Maximize minimum triangle angle.
   /// - Add new triangles to vector tri_vert.
   /// @param[in] quad_vert[] Six vertices shared by the two quadrilaterals.
   /// - Six vertices are listed in clockwise or counter-clockwise order
   ///     around the two quadrilaterals.
   /// - Quadrilaterals share edge (quad_vert[0], quad_vert[3])
   /// @param[in] ivertX Index of possible vertex to use to triangulate
   ///   both quadrilaterals.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_two_quads_max_min_angle
   (const DTYPE dimension,
    const CTYPE * vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivertX,
    const bool flag_reverse_orient,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     const int NUM_VERT_PER_QUAD(4);
     const int NUM_VERT_PER_HEX(6);
     const CTYPE * vcoord0 = vert_coord+quad_vert[0]*dimension;
     const CTYPE * vcoord1 = vert_coord+quad_vert[1]*dimension;
     const CTYPE * vcoord2 = vert_coord+quad_vert[2]*dimension;
     const CTYPE * vcoord3 = vert_coord+quad_vert[3]*dimension;
     const CTYPE * vcoord4 = vert_coord+quad_vert[4]*dimension;
     const CTYPE * vcoord5 = vert_coord+quad_vert[5]*dimension;
     const CTYPE * vcoordX = vert_coord+ivertX*dimension;
     const VTYPE0 quad_vertB[NUM_VERT_PER_QUAD] = 
       { quad_vert[3], quad_vert[4], quad_vert[5], quad_vert[0] };
     CTYPE cos_min_quadA, cos_min_quadB, cos_min_two_quads,cos_min_hex;
     bool flag_diag02_quadA, flag_diag02_quadB;
     bool flag_zero_quadA, flag_zero_quadB;
     bool flag_zero_hex;

     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
        cos_min_quadA, flag_diag02_quadA, flag_zero_quadA);
     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoord3, vcoord4, vcoord5, vcoord0, max_small_magnitude, 
        cos_min_quadB, flag_diag02_quadB, flag_zero_quadB);
     compute_cos_min_hex_tri6_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoord5, 
        vcoordX, max_small_magnitude, cos_min_hex, flag_zero_hex);

     if (cos_min_quadA < cos_min_quadB) 
       { cos_min_two_quads = cos_min_quadB; }
     else 
       { cos_min_two_quads = cos_min_quadA; }

     if (flag_zero_hex || (cos_min_two_quads < cos_min_hex)) {
       // Triangulate quadA and quadB separately using quad diagonals.
       triangulate_quad_using_diagonal
         (quad_vert, flag_diag02_quadA, flag_reverse_orient, tri_vert);
       triangulate_quad_using_diagonal
         (quad_vertB, flag_diag02_quadB, flag_reverse_orient, tri_vert);
     }
     else {
       // Use triangulation into 6 triangles incident on vcoordX.
       triangulate_polygon_with_vertex
         (NUM_VERT_PER_HEX, quad_vert, ivertX, flag_reverse_orient, tri_vert);
     }

   }


   /// Triangulate two quadrilaterals sharing an edge.
   /// - Maximize minimum triangle angle.
   /// - C++ STL vector format for vert_coord.
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, typename VTYPE2,
             typename MTYPE>
   void triangulate_two_quads_max_min_angle
   (const DTYPE dimension,
    const std::vector<CTYPE> & vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivertX,
    const bool flag_reverse_orient,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE2> & tri_vert)
   {
     triangulate_two_quads_max_min_angle
       (dimension, IJK::vector2pointer(vert_coord), quad_vert,
        ivertX, flag_reverse_orient, max_small_magnitude, tri_vert);
   }


   // **************************************************
   /// @name TRIANGULATE THREE QUADRILATERALS
   // **************************************************


   ///@{

   /// Triangulate three quadrilaterals.
   /// - Maximize minimum triangle angle.
   /// - Add new triangles to vector tri_vert.
   /// @param quad_vert[] Eight vertices shared by the three quadrilaterals.
   /// - Quad A has vertices 
   ///     (quad_vert[0], quad_vert[1], quad_vert[6], quad_vert[7])
   /// - Quad B has vertices 
   ///     (quad_vert[1], quad_vert[2], quad_vert[5], quad_vert[6])
   /// - Quad C has vertices 
   ///     (quad_vert[2], quad_vert[3], quad_vert[4], quad_vert[5])
   /// @param ivertX16 Index of vertex splitting (quad_vert[1], quad_vert[6]).
   /// @param ivertX25 Index of vertex splitting (quad_vert[2], quad_vert[5]).
   template <typename DTYPE, typename CTYPE,
             typename VTYPE0, typename VTYPE1, 
             typename VTYPE2, typename VTYPE3,
             typename MTYPE>
   void triangulate_three_quads_max_min_angle
   (const DTYPE dimension,
    const CTYPE * vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivertX16,
    const VTYPE2 ivertX25,
    const bool flag_reverse_orient,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE3> & tri_vert)
   {
     const int NUM_VERT_PER_QUAD(4);
     const int NUM_VERT_PER_HEX(6);
     const CTYPE * vcoord0 = vert_coord+quad_vert[0]*dimension;
     const CTYPE * vcoord1 = vert_coord+quad_vert[1]*dimension;
     const CTYPE * vcoord2 = vert_coord+quad_vert[2]*dimension;
     const CTYPE * vcoord3 = vert_coord+quad_vert[3]*dimension;
     const CTYPE * vcoord4 = vert_coord+quad_vert[4]*dimension;
     const CTYPE * vcoord5 = vert_coord+quad_vert[5]*dimension;
     const CTYPE * vcoord6 = vert_coord+quad_vert[6]*dimension;
     const CTYPE * vcoord7 = vert_coord+quad_vert[7]*dimension;
     const CTYPE * vcoordX16 = vert_coord+ivertX16*dimension;
     const CTYPE * vcoordX25 = vert_coord+ivertX25*dimension;
     const VTYPE0 quadA_vert[NUM_VERT_PER_QUAD] = 
       { quad_vert[0], quad_vert[1], quad_vert[6], quad_vert[7] };
     const VTYPE0 quadB_vert[NUM_VERT_PER_QUAD] = 
       { quad_vert[1], quad_vert[2], quad_vert[5], quad_vert[6] };
     const VTYPE0 quadC_vert[NUM_VERT_PER_QUAD] = 
       { quad_vert[2], quad_vert[3], quad_vert[4], quad_vert[5] };
     CTYPE cos_min_quadA, cos_min_quadB, cos_min_quadC;
     CTYPE cos_min_hexAB, cos_min_hexBC;
     CTYPE cos_min_three_quads, cos_min_hexAB_quadC, cos_min_hexBC_quadA;
     CTYPE cos_min_quadX_lower, cos_min_quadX_upper, cos_min_quadX;
     CTYPE cos_pentagon_triX16, cos_pentagon_triX25;
     bool flag_diag02_quadA, flag_diag02_quadB, flag_diag02_quadC;
     bool flag_diag02_quadX_lower, flag_diag02_quadX_upper;
     bool flag_zero_quadA, flag_zero_quadB, flag_zero_quadC;
     bool flag_zero_quadABC;
     bool flag_zero_hexAB, flag_zero_hexBC;
     bool flag_zero_hexAB_quadC, flag_zero_hexBC_quadA;
     bool flag_zero_quadX_lower, flag_zero_quadX_upper;
     bool flag_zero_pentagon_triX16, flag_zero_pentagon_triX25;
     bool flag_zero_quadX;

     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoord0, vcoord1, vcoord6, vcoord7, max_small_magnitude, 
        cos_min_quadA, flag_diag02_quadA, flag_zero_quadA);
     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoord1, vcoord2, vcoord5, vcoord6, max_small_magnitude, 
        cos_min_quadB, flag_diag02_quadB, flag_zero_quadB);
     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoord2, vcoord3, vcoord4, vcoord5, max_small_magnitude, 
        cos_min_quadC, flag_diag02_quadC, flag_zero_quadC);
     cos_min_three_quads = cos_min_quadA;
     if (cos_min_three_quads < cos_min_quadB)
       { cos_min_three_quads = cos_min_quadB; }
     if (cos_min_three_quads < cos_min_quadC)
       { cos_min_three_quads = cos_min_quadC; }
     flag_zero_quadABC = flag_zero_quadA || flag_zero_quadB || 
       flag_zero_quadC;

     compute_cos_min_hex_tri6_angle
       (dimension, vcoord0, vcoord1, vcoord2, vcoord5, vcoord6, vcoord7, 
        vcoordX16, max_small_magnitude, cos_min_hexAB, flag_zero_hexAB);
     cos_min_hexAB_quadC = cos_min_hexAB;
     if (cos_min_hexAB_quadC < cos_min_quadC)
       { cos_min_hexAB_quadC = cos_min_quadC; }
     flag_zero_hexAB_quadC = flag_zero_hexAB || flag_zero_quadC;

     compute_cos_min_hex_tri6_angle
       (dimension, vcoord1, vcoord2, vcoord3, vcoord4, vcoord5, vcoord6, 
        vcoordX25, max_small_magnitude, cos_min_hexBC, flag_zero_hexBC);
     cos_min_hexBC_quadA = cos_min_hexBC;
     if (cos_min_hexBC_quadA < cos_min_quadA)
       { cos_min_hexBC_quadA = cos_min_quadA; }
     flag_zero_hexBC_quadA = flag_zero_hexBC || flag_zero_quadA;

     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoord1, vcoord2, vcoordX25, vcoordX16, max_small_magnitude, 
        cos_min_quadX_lower, flag_diag02_quadX_lower, 
        flag_zero_quadX_lower);
     compute_cos_min_quad_tri02_tri13_angle
       (dimension, vcoordX16, vcoordX25, vcoord5, vcoord6, max_small_magnitude, 
        cos_min_quadX_upper, flag_diag02_quadX_upper, 
        flag_zero_quadX_upper);

     compute_cos_min_pentagon_triangulation_angle
       (vcoordX16, vcoord6, vcoord7, vcoord0, vcoord1, max_small_magnitude,
        cos_pentagon_triX16, flag_zero_pentagon_triX16);

     compute_cos_min_pentagon_triangulation_angle
      (vcoordX25, vcoord2, vcoord3, vcoord4, vcoord5, max_small_magnitude,
       cos_pentagon_triX25, flag_zero_pentagon_triX25);
 
    cos_min_quadX = cos_min_quadX_lower;
    if (cos_min_quadX < cos_min_quadX_upper)
      { cos_min_quadX = cos_min_quadX_upper; }
    if (cos_min_quadX < cos_pentagon_triX16)
      { cos_min_quadX = cos_pentagon_triX16; }
    if (cos_min_quadX < cos_pentagon_triX25)
      { cos_min_quadX = cos_pentagon_triX25; }
    flag_zero_quadX = flag_zero_quadX_upper || flag_zero_quadX_lower
      || flag_zero_pentagon_triX16 || flag_zero_pentagon_triX25;

    // Default is triangulate quads A, B, C, separately.
    bool flag_tri_quadABC(false);
    bool flag_tri_hexAB_quadC(false), flag_tri_hexBC_quadA(false);
    bool flag_tri_quadX(false);

    if (flag_zero_hexAB_quadC || flag_zero_hexBC_quadA ||
        flag_zero_quadX) { flag_tri_quadABC = true; }
    else {
      if (cos_min_quadX <= cos_min_hexAB_quadC &&
          cos_min_quadX <= cos_min_hexBC_quadA)
        { flag_tri_quadX = true; }
      else if (cos_min_hexAB_quadC < cos_min_hexBC_quadA)
        { flag_tri_hexAB_quadC = true;}
      else
        { flag_tri_hexBC_quadA = true;}

      if (!flag_zero_quadABC) {
        if (cos_min_three_quads <= cos_min_quadX &&
            cos_min_three_quads <= cos_min_hexAB_quadC &&
            cos_min_three_quads <= cos_min_hexBC_quadA) {
          flag_tri_quadX = false;
          flag_tri_hexAB_quadC = false;
          flag_tri_hexBC_quadA = false;
          flag_tri_quadABC = true;
        }
      }
    }

    if (flag_tri_quadABC) {
      // Triangulate quadA, quadB and quadC separately using quad diagonals.
      triangulate_quad_using_diagonal
        (quadA_vert, flag_diag02_quadA, flag_reverse_orient, tri_vert);
      triangulate_quad_using_diagonal
        (quadB_vert, flag_diag02_quadB, flag_reverse_orient, tri_vert);
      triangulate_quad_using_diagonal
        (quadC_vert, flag_diag02_quadC, flag_reverse_orient, tri_vert);
    }
    else if (flag_tri_quadX) {
     const VTYPE0 quadX_lower_vert[NUM_VERT_PER_QUAD] = 
       { quad_vert[1], quad_vert[2], ivertX25, ivertX16 };
     const VTYPE0 quadX_upper_vert[NUM_VERT_PER_QUAD] = 
       { ivertX16, ivertX25, quad_vert[5], quad_vert[6] };

     // Triangulate quadX_lower and quadX_upper
     triangulate_quad_using_diagonal
       (quadX_lower_vert, flag_diag02_quadX_lower, flag_reverse_orient, 
        tri_vert);
     triangulate_quad_using_diagonal
       (quadX_upper_vert, flag_diag02_quadX_upper, flag_reverse_orient, 
        tri_vert);

     // Triangulate quadA with additional vertex at ivertX16.
     triangulate_pentagon
       (ivertX16, quad_vert[6], quad_vert[7], quad_vert[0], quad_vert[1],
        flag_reverse_orient, tri_vert);

     // Triangulate quadC with additional vertex at ivertX25.
     triangulate_pentagon
       (ivertX25, quad_vert[2], quad_vert[3], quad_vert[4], quad_vert[5],
        flag_reverse_orient, tri_vert);
    }
    else if (flag_tri_hexAB_quadC) {
     // Triangulate quadA with additional vertex at ivertX16.
     triangulate_pentagon
       (ivertX16, quad_vert[6], quad_vert[7], quad_vert[0], quad_vert[1],
        flag_reverse_orient, tri_vert);

     // Triangulate quadB with additional vertex at ivertX16.
     triangulate_pentagon
       (ivertX16, quad_vert[1], quad_vert[2], quad_vert[5], quad_vert[6],
        flag_reverse_orient, tri_vert);

     // Triangulate quadC using quad diagonals.
     triangulate_quad_using_diagonal
       (quadC_vert, flag_diag02_quadC, flag_reverse_orient, tri_vert);
    }
    else {
      // Use triangulation hexBC quadA

     // Triangulate quadA using quad diagonals.
      triangulate_quad_using_diagonal
        (quadA_vert, flag_diag02_quadA, flag_reverse_orient, tri_vert);

     // Triangulate quadB with additional vertex at ivertX25.
     triangulate_pentagon
       (ivertX25, quad_vert[5], quad_vert[6], quad_vert[1], quad_vert[2],
        flag_reverse_orient, tri_vert);

     // Triangulate quadC with additional vertex at ivertX25.
     triangulate_pentagon
       (ivertX25, quad_vert[2], quad_vert[3], quad_vert[4], quad_vert[5],
        flag_reverse_orient, tri_vert);
    }
    
  }


  /// Triangulate three quadrilaterals.
  /// - Maximize minimum triangle angle.
  /// - C++ STL vector format for vert_coord.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, 
            typename VTYPE2, typename VTYPE3,
            typename MTYPE>
  void triangulate_three_quads_max_min_angle
   (const DTYPE dimension,
    const std::vector<CTYPE> & vert_coord,
    const VTYPE0 * quad_vert,
    const VTYPE1 ivertX16,
    const VTYPE2 ivertX25,
    const bool flag_reverse_orient,
    const MTYPE max_small_magnitude,
    std::vector<VTYPE3> & tri_vert)
  {
    triangulate_three_quads_max_min_angle
      (dimension, IJK::vector2pointer(vert_coord), quad_vert,
       ivertX16, ivertX25, flag_reverse_orient, max_small_magnitude, 
       tri_vert);
  }



  // **************************************************
  /// @name TRIANGULATE QUAD-TRIANGLE-QUAD
  // **************************************************


  ///@{

  /// Triangulate quad-triangle-quad.
  /// - Maximize minimum triangle angle.
  /// - Add new triangles to vector tri_vert.
  /// @param poly_vert[] Seven vertices shared by the triangle and two quads.
  /// - Quad A has vertices 
  ///     (poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6])
  /// - Triangle B has vertices (poly_vert[1], poly_vert[4], poly_vert[5])
  /// - Quad C has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4])
  /// @param ivertX15 Index of vertex splitting (poly_vert[1], poly_vert[5]).
  /// @param ivertX14 Index of vertex splitting (poly_vert[1], poly_vert[4]).
  ///  - Note: ivertX15 precedes ivertX14 in the argument list because
  ///    ivertX15 is on the boundary of quadA.
  /// @param[out] is_triangleB_split If true, triangleB is split into
  ///    new triangles which are added to tri_vert[].
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, 
            typename VTYPE2, typename VTYPE3,
            typename MTYPE>
  void triangulate_quad_triangle_quad_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vert_coord,
   const VTYPE0 * poly_vert,
   const VTYPE1 ivertX15,
   const VTYPE2 ivertX14,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE3> & tri_vert,
   bool & is_triangleB_split)
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_HEX(6);
    const CTYPE * vcoord0 = vert_coord+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vert_coord+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vert_coord+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vert_coord+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vert_coord+poly_vert[4]*dimension;
    const CTYPE * vcoord5 = vert_coord+poly_vert[5]*dimension;
    const CTYPE * vcoord6 = vert_coord+poly_vert[6]*dimension;
    const CTYPE * vcoordX14 = vert_coord+ivertX14*dimension;
    const CTYPE * vcoordX15 = vert_coord+ivertX15*dimension;
    const VTYPE0 quadA_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6] };
    const VTYPE0 triB_vert[NUM_VERT_PER_TRIANGLE] = 
      { poly_vert[1], poly_vert[4], poly_vert[5] };
    const VTYPE0 quadC_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4] };
    CTYPE cos_min_quadA, cos_min_triB, cos_min_quadC;
    CTYPE cos_min_pentagonAB, cos_min_pentagonBC;
    CTYPE cos_min_no_edge_split;
    CTYPE cos_min_pentagonAB_quadC, cos_min_pentagonBC_quadA;
    CTYPE cos_min_quadX_upper, cos_min_quadX;
    CTYPE cos_pentagon_triX14, cos_pentagon_triX15;
    bool flag_diag02_quadA, flag_diag02_quadC;
    bool flag_diag02_quadX_upper;
    bool flag_zero_quadA, flag_zero_triB, flag_zero_quadC;
    bool flag_zero_no_edge_split;
    bool flag_zero_pentagonAB, flag_zero_pentagonBC;
    bool flag_zero_pentagonAB_quadC, flag_zero_pentagonBC_quadA;
    bool flag_pentagon_quadX_upper;
    bool flag_zero_pentagon_triX14, flag_zero_pentagon_triX15;
    bool flag_zero_quadX_upper, flag_zero_quadX;

    // Initialize
    is_triangleB_split = false;

    compute_cos_min_quad_tri02_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord5, vcoord6, max_small_magnitude, 
       cos_min_quadA, flag_diag02_quadA, flag_zero_quadA);
    compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord4, vcoord5, max_small_magnitude, 
       cos_min_triB, flag_zero_triB);
    compute_cos_min_quad_tri02_tri13_angle
      (dimension, vcoord1, vcoord2, vcoord3, vcoord4, max_small_magnitude, 
       cos_min_quadC, flag_diag02_quadC, flag_zero_quadC);
    cos_min_no_edge_split = cos_min_quadA;
    if (cos_min_no_edge_split < cos_min_triB)
      { cos_min_no_edge_split = cos_min_triB; }
    if (cos_min_no_edge_split < cos_min_quadC)
      { cos_min_no_edge_split = cos_min_quadC; }
    flag_zero_no_edge_split = flag_zero_quadA || flag_zero_triB || 
      flag_zero_quadC;

    compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord4, vcoord5, vcoord6, vcoordX15, 
       max_small_magnitude, cos_min_pentagonAB, flag_zero_pentagonAB);
    cos_min_pentagonAB_quadC = cos_min_pentagonAB;
    if (cos_min_pentagonAB_quadC < cos_min_quadC)
      { cos_min_pentagonAB_quadC = cos_min_quadC; }
    flag_zero_pentagonAB_quadC = flag_zero_pentagonAB || flag_zero_quadC;

    compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord1, vcoord2, vcoord3, vcoord4, vcoord5, vcoordX14, 
       max_small_magnitude, cos_min_pentagonBC, flag_zero_pentagonBC);
    cos_min_pentagonBC_quadA = cos_min_pentagonBC;
    if (cos_min_pentagonBC_quadA < cos_min_quadA)
      { cos_min_pentagonBC_quadA = cos_min_quadA; }
    flag_zero_pentagonBC_quadA = flag_zero_pentagonBC || flag_zero_quadA;

    compute_cos_min_quad_tri02_tri13_angle
      (dimension, vcoordX15, vcoordX14, vcoord4, vcoord5, max_small_magnitude, 
       cos_min_quadX_upper, flag_diag02_quadX_upper, 
       flag_zero_quadX_upper);

    compute_cos_min_pentagon_triangulation_angle
      (vcoordX15, vcoord5, vcoord6, vcoord0, vcoord1, max_small_magnitude,
       cos_pentagon_triX15, flag_zero_pentagon_triX15);

    compute_cos_min_pentagon_triangulation_angle
      (vcoordX14, vcoord1, vcoord2, vcoord3, vcoord4, max_small_magnitude,
       cos_pentagon_triX14, flag_zero_pentagon_triX14);
 
    cos_min_quadX = cos_min_quadX_upper;
    if (cos_min_quadX < cos_pentagon_triX15)
      { cos_min_quadX = cos_pentagon_triX15; }
    if (cos_min_quadX < cos_pentagon_triX14)
      { cos_min_quadX = cos_pentagon_triX14; }
    flag_zero_quadX = flag_zero_quadX_upper|| flag_zero_pentagon_triX15 || 
      flag_zero_pentagon_triX14;

    // Default is do not split edges (1,4) or (1,5).
    bool flag_tri_no_edge_split(false);
    bool flag_tri_pentagonAB_quadC(false), flag_tri_pentagonBC_quadA(false);
    bool flag_tri_quadX(false);

    if (flag_zero_pentagonAB_quadC || flag_zero_pentagonBC_quadA ||
        flag_zero_quadX) { flag_tri_no_edge_split = true; }
    else {
      if (cos_min_quadX <= cos_min_pentagonAB_quadC &&
          cos_min_quadX <= cos_min_pentagonBC_quadA)
        { flag_tri_quadX = true; }
      else if (cos_min_pentagonAB_quadC < cos_min_pentagonBC_quadA)
        { flag_tri_pentagonAB_quadC = true;}
      else
        { flag_tri_pentagonBC_quadA = true;}

      if (!flag_zero_no_edge_split) {
        if (cos_min_no_edge_split <= cos_min_quadX &&
            cos_min_no_edge_split <= cos_min_pentagonAB_quadC &&
            cos_min_no_edge_split <= cos_min_pentagonBC_quadA) {
          flag_tri_quadX = false;
          flag_tri_pentagonAB_quadC = false;
          flag_tri_pentagonBC_quadA = false;
          flag_tri_no_edge_split = true;
        }
      }
    }

    if (flag_tri_no_edge_split) {
      // Triangulate quadA and quadC separately using quad diagonals.
      triangulate_quad_using_diagonal
        (quadA_vert, flag_diag02_quadA, flag_reverse_orient, tri_vert);
      triangulate_quad_using_diagonal
        (quadC_vert, flag_diag02_quadC, flag_reverse_orient, tri_vert);
    }
    else if (flag_tri_quadX) {
      const VTYPE0 quadX_upper_vert[NUM_VERT_PER_QUAD] = 
        { ivertX15, ivertX14, poly_vert[4], poly_vert[5] };

      // Triangulate quadX_upper
      triangulate_quad_using_diagonal
        (quadX_upper_vert, flag_diag02_quadX_upper, flag_reverse_orient, 
         tri_vert);

      // Add triangle below quadX_upper
      add_triangle_vertices
        (poly_vert[1], ivertX14, ivertX15, flag_reverse_orient, tri_vert);

      // Triangulate quadA with additional vertex at ivertX16.
      triangulate_pentagon
        (ivertX15, poly_vert[5], poly_vert[6], poly_vert[0], poly_vert[1],
         flag_reverse_orient, tri_vert);

      // Triangulate quadC with additional vertex at ivertX14.
      triangulate_pentagon
        (ivertX14, poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4],
         flag_reverse_orient, tri_vert);

      is_triangleB_split = true;
    }
    else if (flag_tri_pentagonAB_quadC) {

      // Triangulate pentagon AB with additional vertex at ivertX15.
      triangulate_pentagon_with_vertex
        (poly_vert[0], poly_vert[1], poly_vert[4], poly_vert[5], poly_vert[6],
         ivertX15, flag_reverse_orient, tri_vert);

      // Triangulate quadC using quad diagonals.
      triangulate_quad_using_diagonal
        (quadC_vert, flag_diag02_quadC, flag_reverse_orient, tri_vert);

      is_triangleB_split = true;
    }
    else {
      // Use triangulation pentagonBC quadA

      // Triangulate quadA using quad diagonals.
      triangulate_quad_using_diagonal
        (quadA_vert, flag_diag02_quadA, flag_reverse_orient, tri_vert);

      // Triangulate pentagonBC with additional vertex at ivertX14.
      triangulate_pentagon_with_vertex
        (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4], poly_vert[5],
         flag_reverse_orient, tri_vert);

      is_triangleB_split = true;
    }
    
  }


  /// Triangulate quad, triangle, quad.
  /// - Maximize minimum triangle angle.
  /// - C++ STL vector format for vert_coord.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, 
            typename VTYPE2, typename VTYPE3,
            typename MTYPE>
  void triangulate_quad_triangle_quad_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vert_coord,
   const VTYPE0 * poly_vert,
   const VTYPE1 ivertX15,
   const VTYPE2 ivertX14,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE3> & tri_vert,
   bool & is_triangleB_split)
  {
    triangulate_quad_triangle_quad_max_min_angle
      (dimension, IJK::vector2pointer(vert_coord), poly_vert,
       ivertX15, ivertX14, flag_reverse_orient, max_small_magnitude, 
       tri_vert, is_triangleB_split);
  }


  // **************************************************
  //! @name ENVELOPE
  // **************************************************

  ///@{

  /// Return true if diagonal (w0,w2) is in envelope.
  /// @param v0[] Grid vertex v0.
  /// @param v1[] Grid vertex v1. (v0,v1) is a grid edge.
  /// @param w0[] Isosurface vertex w0.
  /// @param w1[] Isosurface vertex w1.
  /// @param w2[] Isosurface vertex w2.
  /// @param w3[] Isosurface vertex w3.
  ///        (w0,w1,w2,w3) is an isosurface quadrilateral q dual to (v0,v1).
  ///        Vertices (w0,w1,w2,w3) are listed in order around q.
  /// @pre  Orientation of (w0,w1,v0,v1) is positive.
  /// @param epsilon Tolerance.  
  ///        Determinants less than epsilon are considered equivalent to 0.
  template <typename VCOORD_TYPE, typename WCOORD_TYPE,
            typename EPSILON_TYPE>
  bool is_in_envelope_3D
  (const VCOORD_TYPE v0[3], const VCOORD_TYPE v1[3],
   const WCOORD_TYPE w0[3], const WCOORD_TYPE w1[3],
   const WCOORD_TYPE w2[3], const WCOORD_TYPE w3[3],
   const EPSILON_TYPE epsilon)
  {
    double D;

    IJK::determinant_point_3D(w0, w2, v1, w1, D);
    if (D < epsilon) { return(false); }
    IJK::determinant_point_3D(w0, w2, v0, w1, D);
    if (-D < epsilon) { return(false); }
    IJK::determinant_point_3D(w0, w2, v1, w3, D);
    if (-D < epsilon) { return(false); }
    IJK::determinant_point_3D(w0, w2, v0, w3, D);
    if (D < epsilon) { return(false); }

    return(true);
  }

  ///@}

}

#endif
