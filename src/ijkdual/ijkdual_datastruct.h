/// @file ijkdual_datastruct.h
/// @brief Data structure definitions for ijkdual.
/// - Version 0.3.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2026 Rephael Wenger

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

#ifndef _IJKDUAL_DATASTRUCT_H_
#define _IJKDUAL_DATASTRUCT_H_

#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkisopoly.tpp"
#include "ijkscalar_grid.tpp"
#include "ijkmerge.tpp"


#include "ijkdualtable.tpp"

#include "ijkdual_types.h"
#include "ijkdual_constants.h"

namespace IJKDUAL {

  class MERGE_DATA;
  class MULTIRES_GRID;


  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef IJK::GRID_PLUS
  <GRID_CUBEV_BITSET_SIZE, int, AXIS_SIZE_TYPE, VERTEX_INDEX, VERTEX_INDEX> 
    GRID_PLUS;                     ///< Regular grid.
  typedef IJK::GRID_SPACING<COORD_TYPE, GRID_PLUS>
    DUALISO_GRID;                  ///< Grid and spacing information.
  typedef IJK::SCALAR_GRID_BASE<DUALISO_GRID, SCALAR_TYPE> 
    DUALISO_SCALAR_GRID_BASE;      ///< Marching Cubes base scalar grid.
  typedef IJK::SCALAR_GRID_WRAPPER<DUALISO_GRID, SCALAR_TYPE>
    DUALISO_SCALAR_GRID_WRAPPER;   ///< Marching Cubes scalar grid wrapper.
  typedef IJK::SCALAR_GRID<DUALISO_GRID, SCALAR_TYPE> 
    DUALISO_SCALAR_GRID;           ///< Marching Cubes scalar grid.


  // **************************************************
  // ISOSURFACE LOOKUP TABLES
  // **************************************************

  typedef IJKDUALTABLE::ISODUAL_TABLE_ENTRY<int,int,int> ISODUAL_TABLE_ENTRY;

  typedef IJKDUALTABLE::ISODUAL_AMBIG_TABLE_ENTRY<int,int,int,FACET_BITS_TYPE> 
  ISODUAL_AMBIG_TABLE_ENTRY;

  typedef IJKDUALTABLE::ISODUAL_CUBE_TABLE
  <int,int,TABLE_INDEX,ISODUAL_TABLE_ENTRY>
  ISODUAL_CUBE_TABLE;

  typedef IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG
  <int,int,TABLE_INDEX,ISODUAL_AMBIG_TABLE_ENTRY>
  ISODUAL_CUBE_TABLE_AMBIG;


  // **************************************************
  // CUBE DATA STRUCTURES
  // **************************************************

  typedef IJK::CUBE_FACE_INFO<int, int, int>
    DUALISO_CUBE_FACE_INFO;


  // **************************************************
  // DUAL ISOSURFACE VERTICES
  // **************************************************

  typedef IJK::DUAL_ISOVERT
  <ISO_VERTEX_INDEX, FACET_VERTEX_INDEX, TABLE_INDEX,
   ISO_VERTEX_INDEX>
  DUAL_ISOVERT;

  typedef IJK::DUAL_ISOVERT_BFLAG
  <ISO_VERTEX_INDEX, FACET_VERTEX_INDEX, TABLE_INDEX,
   ISO_VERTEX_INDEX>
  DUAL_ISOVERT_BFLAG;


  // **************************************************
  // GRID EDGE
  // **************************************************

  typedef IJK::GRID_EDGE<VERTEX_INDEX, DIRECTION_TYPE> GRID_EDGE_TYPE;

  /// Array of grid edges.
  typedef std::vector<GRID_EDGE_TYPE> GRID_EDGE_ARRAY;


  // ***************************************************************
  // ACTIVE GRID CUBE
  // ***************************************************************  

  typedef IJK::ACTIVE_GRID_CUBE_X_T<ISO_VERTEX_INDEX, TABLE_INDEX, 
                                    FACET_VERTEX_INDEX, ISO_VERTEX_INDEX>
  GRID_CUBE_DATA;
  

  // ***************************************************************
  // DUAL CONTOURING ISOSURFACE CLASS
  // ***************************************************************

  //! Dual contouring isosurface.
  //! Representation of isosurface returned by Dual Contouring algorithm.
  //! @tparam ISOPOLY_INFO_TYPE Class storing isosurface polytope information.
  //!   Should be class derived from IJK::GRID_EDGE.
  template <typename ISOPOLY_INFO_TYPE>
  class DUAL_ISOSURFACE_BASE {

  protected:
    int dimension;
    VERTEX_INDEX numv_per_isopoly;

    void Init(const int dimension, const VERTEX_INDEX numv_per_isopoly);

  public:
    /// List of isosurface polytope vertices.
    VERTEX_INDEX_ARRAY isopoly_vert;

    std::vector<ISOPOLY_INFO_TYPE> isopoly_info;

    /// List of vertex coordinates.
    COORD_ARRAY vertex_coord;

    /// Index of first isosurface vertex dual to an isosurface polytope.
    /// - i'th isosurface vertex after first_isov_dual_to_iso_poly
    ///   is dual to i'th isosurface polytope.
    ISO_VERTEX_INDEX first_isov_dual_to_iso_poly;

    /// List of vertices of each triangle in triangle mesh.
    /// Only for 2D surface embedded in 3D.
    /// Not always set.
    VERTEX_INDEX_ARRAY tri_vert;

    /// cube_containing_isopoly[i] = cube containing i'th isopoly.
    /// Not always set.
    VERTEX_INDEX_ARRAY cube_containing_isopoly;

  public:
    DUAL_ISOSURFACE_BASE
      (const int dimension, const VERTEX_INDEX numv_per_isopoly)
      { Init(dimension, numv_per_isopoly); }

    int Dimension() const
    { return(dimension); }

    VERTEX_INDEX NumVerticesPerIsoPoly() const
      { return(numv_per_isopoly); };

    VERTEX_INDEX NumIsoPoly() const
      { return(isopoly_vert.size()/NumVerticesPerIsoPoly()); };

    VERTEX_INDEX NumIsoVert() const
    { return(vertex_coord.size()/dimension); }

    void Clear();

    /// Return false and set error message in error if isopoly_info.size() 
    ///   does not match number of poly.
    bool CheckIsopolyInfo(IJK::ERROR & error) const;
  };

  typedef DUAL_ISOSURFACE_BASE<GRID_EDGE_TYPE> DUAL_ISOSURFACE;



  // **************************************************
  // DUAL CONTOURING INPUT DATA AND DATA STRUCTURES
  // **************************************************

  /// dual contouring flags.
  class DUALISO_DATA_FLAGS {

  protected:
    void Init();

  public:
    // control parameters
    INTERPOLATION_TYPE interpolation_type;
    VERTEX_POSITION_METHOD vertex_position_method;
    bool use_triangle_mesh;
    QUAD_TRI_METHOD quad_tri_method;
    QUAD_EDGE_INTERSECTION_METHOD quad_edge_intersection_method;

    /// If true, constructing an interval volume.
    bool flag_interval_volume;

    /// If true, allow triangulation of hexagons into 4 triangles
    ///  by adding an extra vertex.
    bool flag_tri4_quad;

    /// If true, allow triangulation of pentagons into 5 triangles
    ///   by adding an extra vertex.
    /// Used when isosurface construction combines quad and triangles
    ///   into pentagons.
    bool flag_tri5_pentagon;

    /// Position method for additional vertex in quadrilateral triangulations.
    TRI4_POSITION_METHOD tri4_position_method;

    /// If true, allow multiple isosurface vertices in a single cube
    ///   to represent multiple isosurface patches.
    bool allow_multiple_iso_vertices;

    /// If true, split non-manifold isosurface vertex pairs.
    /// Isosurface will be a manifold is isovalue does not
    ///   equal any grid scalar value.
    bool flag_split_non_manifold;

    /// If true, select the cube with split isosurface patch,
    /// when adjacent cubes share an ambiguous facet and one will have
    /// one isosurface vertex while the other has two isosurface vertices.
    bool flag_select_split;

    /// If true, select configurations of ambiguous cubes to increase
    /// connectivity of the isosurface.
    bool flag_connect_ambiguous;

    /// If true, isosurfaced patches separate negative vertices.
    bool flag_separate_neg;

    /// If true, all isosurface quadrilaterals are dual to grid edges.
    bool flag_iso_quad_dual_to_grid_edges;

    /// If true, call dual_collapse to merge close isosurface vertices.
    bool flag_dual_collapse;

    /// Parameter for selecting triangulation based on distance to facets.
    /// If all vertices of quadrilateral q are distance at least
    ///   min_distance_use_tri4 to the facets incident on the grid edge
    ///   dual to q, then split q into 4 triangles.
    /// Split q into 4 triangles by adding a vertex on the grid edge
    ///   dual to the triangle.
    COORD_TYPE min_distance_use_tri4;

    /// Parameter for selecting triangulation based on distance to facets.
    /// If all vertices of quadrilateral q are distance at least
    ///   min_distance_allow_tri4 to the facets incident on the grid edge
    ///   dual to q, then consider triangulation of q into 4 triangles.
    /// Split q into 4 triangles by adding a vertex on the grid edge
    ///   dual to the triangle.
    COORD_TYPE min_distance_allow_tri4;

    /// Vectors with magnitude less than or equal to max_small_magnitude
    ///   are treated as zero vectors.
    COORD_TYPE max_small_magnitude;

  public:
    DUALISO_DATA_FLAGS() { Init(); };
    ~DUALISO_DATA_FLAGS() { Init(); };

    /// Set
    void Set(const DUALISO_DATA_FLAGS & data_flags);

    // Get functions
    bool AllowMultipleIsoVertices() const
      { return(allow_multiple_iso_vertices); }
    bool CollapseFlag() const
      { return(flag_dual_collapse); }
    bool SplitNonManifoldFlag() const
      { return(flag_split_non_manifold); }
    bool SelectSplitFlag() const
      { return(flag_select_split); }
    bool ConnectAmbiguousFlag() const
      { return(flag_connect_ambiguous); }
    bool SeparateNegFlag() const
      { return(flag_separate_neg); }
    COORD_TYPE MaxSmallMagnitude() const
      { return(max_small_magnitude); }
    bool UseTriangleMesh() const
      { return(use_triangle_mesh); }
    QUAD_TRI_METHOD QuadTriangulationMethod() const
      { return(quad_tri_method); }
    COORD_TYPE MinDistanceUseTri4() const
      { return(min_distance_use_tri4); }
    COORD_TYPE MinDistanceAllowTri4() const
      { return(min_distance_allow_tri4); }

    /// Return interpolation type.
    INTERPOLATION_TYPE InterpolationType() const
      { return(interpolation_type); };

    /// Return isosurface vertex position method.
    VERTEX_POSITION_METHOD VertexPositionMethod() const
      { return(vertex_position_method); };

    /// Return false if flag_iso_quad_dual_to_grid_edges is false.
    bool CheckIsoQuadDualToGridEdges
    (const char * proc_name, IJK::ERROR & error) const;
  };


  /// Input data to Dual Contouring and related algorithms
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  class DUALISO_DATA_BASE:public DATA_FLAGS_TYPE {

  protected:
    SCALAR_GRID_TYPE scalar_grid;

    // flags
    bool is_scalar_grid_set;

    void Init();
    void FreeAll();

  public:
    DUALISO_DATA_BASE() { Init(); };
    ~DUALISO_DATA_BASE() { FreeAll(); };

    // Set functions
    template <typename SCALAR_GRID_BASE_TYPE>
    void CopyScalarGrid             /// Copy scalar_grid to DUALISO_DATA
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2);
    template <typename SCALAR_GRID_BASE_TYPE>
    void SubsampleScalarGrid        /// Subsample scalar_grid.
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
       const int subsample_resolution);
    template <typename SCALAR_GRID_BASE_TYPE>
    void SupersampleScalarGrid      /// Supersample scalar_grid.
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
       const int supersample_resolution);
    void SetInterpolationType       /// Set type of interpolation.
      (const INTERPOLATION_TYPE interpolation_type);
    void SetVertexPositionMethod    /// Set isosurface vertex position method.
      (const VERTEX_POSITION_METHOD vertex_position_method);
    void SetAllowMultipleIsoVertices /// Set flag for multiple iso vertices.
      (const bool flag);
    void SetConvertQuadToTri
      (const QUAD_TRI_METHOD quad_tri_method);
    void UnsetConvertQuadToTri();

    /// Copy, subsample or supersample scalar grid.
    /// Precondition: flag_subsample and flag_supersample are not both true.
    template <typename SCALAR_GRID_BASE_TYPE>
    void SetScalarGrid
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
       const bool flag_subsample, const int subsample_resolution,
       const bool flag_supersample, const int supersample_resolution);

    // Get functions
    bool IsScalarGridSet() const     /// Return true if scalar grid is set.
      { return(is_scalar_grid_set); };
    const DUALISO_SCALAR_GRID_BASE & ScalarGrid() const /// Return scalar_grid
      { return(scalar_grid); };

    /// Check data structure.
    /// Return true if no errors found.
    bool Check(IJK::ERROR & error) const;
  };

  typedef DUALISO_DATA_BASE<DUALISO_SCALAR_GRID, DUALISO_DATA_FLAGS> 
  DUALISO_DATA;


  // **************************************************
  // DUALISO TIME
  // **************************************************

  /// dual contouring time.
  /// Uses system clock() function to determine time.  
  /// System clock() function should return cpu time, 
  ///   but may return wall time.
  class DUALISO_TIME {

  public:
    // all times are in seconds

    /// time to create data structure for faster isosurface extraction
    float preprocessing;  

    float extract;    ///< time to extract isosurface mesh
    float merge;      ///< time to merge identical vertices
    float position;   ///< time to position isosurface vertices
    float collapse;   ///< time to collapse isosurface vertices
    float total;      /// extract_time+merge_time+position_time+collapse_time

    DUALISO_TIME();
    void Clear();
    void Add(const DUALISO_TIME & dualiso_time);
  };

  // **************************************************
  // GRID INFO
  // **************************************************

  /// Regular grid information.
  class GRID_INFO {

  public:
    GRID_INFO();                    ///< Constructor.

    VERTEX_INDEX num_cubes;         ///< Number of grid cubes.

    void Clear();                   ///< Clear all data.
  };

  // **************************************************
  // SCALAR INFO
  // **************************************************

  /// Scalar grid information.
  class SCALAR_INFO {

  protected:
    int dimension;

    void Copy(const SCALAR_INFO & info);  ///< Copy scalar info.
    void Init(const int dimension);       ///< Initialize scalar info.
    void FreeAll();                       ///< Free all memory.

  public:
    SCALAR_INFO() { Init(3); };
    SCALAR_INFO(const int dimension) { Init(dimension); };
    ~SCALAR_INFO();
    SCALAR_INFO(const SCALAR_INFO & info) ///< Copy constructor. 
      { Copy(info); };       
    const SCALAR_INFO & operator =        ///< Copy assignment.
      (const SCALAR_INFO & right);

    VERTEX_INDEX num_non_empty_cubes;     ///< Number of cubes containing an isosurface simplex.

    VERTEX_INDEX num_bipolar_edges; ///< Number of bipolar edges.
    //   (number of edges with some vertex value less than the isovalue
    //    and some vertex value greater than or equal to the isovalue)


    // Set functions.
    void SetDimension(const int dimension); ///< Set dimension.

    // Get functions.
    int Dimension() const { return(dimension); };

    void Clear();     // clear all data
  };


  // **************************************************
  // MULTI ISOV INFO
  // **************************************************

  /// Information on mutiple isosurface vertices per cube
  class MULTI_ISOV_INFO {

  public:
    MULTI_ISOV_INFO();

    VERTEX_INDEX num_cubes_single_isov;
    VERTEX_INDEX num_cubes_multi_isov;
    VERTEX_INDEX num_non_manifold_split;
    VERTEX_INDEX num_1_2_changed;
    VERTEX_INDEX num_connect_changed;
    VERTEX_INDEX num_ambig_ridge_cubes_changed;
    VERTEX_INDEX num_non_ambig_ridge_cubes_changed;

    
    void Clear();     // clear all data
  };


  // **************************************************
  // DUALISO INFO
  // **************************************************

  /// dual contouring information.
  /// Statistical and timing information from the Dual Contouring algorithm.
  class DUALISO_INFO {

  public:
    GRID_INFO grid;
    SCALAR_INFO scalar;
    DUALISO_TIME time;
    MULTI_ISOV_INFO multi_isov;

    DUALISO_INFO();
    DUALISO_INFO(const int dimension);

    void Clear();     // clear all data
  };

  // **************************************************
  // MERGE DATA
  // **************************************************

  /// Internal data structure for merge_identical_vertices 
  class MERGE_DATA: public IJK::INTEGER_LIST<MERGE_INDEX, MERGE_INDEX> {

  protected:
    MERGE_INDEX num_edges;             ///< Number of edges.
    MERGE_INDEX num_vertices;          ///< Number of vertices.
    MERGE_INDEX num_obj_per_vertex;    ///< Number of objects per vertex.
    MERGE_INDEX num_obj_per_edge;      ///< Number of objects per edge.
    MERGE_INDEX num_obj_per_grid_vertex; ///< Number of objects per grid vertex.
    MERGE_INDEX vertex_id0;            ///< First vertex identifier.

    /// Initialize.
    void Init(const int dimension, const AXIS_SIZE_TYPE * axis_size,
              const MERGE_INDEX num_obj_per_vertex,
              const MERGE_INDEX num_obj_per_edge);

  public:
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size)
      { Init(dimension, axis_size, 0, 1); };
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size,
               const MERGE_INDEX num_obj_per_vertex, 
               const MERGE_INDEX num_obj_per_edge)
      { Init(dimension, axis_size, num_obj_per_vertex, num_obj_per_edge); };

    // get functions
    MERGE_INDEX NumEdges() const        /// Number of edges.
      { return(num_edges); };
    MERGE_INDEX NumVertices() const     /// Number of vertices.
      { return(num_vertices); };
    MERGE_INDEX NumObjPerVertex() const { return(num_obj_per_vertex); };
    MERGE_INDEX NumObjPerEdge() const { return(num_obj_per_edge); };
    MERGE_INDEX NumObjPerGridVertex() const
      { return(num_obj_per_grid_vertex); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv) const { return(vertex_id0 + iv); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv, const MERGE_INDEX j) const 
      { return(vertex_id0 + iv + j * num_vertices); };
    MERGE_INDEX EdgeIdentifier         /// Edge identifier.
      (const MERGE_INDEX ie) const { return(ie); };
    MERGE_INDEX EdgeIdentifier         /// edge identifier
      (const MERGE_INDEX ie, const MERGE_INDEX j) const 
      { return(ie + j * num_edges); };

    /// Get first endpoint of edge containing isosurface vertex isov.
    inline VERTEX_INDEX GetFirstEndpoint(const MERGE_INDEX isov) const
      { return(isov/NumObjPerGridVertex()); };

    /// Get direction of edge containing isosurface vertex isov.
    inline MERGE_INDEX GetEdgeDir(const MERGE_INDEX isov) const
      { return(isov%NumObjPerGridVertex()); };

    bool Check(IJK::ERROR & error) const;     ///< Check allocated memory.
  };

  /// Merge data structure for isosurface vertices.
  /// Vertices are identified by a single integer.
  class ISO_MERGE_DATA: public MERGE_DATA {
  public:
    ISO_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 1, 0) {};
  };

  // **************************************************
  // MESH VERTEX LIST
  // **************************************************

  /// List of mesh vertices.
  template <class DTYPE, class CTYPE, class STYPE>
    class MESH_VERTEX_LIST {

    protected:
    DTYPE dimension;        ///< Dimension.

    /// Coordinates of mesh vertices.
    /// coord[i*dimension+j] = j'th coordinate of i'th mesh vertex.
    std::vector<CTYPE> coord;
    std::vector<STYPE> scalar;   ///< scalar[i] = scalar value of i'th mesh vertex.

    public:
    MESH_VERTEX_LIST(const DTYPE dim) { this->dimension = dim; };

    // Get functions.
    DTYPE Dimension() const { return(dimension); } ///< Return dimension.

    /// Return j'th coordinate of i'th mesh vertex.
    CTYPE Coord(const int i, const int j) const
      { return(coord[i*Dimension()+j]); };

    /// Copy coordinates of \a i'th mesh vertex into array \a coord[].
    /// @pre coord[] is preallocated to length at least \a dimension.
    template <class CTYPE2>
      void CopyCoord(const int i, CTYPE2 * coord2) const
      {
        const CTYPE * coord_begin = &(coord[i*Dimension()]);
        std::copy(coord_begin, coord_begin+Dimension(), coord2);
      }

    /// Return scalar value of i'th mesh vertex.
    STYPE Scalar(const int i) const { return(scalar[i]); };

    int NumVertices() const     ///< Return number of vertices in list.
      { return(scalar.size()); };

    // Set functions.

    void Add()                 ///< Add a vertex to the list.
      {
        const int k = coord.size();
        coord.resize(k+Dimension());
        scalar.push_back(0);
      }

    void SetScalar(const int i, const STYPE s)  ///< Set scalar of vertex i.
      { scalar[i] = s; };

    /// Set j'th coordinate of vertex i.
    void SetCoord(const int i, const int j, const CTYPE c)
      { coord[i*Dimension()+j] = c; };

    /// Set coordinates of vertex i
    template <class CTYPE2>
      void SetCoord(const int i, CTYPE2 * coord2)
      {
        for (DTYPE d = 0; d < Dimension(); d++)
          { SetCoord(i, d, coord2[d]); };
      }

  };

  typedef MESH_VERTEX_LIST<int, COORD_TYPE, SCALAR_TYPE> 
    DUALISO_MESH_VERTEX_LIST;    ///< Dual Contouring mesh vertex list.


  // *******************************************************************
  // DUAL_ISOSURFACE_BASE CLASS MEMBER FUNCTIONS
  // *******************************************************************

  template <typename ISOPOLY_INFO_TYPE>
  void DUAL_ISOSURFACE_BASE<ISOPOLY_INFO_TYPE>::Init
  (const int dimension, const VERTEX_INDEX numv_per_isopoly)
  {
    this->dimension = dimension;
    this->numv_per_isopoly = numv_per_isopoly;
    first_isov_dual_to_iso_poly = 0;
  }


  template <typename ISOPOLY_INFO_TYPE>
  void DUAL_ISOSURFACE_BASE<ISOPOLY_INFO_TYPE>::Clear()
  {
    isopoly_vert.clear();
    vertex_coord.clear();
    first_isov_dual_to_iso_poly = 0;
  }


  template <typename ISOPOLY_INFO_TYPE>
  bool DUAL_ISOSURFACE_BASE<ISOPOLY_INFO_TYPE>::
  CheckIsopolyInfo(IJK::ERROR & error) const
  {
    const int numv_per_iso_poly = NumVerticesPerIsoPoly();
    const int num_poly = isopoly_vert.size()/numv_per_iso_poly;

    if (isopoly_info.size() != num_poly) {
      if (isopoly_info.size() == 0) {
        error.AddMessage
          ("Programming error. Array DUAL_ISOSURFACE_BASE::isopoly_info not allocated.");
        return(false);
      }
      else {
        error.AddMessage
          ("Programming error. Size of array DUAL_ISOSURFACE_BASE::isopoly_info");
        error.AddMessage
          ("  does not match number of isosurface polytopes stored in isopoly_vert.");
        error.AddMessage("  Number of isosurface polytopes: ", num_poly, "");
        error.AddMessage
          ("  Size of array DUAL_ISOSURFACE_BASE::isopoly_info: ", 
           isopoly_info.size(), "");
        return(false);
      }
    }

    return(true);
  }

  // *******************************************************************
  // DUALISO_DATA_BASE CLASS MEMBER FUNCTIONS
  // *******************************************************************

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::Init()
  {
    is_scalar_grid_set = false;
  }

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE, DATA_FLAGS_TYPE>::FreeAll()
  {
    is_scalar_grid_set = false;
  }


  // Copy scalar grid
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE, DATA_FLAGS_TYPE>::
  CopyScalarGrid(const SCALAR_GRID_BASE_TYPE & scalar_grid2)
  {
    scalar_grid.Copy(scalar_grid2);
    scalar_grid.SetSpacing(scalar_grid2.SpacingPtrConst());
    is_scalar_grid_set = true;
  }

  // Subsample scalar grid
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SubsampleScalarGrid(const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
                      const int subsample_resolution)
  {
    const int dimension = scalar_grid2.Dimension();
    IJK::ARRAY<COORD_TYPE> spacing(dimension);

    scalar_grid.Subsample(scalar_grid2, subsample_resolution);

    IJK::copy_coord(dimension, scalar_grid2.SpacingPtrConst(),
                    spacing.Ptr());
    IJK::multiply_coord
      (dimension, subsample_resolution, spacing.PtrConst(), spacing.Ptr());
    scalar_grid.SetSpacing(spacing.PtrConst());

    is_scalar_grid_set = true;
  }

  // Supersample scalar grid
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SupersampleScalarGrid
  (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
   const int supersample_resolution)
  {
    const int dimension = scalar_grid2.Dimension();
    IJK::ARRAY<COORD_TYPE> spacing(dimension);
    scalar_grid.Supersample(scalar_grid2, supersample_resolution);

    IJK::copy_coord(dimension, scalar_grid2.SpacingPtrConst(),
                    spacing.Ptr());
    IJK::divide_coord
      (dimension, supersample_resolution, spacing.PtrConst(), spacing.Ptr());
    scalar_grid.SetSpacing(spacing.PtrConst());

    is_scalar_grid_set = true;
  }

  // Copy, subsample or supersample scalar grid.
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::SetScalarGrid
  (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
   const bool flag_subsample, const int subsample_resolution,
   const bool flag_supersample, const int supersample_resolution)
  {
    IJK::PROCEDURE_ERROR error("DUALISO_DATA_BASE::SetScalarGrid");

    if (flag_subsample && flag_supersample) {
      error.AddMessage
        ("Scalar grid cannot both be subsampled and supersampled.");
      throw error;
    }
  
    if (flag_subsample) {
      // subsample grid
      SubsampleScalarGrid(scalar_grid2, subsample_resolution);
    }
    else if (flag_supersample) {
      // supersample grid
      SupersampleScalarGrid(scalar_grid2, supersample_resolution);
    }
    else {
      CopyScalarGrid(scalar_grid2);
    };
  }

  // Set type of interpolation
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetInterpolationType
  (const INTERPOLATION_TYPE interpolation_type)
  {
    this->interpolation_type = interpolation_type;
  }

  // Set isosurface vertex position method.
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetVertexPositionMethod    
  (const VERTEX_POSITION_METHOD vertex_position_method)
  {
    this->vertex_position_method = vertex_position_method;
  }

  // Set flag for allowing multiple isosurface vertices in a single cube.
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetAllowMultipleIsoVertices(const bool flag)
  {
    this->allow_multiple_iso_vertices = flag;
  }

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetConvertQuadToTri(const QUAD_TRI_METHOD method)
  {
    this->use_triangle_mesh = true;
    this->quad_tri_method = method;
  }

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  UnsetConvertQuadToTri()
  {
    this->use_triangle_mesh = false;
  }


  /// Check data structure
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  bool DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  Check(IJK::ERROR & error) const
  {
    IJK::ERROR error2;

    if (!IsScalarGridSet()) {
      error.AddMessage("Scalar grid is not set.");
      return(false);
    }

    return(true);
  }

}

#endif
