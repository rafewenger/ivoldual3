/// @file ivoldual_datastruct.h
/// @brief Data structure definitions for ivoldual.
/// - Copies types from ijkdual.h.
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

#ifndef _IVOLDUAL_DATASTRUCT_
#define _IVOLDUAL_DATASTRUCT_

#include "ijkcube.tpp"
#include "ijkmesh_datastruct.tpp"
#include "ijkdual_mesh.tpp"

#include "ijkdual_datastruct.h"
#include "ivoldual_types.h"
#include "ivoldualtable.h"


namespace IVOLDUAL {

  // **************************************************
  // INTERVAL VOLUME GRID VERTEX ENCODING
  // **************************************************

  /// Type for encoding of grid vertices.
  /// - Encoded vertices have scalar values 0,1,2 or 3.
  typedef unsigned char GRID_VERTEX_ENCODING;

  // **************************************************
  // GRID DATA STRUCTURES
  // **************************************************

  typedef IJKDUAL::DUALISO_GRID DUALISO_GRID;
  typedef IJKDUAL::DUALISO_SCALAR_GRID_BASE DUALISO_SCALAR_GRID_BASE;
  typedef IJKDUAL::DUALISO_SCALAR_GRID DUALISO_SCALAR_GRID;

  class IVOLDUAL_SCALAR_GRID: public DUALISO_SCALAR_GRID {

  public:
    IVOLDUAL_SCALAR_GRID() {};

    // Subdivide \a scalar grid.
    void Subdivide
    (const DUALISO_SCALAR_GRID_BASE & scalar_grid2, const int subdivide_period, 
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    // Interpolate between subdivide vertices.
    void SubdivideInterpolate
    (const int subdivide_period, const SCALAR_TYPE isovalue0, 
     const SCALAR_TYPE isovalue1); 

    int EliminateAmbiguity
      (const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    void EliminateDiagonalPlus
      (const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1, 
       VTYPE icube, int caseID, int& num_change);

    void EliminateDiagonalMinus
      (const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1, 
       VTYPE icube, int caseID, int& num_change);

    void AddOuterLayer();
  };

  /// Type of grid encoding grid vertices.
  /// - 0: Below lower isovalue.
  /// - 1: Between lower and upper isovalue.
  /// - 2: Between lower and upper isovalue.
  /// - 3: Above upper isovalue.
  typedef IJK::SCALAR_GRID<DUALISO_GRID, GRID_VERTEX_ENCODING> 
  IVOLDUAL_ENCODED_GRID;

  /// Type of cube with face info.
  typedef IJK::CUBE_FACE_INFO<int,int,int> CUBE_FACE_INFO;

  /// Directions of cube center to cube vertices.
  typedef IJK::CUBE_CENTER_TO_VERTEX_DIRECTIONS_3D<int,int,COORD_TYPE>
  CUBE_DIRECTIONS;


  // **************************************************
  // DUAL INTERVAL VOLUME VERTICES
  // **************************************************

  /// Information about an interval volume vertex.
  class DUAL_IVOLVERT: public IJKDUAL::DUAL_ISOVERT {

  public:

    typedef IVOLDUAL_TABLE_VERTEX_INFO::DIR_BITS_TYPE DIR_BITS_TYPE;

    /// Bit flag for connection direction.
    /// If bit i is 1, then vertex connects across facet i.
    DIR_BITS_TYPE connect_dir;

    /// Number of hexahedra incident on the vertex in the lookup table.
    /// - If there are no missing ivol hexahedra, then num_incident_hex
    ///   equals the number of incident hexahedra in the mesh.
    DEGREE_TYPE num_incident_hex;

    /// Number of isosurface quadrilaterals incident on the vertex 
    ///   in the lookup table.
    /// - If there are no missing ivol hexahedra, then num_incident_iso_quad
    ///   equals the number of incident isosurface quadrilaterals in the mesh.
    DEGREE_TYPE num_incident_iso_quad;

    // *** CHANGE TO flag_on_lower_isosurface ***
    /// True if vertex is on the lower isosurface.
    bool flag_lower_isosurface;

    // *** CHANGE TO flag_on_upper_isosurface ***
    /// True if vertex is on the upper isosurface.
    /// - A vertex cannot be in both the lower and upper isosurface.
    /// - A vertex could be in neither lower nor upper isosurface.
    bool flag_upper_isosurface;

    /// If true, the mesh is missing some interval volume hexahedra
    ///   which would normally be incident on the vertex.
    /// The hexahedra are missing because they are dual to grid edges
    ///   or vertices which are on the grid boundary.
    bool flag_missing_ivol_hexahedra;

    /// Separation vertex.
    VERTEX_INDEX separation_vertex;

    /// Separation edge direction.
    /// Separation edge has lower/left endpoint separation_vertex.
    VERTEX_INDEX separation_edge_direction;

    /// If ivol vertex is doubly connected, then ivol vertex is
    ///   doubly connected across facet doubly_connected_facet.
    FACET_INDEX doubly_connected_facet;

    bool is_doubly_connected;        ///< True, if vertex is doubly connected.
    bool in_loop;                    ///< True, if vertex is in loop.
    bool in_box;                     ///< True, if vertex is in box.
    bool in_pseudobox;               ///< True, if vertex is in pseudobox.
    bool in_thin_region;             ///< True, if vertex is in thin region.

    /// If true, adjacent to upper isosurface.
    bool flag_adjacent_to_upper_isosurface;

    /// If true, adjacent to lower isosurface.
    bool flag_adjacent_to_lower_isosurface;

    /// Map for degenerate mesh.
    ISO_VERTEX_INDEX map_to;

  protected:
    void Init();                    ///< Initialization function.

  public:
    DUAL_IVOLVERT() { Init(); }

    /// Return number of incident hexahedra in the lookup table.
    DEGREE_TYPE NumIncidentHex() const
    { return(num_incident_hex); }

    /// Return number of incident isosurface quadrilaterals in the lookup table.
    DEGREE_TYPE NumIncidentIsoQuad() const
    { return(num_incident_iso_quad); }

    /// Return separation vertex.
    VERTEX_INDEX SeparationVertex() const
    { return(separation_vertex); }
  };
  typedef std::vector<DUAL_IVOLVERT> DUAL_IVOLVERT_ARRAY;
  

  // ***************************************************************
  // ACTIVE GRID CUBE
  // ***************************************************************  

  typedef IJK::ACTIVE_GRID_CUBE_WITH_COORD_T
  <3, ISO_VERTEX_INDEX, TABLE_INDEX, FACET_VERTEX_INDEX, ISO_VERTEX_INDEX,
   GRID_COORD_TYPE> ACTIVE_GRID_CUBE;


  // **************************************************
  // INTERVAL VOLUME POLY INFO
  // **************************************************

  class IVOLDUAL_POLY_INFO {

  protected:
    void Init();

  public:

    /// Constructor.
    IVOLDUAL_POLY_INFO() { Init(); };

    /// If true, dual to grid edge.
    /// Else dual to grid vertex.
    bool flag_dual_to_edge;

    /// If true, reverse the orientation of the polytope.
    bool flag_reverse_orient;

    /// If true, the polytope is subdivided into children polytopes.
    bool flag_subdivide_hex;

    /// Grid vertex or lower/leftmost endoint of grid edge.
    ISO_VERTEX_INDEX v0;

    /// Direction of dual edge if dual to grid edge.
    int edge_direction;

    void SetDualToEdge
    (const ISO_VERTEX_INDEX iv0, const int edge_direction)
    {
      flag_dual_to_edge = true;
      v0 = iv0;
      this->edge_direction = edge_direction;
      flag_reverse_orient = false;
    }

    void SetDualToEdge
    (const ISO_VERTEX_INDEX iv0, const int edge_direction,
     const bool flag_reverse_orient)
    {
      flag_dual_to_edge = true;
      v0 = iv0;
      this->edge_direction = edge_direction;
      this->flag_reverse_orient = flag_reverse_orient;
    }

    void SetDualToVertex(const ISO_VERTEX_INDEX iv0)
    {
      flag_dual_to_edge = false;
      v0 = iv0;
    }

  };


  /// Array of IVOLDUAL_POLY_INFO.
  typedef std::vector<IVOLDUAL_POLY_INFO> IVOLDUAL_POLY_INFO_ARRAY;


  // **************************************************
  // DUAL CONTOURING INTERVAL VOLUME FLAGS
  // **************************************************

  /// Data structure for flags controlling interval volume construction.<br>
  /// To add a flag:
  /// - Add a field in this class;
  /// - Initialize the field in IVOLDUAL_DATA_FLAGS::Init()
  ///   in file ivoldual_datastruct.cxx;
  /// - Add a command line option type to OPTION_TYPE in file ivoldualIO.cxx;
  /// - Add a command line option using function option.AddOptionNoArg(),
  ///   or option.AddOption1Arg, in procedure set_command_line_options 
  ///   in file ivoldualIO.cxx.  Usage message and the help message are
  ///   passed as parameters in those commands;
  /// - Add a case to the switch statement in procedure process_option
  ///   that sets the appropriate IVOLDUAL_DATA_FLAGS field when the
  ///   option is invoked on the command line.
  class IVOLDUAL_DATA_FLAGS:public IJKDUAL::DUALISO_DATA_FLAGS {

  protected:
    void Init();

  public:

    /// If true, split vertices in ambiguous cubes sharing facets.
    bool flag_split_ambig_pairs;
    bool flag_split_ambig_pairsB;
    bool flag_split_ambig_pairsC;
    bool flag_split_ambig_pairsD;

    /// Hexahedral triangulation method.
    HEX_TRI_METHOD hex_tri_method;

    /// If true, add extra vertices dual to iso hex for triangulation.
    bool flag_add_isov_dual_to_hexahedra;

    /// Hexahedra orientation.
    /// - If true, orient hexahedra so that facet normals point in.
    bool flag_orient_in;

    /// If true, change vertices values to remove ambiguous facets.
    bool flag_rm_non_manifold;

    /// Flags for smoothing.
    bool flag_lsmooth_elength;
    bool flag_lsmooth_jacobian;
    bool flag_gsmooth_jacobian;

    int lsmooth_elength_iter;
    int lsmooth_jacobian_iter;
    int gsmooth_jacobian_iter;

    float elength_threshold;
    float jacobian_threshold;

    bool flag_split_hex;
    bool flag_collapse_hex;

    float split_hex_threshold;
    float collapse_hex_threshold;

    /// If true, expand thin regions.
    ///   Move isosurface vertices in thin regions away from cube facets.
    bool flag_expand_thin_regions;

    /// Separation distance for thin regions.
    COORD_TYPE thin_separation_distance;

    /// Default encoded value (1 or 2) for vertices in the interior
    ///   of the interval volume.
    GRID_VERTEX_ENCODING default_interior_code;

    /// Set interior code based on scalar values.
    /// - Vertices with scalar value between isovalue0 
    ///     and (isovalue0+isovalue1)/2 receive encoded value 1.
    /// - Vertices with scalar value between (isovalue0+isovalue1)/2 
    ///     and isovalue1 receive encoded value 2.
    /// - Used for case analysis.
    bool flag_set_interior_code_from_scalar;

  public:

    /// Constructor.
    IVOLDUAL_DATA_FLAGS() { Init(); };

    /// Set flags.
    void Set(const IVOLDUAL_DATA_FLAGS & data_flags)
    { *this = data_flags; };
  };


  // **************************************************
  // DUAL CONTOURING INTERVAL VOLUME
  // **************************************************

  /// Representation of dDual contouring interval volume created
  /// by dual contouring interval volume algorithm.
  class DUAL_INTERVAL_VOLUME:
    public IJKDUAL::DUAL_ISOSURFACE_BASE<IVOLDUAL_POLY_INFO> {

  public:
    /// List of interval volume vertices.
    DUAL_IVOLVERT_ARRAY ivolv_list;

  public:
    DUAL_INTERVAL_VOLUME
    (const int dimension, const VERTEX_INDEX numv_per_ivolpoly):
      DUAL_ISOSURFACE_BASE<IVOLDUAL_POLY_INFO>(dimension, numv_per_ivolpoly)
    {};
  };


  // **************************************************
  // IVOLDUAL INPUT DATA AND DATA STRUCTURES
  // **************************************************

  /// Input data to ivoldual.
  class IVOLDUAL_DATA:
    public IJKDUAL::DUALISO_DATA_BASE<IVOLDUAL_SCALAR_GRID,IVOLDUAL_DATA_FLAGS> 
  {
  public:
    IVOLDUAL_DATA() {}; 

    /// Copy, subsample, supersample or subdivide scalar grid.
    /// @pre At most one of flag_subsample, flag_supersample or
    ///   flag_subdivide may be true.
    void SetScalarGrid
      (const DUALISO_SCALAR_GRID_BASE & scalar_grid2, 
       const bool flag_subsample, const int subsample_resolution,
       const bool flag_supersample, const int supersample_resolution, 
       const bool flag_subdivide, const bool flag_rm_diag_ambig,
       const bool flag_add_outer_layer, 
       const GRID_VERTEX_ENCODING interior_code,
       const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    void SubdivideScalarGrid      /// Subdivide scalar_grid.
      (const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    void EvaluateCubeCenter      /// Subdivide scalar_grid.
      (const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    bool EvaluateSubdivideCenter
      (int corner[], int edge[], int icenter,
       const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    bool CheckManifold
      (int corner[], int edge[], int icenter,
       const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1);

    void EliminateAmbigFacets
      (const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1, 
       VERTEX_INDEX & changes_of_ambiguity);

    int RmDiagonalAmbig
      (const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1, 
       int caseID);

    int symbol
      (const int cur, const SCALAR_TYPE v0, const SCALAR_TYPE v1);

  };

  // **************************************************
  // DUALISO TIME
  // **************************************************

  typedef IJKDUAL::DUALISO_TIME DUALISO_TIME;


  // **************************************************
  // DUALISO INFO
  // **************************************************

  class IVOLDUAL_INFO: public IJKDUAL::DUALISO_INFO {

  public:
    IVOLDUAL_INFO();
    IVOLDUAL_INFO(const int dimension);

    int num_non_manifold_changes;
    int num_vertices_moved_in_expand_thin;

    void Clear(); // clear all data
  };

  // **************************************************
  // HEXAHEDRAL MESH
  // **************************************************

  typedef IJK::POLYMESH<VERTEX_INDEX, int> IVOL_POLYMESH;

  typedef IJKDUAL::VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST
  <int,ISO_VERTEX_INDEX,FACET_INDEX,int>
  IVOL_VERTEX_ADJACENCY_LIST;

  typedef IJK::VERTEX_POLY_INCIDENCE_WITH_VLOC<int, int, int>
  IVOL_VERTEX_POLY_INCIDENCE;

  // **************************************************
  // MERGE DATA
  // **************************************************

  typedef IJKDUAL::MERGE_DATA MERGE_DATA;


  // **************************************************
  // INTERVAL VOLUME MESH VERT INFO
  // **************************************************

  class MESH_VERT_INFO {

  protected:
    void Init();

  public:

    /// Constructor.
    MESH_VERT_INFO() { Init(); };

    /// If false, the vertex is created during subdivision.
    /// Else the vertex is dual to cube
    bool is_dual;

    /// Dual cube index
    VERTEX_INDEX cube_index;

  };

  /// Array of IVOLDUAL_POLY_INFO.
  typedef std::vector<MESH_VERT_INFO> MESH_VERT_INFO_ARRAY;


  // **************************************************
  // INTERVAL VOLUME MESH POLY INFO
  // **************************************************

  class MESH_POLY_INFO {

  protected:
    void Init();

  public:

    /// Constructor.
    MESH_POLY_INFO() { Init(); };

    /// If true, the polytope is subdivided.
    bool is_subdivided;

    /// If true, the polytope is a dual polytope.
    bool is_dual;

  };

  /// Array of IVOLDUAL_POLY_INFO.
  typedef std::vector<MESH_POLY_INFO> MESH_POLY_INFO_ARRAY;

}

#endif
