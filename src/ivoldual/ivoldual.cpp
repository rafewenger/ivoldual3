/// @file ivoldual.cpp
/// @brief Generate interval volume in arbitrary dimensions using dual contouring.
/// - Version 0.3.1

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


#include "ijkisopoly.tpp"
#include "ijkinterpolate.tpp"
#include "ijkmesh.tpp"
#include "ijkmesh_datastruct.tpp"
#include "ijktime.tpp"

#include "ijkdual_extract.tpp"

#include "ivoldual_ivolpoly.tpp"

#include "ivoldual.h"
#include "ivoldual_datastruct.h"
#include "ivoldual_compute.h"
#include "ivoldual_query.h"
#include "ivoldual_move.h"
#include "ivoldual_reposition.h"
#include "ivoldual_divide_hex.h"

#include "ijktriangulate.txx"

using namespace IJK;
using namespace IVOLDUAL;

// *** DEBUG ***
#include "ijkprint.tpp"


// **************************************************
// DUAL CONTOURING INTERVAL VOLUME
// **************************************************

// Construct interval volume using dual contouring.
void IVOLDUAL::dual_contouring_interval_volume
(const IVOLDUAL_DATA & ivoldual_data, 
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 DUAL_INTERVAL_VOLUME & dual_interval_volume, IVOLDUAL_INFO & dualiso_info)
{
  const int dimension = ivoldual_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = ivoldual_data.ScalarGrid().AxisSize();
  PROCEDURE_ERROR error("dual_contouring_interval_volume");

  clock_t t_start = clock();

  if (!ivoldual_data.Check(error)) { throw error; };

  dual_interval_volume.Clear();
  dualiso_info.time.Clear();

  IJKDUAL::ISO_MERGE_DATA merge_data(dimension, axis_size);

  dual_contouring_interval_volume
    (ivoldual_data.ScalarGrid(), isovalue0, isovalue1, ivoldual_data,
     dual_interval_volume.isopoly_vert, dual_interval_volume.isopoly_info, 
     dual_interval_volume.ivolv_list,
     dual_interval_volume.vertex_coord, merge_data, dualiso_info);

  // store times
  clock_t t_end = clock();
  clock2seconds(t_end-t_start, dualiso_info.time.total);
}


// Construct interval volume using dual contouring.
// - Returns list of interval volume polytope vertices
//   and list of interval volume vertex coordinates.
// - Version which creates ivoldual_table and ivolv_list.
void IVOLDUAL::dual_contouring_interval_volume
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 const IVOLDUAL_DATA_FLAGS & param,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 IVOLDUAL_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const bool flag_separate_neg = param.SeparateNegFlag();

  IVOLDUAL_CUBE_TABLE
    ivoldual_table(dimension, flag_separate_neg);
  DUAL_IVOLVERT_ARRAY ivolv_list;

  dual_contouring_interval_volume
    (scalar_grid, isovalue0, isovalue1, ivoldual_table, param, ivolpoly_vert, 
     ivolpoly_info, ivolv_list, vertex_coord, merge_data, dualiso_info);
}


// Construct interval volume using dual contouring.
// - Returns list of interval volume polytope vertices
//   and list of interval volume vertex coordinates.
// - Version which creates ivoldual_table.
void IVOLDUAL::dual_contouring_interval_volume
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 const IVOLDUAL_DATA_FLAGS & param,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 IVOLDUAL_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const bool flag_separate_neg = param.SeparateNegFlag();

  IVOLDUAL_CUBE_TABLE
    ivoldual_table(dimension, flag_separate_neg);

  dual_contouring_interval_volume
    (scalar_grid, isovalue0, isovalue1, ivoldual_table, param, ivolpoly_vert, 
     ivolpoly_info, ivolv_list, vertex_coord, merge_data, dualiso_info);
}


// Construct interval volume using dual contouring.
// - Returns list of interval volume polytope vertices
//   and list of interval volume vertex coordinates.
// - Version which creates cube_ivolv_list.
void IVOLDUAL::dual_contouring_interval_volume
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const IVOLDUAL_DATA_FLAGS & param,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord,
 MERGE_DATA & merge_data, 
 IVOLDUAL_INFO & dualiso_info)
{
  std::vector<GRID_CUBE_DATA> cube_ivolv_list;

  dual_contouring_interval_volume
    (scalar_grid, isovalue0, isovalue1, ivoldual_table, param, ivolpoly_vert,
     cube_ivolv_list, ivolv_list, ivolpoly_info, vertex_coord, 
     merge_data, dualiso_info);
}


// Construct interval volume using dual contouring.
// Returns list of isosurface polytope vertices
//   and list of isosurface vertex coordinates
void IVOLDUAL::dual_contouring_interval_volume
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
 IVOLDUAL_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const bool flag_separate_neg = param.SeparateNegFlag();
  const bool flag_always_separate_opposite(true);
  const GRID_VERTEX_ENCODING default_interior_code =
    param.default_interior_code;
  const VERTEX_INDEX num_grid_vertices = scalar_grid.NumVertices();
  int num_non_manifold_split(0);
  IJK::ARRAY<VERTEX_INDEX> index_to_cube_list
    (num_grid_vertices, num_grid_vertices);
  IVOL_POLYMESH polymesh;
  IVOL_VERTEX_ADJACENCY_LIST vertex_adjacency_list(dimension);
  CUBE_FACE_INFO cube_info(dimension);
  IJK::PROCEDURE_ERROR error("dual_contouring_interval_volume");
  clock_t t0, t1, t2;

  if (scalar_grid.Dimension() != ivoldual_table.Dimension()) {
    error.AddMessage
      ("Programming error.  Incorrect isodual table dimension.");
    error.AddMessage
      ("  Interval volume table dimension does not match scalar grid dimension.");
    error.AddMessage
      ("    Interval volume table dimension: ", ivoldual_table.Dimension(), "");
    error.AddMessage
      ("    Scalar grid dimension: ", scalar_grid.Dimension(), "");
    throw error;
  }

  t0 = clock();

  ivolpoly_vert.clear();
  dualiso_info.time.Clear();

  IVOLDUAL_ENCODED_GRID encoded_grid;
  if (param.flag_set_interior_code_from_scalar) {
    encode_grid_vertices_set_interior_from_scalar
      (scalar_grid, isovalue0, isovalue1, encoded_grid, dualiso_info);
  }
  else {
    encode_grid_vertices
      (scalar_grid, isovalue0, isovalue1, default_interior_code, 
       encoded_grid, dualiso_info);
  }

  std::vector<ISO_VERTEX_INDEX> ivolpoly;
  std::vector<POLY_VERTEX_INDEX> poly_vertex;
  extract_dual_ivolpoly
    (encoded_grid, ivolpoly, poly_vertex, ivolpoly_info, dualiso_info);
  t1 = clock();

  std::vector<ISO_VERTEX_INDEX> cube_list;
  std::vector<ISO_VERTEX_INDEX> ivolpoly_cube;
  merge_identical(ivolpoly, cube_list, ivolpoly_cube, merge_data);
  t2 = clock();

  set_grid_cube_indices(cube_list, cube_ivolv_list);
  set_grid_coord(scalar_grid, cube_ivolv_list);

  set_cube_ivoltable_info
    (encoded_grid, ivoldual_table, cube_ivolv_list);
  IJK::set_index_to_cube_list(cube_ivolv_list, index_to_cube_list);

  if (param.flag_split_ambig_pairsB) {
    // *** Probably not necessary
    split_non_manifold_ivolv_pairs_ambigB
      (encoded_grid, ivoldual_table, cube_ivolv_list, num_non_manifold_split);
  }
  else if (param.flag_split_ambig_pairsC) {
    split_non_manifold_ivolv_pairs_ambigC
      (encoded_grid, ivoldual_table, cube_ivolv_list, num_non_manifold_split);
  }
  else if (param.flag_split_ambig_pairs) {
    split_non_manifold_ivolv_pairs_ambig
      (encoded_grid, ivoldual_table, cube_ivolv_list, num_non_manifold_split);
  }
  else if (param.flag_split_ambig_pairsD) {
    split_non_manifold_ivolv_pairs_ambigD
      (encoded_grid, ivoldual_table, cube_ivolv_list, num_non_manifold_split);
  }

  VERTEX_INDEX num_split;
  split_dual_ivolvert
    (ivoldual_table, ivolpoly_cube, poly_vertex, ivolpoly_info, 
     cube_ivolv_list, ivolv_list, ivolpoly_vert, num_split);

  IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG isodual_table
    (dimension, flag_separate_neg, flag_always_separate_opposite);
  position_all_dual_ivol_vertices
    (scalar_grid, ivoldual_table, isodual_table, isovalue0, isovalue1, 
     ivolv_list, vertex_coord);

  polymesh.AddPolytopes(ivolpoly_vert, cube_info.NumVertices());
  vertex_adjacency_list.SetFromMeshOfCubes(polymesh, cube_info);
  vertex_adjacency_list.SetAllDualFacetsFromVertexAndCubeLists
    (ivolv_list, cube_ivolv_list);

  set_ivol_vertex_info
    (scalar_grid, ivoldual_table, ivolpoly_vert, 
     cube_ivolv_list, index_to_cube_list.PtrConst(),
     vertex_adjacency_list, ivolv_list);

  // Expand thin regions.
  if (param.flag_expand_thin_regions) {
    int num_moved;
    expand_thin_regions
      (scalar_grid, cube_ivolv_list, ivolv_list,
       param.thin_separation_distance, vertex_coord, num_moved);
  }

  // Polytopes dual to vertex.
  const int NUM_VERT_PER_HEXAHEDRON(8);
  IJK::POLYMESH_DATA<VERTEX_INDEX,int, 
    IJK::HEX_TRIANGULATION_INFO<char,char>> hex_data;
  hex_data.AddPolytopes(ivolpoly_vert, NUM_VERT_PER_HEXAHEDRON);
  IJK::VERTEX_POLY_INCIDENCE<int,int> vertex_poly_incidence(hex_data);

  // Split or Collapse hexahedron to improve Jacobian.
  if (param.flag_split_hex) {
    split_hex
    (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, ivolv_list,
     ivolpoly_info, vertex_coord, param.split_hex_threshold);
  }
  if (param.flag_collapse_hex) {
    collapse_hex
    (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, ivolv_list, 
     ivolpoly_info, vertex_coord, param.collapse_hex_threshold);
  }
  // Edge length improvement.
  if (param.flag_lsmooth_elength) {
    laplacian_smooth_elength
    (ivoldual_table, vertex_adjacency_list, ivolv_list, 
     vertex_coord, param.elength_threshold, param.lsmooth_elength_iter);         
  }

  // Jacobian improvement.
  if (param.flag_lsmooth_jacobian) {
    laplacian_smooth_jacobian
    (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, vertex_poly_incidence, ivolv_list, 
     vertex_coord, param.jacobian_threshold, param.lsmooth_jacobian_iter);
  } 
  else if (param.flag_gsmooth_jacobian) {
    gradient_smooth_jacobian
    (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, vertex_poly_incidence, ivolv_list, 
     ivolpoly_info, vertex_coord, param.jacobian_threshold, param.gsmooth_jacobian_iter);
  } 


  if (param.flag_orient_in) {
    const int num_vert_per_cube_facet =  
      compute_num_cube_facet_vertices(dimension);
    IJK::reverse_orientations_cube_list
      (ivolpoly_vert, num_vert_per_cube_facet);
  }

  dualiso_info.scalar.num_non_empty_cubes = cube_list.size();
  dualiso_info.multi_isov.num_cubes_multi_isov = num_split;
  dualiso_info.multi_isov.num_cubes_single_isov =
    cube_list.size() - num_split;
  dualiso_info.multi_isov.num_non_manifold_split = num_non_manifold_split;

  // store times
  IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  IJK::clock2seconds(t2-t1, dualiso_info.time.merge);

}


// **************************************************
// ENCODE VERTICES
// **************************************************

void IVOLDUAL::encode_grid_vertices
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 const GRID_VERTEX_ENCODING default_interior_code,
 IVOLDUAL_ENCODED_GRID & encoded_grid,
 IVOLDUAL_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  encoded_grid.SetSize(dimension, axis_size);

  for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
    if (scalar_grid.Scalar(iv) < isovalue0) 
      { encoded_grid.Set(iv, 0); }
    else if (scalar_grid.Scalar(iv) > isovalue1) 
      { encoded_grid.Set(iv, 3); }
    else {
      // Set all vertices with scalar value 
      //   in range [isovalue0,isovalue1] to default_interior_code
      encoded_grid.Set(iv,default_interior_code);
    }
  }
}


// Encode grid vertices. Set interior codes based on scalar grid.
// Used for case analysis.
void IVOLDUAL::encode_grid_vertices_set_interior_from_scalar
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0,  const SCALAR_TYPE isovalue1, 
 IVOLDUAL_ENCODED_GRID & encoded_grid,
 IVOLDUAL_INFO & dualiso_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE isovalue_average = (isovalue0+isovalue1)/2.0;

  encoded_grid.SetSize(dimension, axis_size);

  for (VERTEX_INDEX iv = 0; iv < scalar_grid.NumVertices(); iv++) {
    if (scalar_grid.Scalar(iv) < isovalue0) 
      { encoded_grid.Set(iv, 0); }
    else if (scalar_grid.Scalar(iv) > isovalue1) 
      { encoded_grid.Set(iv, 3); }
    else if (scalar_grid.Scalar(iv) < isovalue_average) 
      { encoded_grid.Set(iv, 1); }
    else {
      // (scalar_grid.Scalar(iv) > isovalue_average)
      encoded_grid.Set(iv, 2);
    }
  }
}


// **************************************************
// SET IVOLTABLE INFO FOR EACH ACTIVE GRID CUBE
// **************************************************

namespace {

  TABLE_INDEX compute_table_index_from_encoded_grid
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   const int num_vertex_types,
   const VERTEX_INDEX icube)
  {
    const int num_cube_vertices = encoded_grid.NumCubeVertices();
    TABLE_INDEX table_index = 0;
    TABLE_INDEX factor = 1;
    for (int j = 0; j < num_cube_vertices; j++) {
      const int j2 = num_cube_vertices-1-j;
      const VERTEX_INDEX jv2 = encoded_grid.CubeVertex(icube, j2);
      table_index = 
        (table_index*num_vertex_types) + encoded_grid.Scalar(jv2);
    }

    return(table_index);
  }

}

void IVOLDUAL::set_cube_ivoltable_info
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 std::vector<GRID_CUBE_DATA> & cube_ivolv_list)
{
  const int dimension = encoded_grid.Dimension();
  const int num_vertex_types = ivoldual_table.NumVertexTypes();

  for (int i = 0; i < cube_ivolv_list.size(); i++) {
    const VERTEX_INDEX cube_index = cube_ivolv_list[i].cube_index;
    const TABLE_INDEX table_index =
      compute_table_index_from_encoded_grid
      (encoded_grid, num_vertex_types, cube_index);

    cube_ivolv_list[i].table_index = table_index;
    cube_ivolv_list[i].num_isov = ivoldual_table.NumIsoVertices(table_index);
  }
}


// **************************************************
// SET IVOL VERTEX INFORMATION
// **************************************************

void IVOLDUAL::set_ivol_vertex_info
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const std::vector<ISO_VERTEX_INDEX> & poly_vert,
 const std::vector<GRID_CUBE_DATA> & cube_list,
 const VERTEX_INDEX index_to_cube_list[],
 const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  typedef IVOLDUAL_TABLE_VERTEX_INFO::CUBE_VERTEX_TYPE CUBE_VERTEX_TYPE;
  typedef IVOLDUAL_TABLE_VERTEX_INFO::CUBE_EDGE_TYPE CUBE_EDGE_TYPE;

  const int DIM3(3);
  const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
    IVOLDUAL_TABLE_VERTEX_INFO::UNDEFINED_CUBE_VERTEX;
  const CUBE_EDGE_TYPE UNDEFINED_CUBE_EDGE = 
    IVOLDUAL_TABLE_VERTEX_INFO::UNDEFINED_CUBE_EDGE;
  const CUBE_FACE_INFO cube(DIM3);


  for (IVOL_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) {
    const TABLE_INDEX table_index = ivolv_list[ivolv].table_index;
    const FACET_VERTEX_INDEX patch_index = ivolv_list[ivolv].patch_index;
    const VERTEX_INDEX cube_index = ivolv_list[ivolv].cube_index;

    // Initialize.
    ivolv_list[ivolv].separation_vertex = grid.NumVertices();
    ivolv_list[ivolv].separation_edge_direction = 0;

    ivolv_list[ivolv].connect_dir = 
      ivoldual_table.VertexInfo(table_index, patch_index).connect_dir;

    ivolv_list[ivolv].num_incident_hex = 
      ivoldual_table.VertexInfo(table_index, patch_index).NumIncidentPoly();

    ivolv_list[ivolv].num_incident_iso_quad = 
      ivoldual_table.VertexInfo(table_index, patch_index).NumIncidentIsoPoly();

    ivolv_list[ivolv].flag_lower_isosurface = 
      ivoldual_table.OnLowerIsosurface(table_index, patch_index);

    ivolv_list[ivolv].flag_upper_isosurface = 
      ivoldual_table.OnUpperIsosurface(table_index, patch_index);

    const CUBE_VERTEX_TYPE icorner =
      ivoldual_table.VertexInfo(table_index, patch_index).separation_vertex;

    if (icorner != UNDEFINED_CUBE_VERTEX) {
      ivolv_list[ivolv].separation_vertex = 
        grid.CubeVertex(cube_index, icorner);
    }

    const CUBE_EDGE_TYPE iedge =
      ivoldual_table.VertexInfo(table_index, patch_index).separation_edge;

    if (iedge != UNDEFINED_CUBE_EDGE) {
      const CUBE_VERTEX_TYPE icorner = cube.EdgeEndpoint(iedge, 0);
      ivolv_list[ivolv].separation_vertex = 
        grid.CubeVertex(cube_index, icorner);
      ivolv_list[ivolv].separation_edge_direction = cube.EdgeDir(iedge);
    }

    ivolv_list[ivolv].is_doubly_connected =
      ivoldual_table.VertexInfo(table_index, patch_index).is_doubly_connected;
    ivolv_list[ivolv].doubly_connected_facet =
      ivoldual_table.VertexInfo(table_index, patch_index).
      doubly_connected_facet;
    ivolv_list[ivolv].map_to = ivolv;
  }



  determine_ivol_vertices_missing_incident_hex
    (ivoldual_table, poly_vert, ivolv_list);

  determine_ivol_vertices_in_isosurface_boxes
    (vertex_adjacency_list, ivolv_list);

  determine_ivol_vertices_in_isosurface_pseudoboxes
    (vertex_adjacency_list, cube_list, ivolv_list);

  determine_ivol_vertices_adjacent_to_upper_and_lower_isosurfaces
    (vertex_adjacency_list, ivolv_list);

  determine_thin_regions
    (grid, cube_list, index_to_cube_list, vertex_adjacency_list, ivolv_list);

}


void IVOLDUAL::determine_ivol_vertices_missing_incident_hex
(const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const std::vector<ISO_VERTEX_INDEX> & poly_vert,
 DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  std::vector<ISO_VERTEX_INDEX> 
    num_hex_incident_on_ivolv(ivolv_list.size(), 0);

  for (int i = 0; i < poly_vert.size(); i++) {
    ISO_VERTEX_INDEX ivolv = poly_vert[i];
    num_hex_incident_on_ivolv[ivolv]++;
  }

  for (ISO_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) {
    const TABLE_INDEX table_index = ivolv_list[ivolv].table_index;
    const FACET_VERTEX_INDEX patch_index = ivolv_list[ivolv].patch_index;
    const int num_incident_hex = 
      ivoldual_table.VertexInfo(table_index, patch_index).NumIncidentPoly();

    if (num_incident_hex == num_hex_incident_on_ivolv[ivolv]) {
      ivolv_list[ivolv].flag_missing_ivol_hexahedra = false;
    }
    else {
      ivolv_list[ivolv].flag_missing_ivol_hexahedra = true;
    }
  }

}


void IVOLDUAL::determine_ivol_vertices_in_isosurface_boxes
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  const int NUM_CUBE_VERTICES(8);
  VERTEX_INDEX box_corner[NUM_CUBE_VERTICES];

  for (IVOL_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) 
    { ivolv_list[ivolv].in_box = false; }

  for (IVOL_VERTEX_INDEX ivolv = 0; 
       ivolv < vertex_adjacency_list.NumVertices(); ivolv++) {

    if (is_ivolv_lower_left_vertex_of_isosurface_box
        (vertex_adjacency_list, ivolv_list, ivolv, box_corner)) {

      for (int j = 0; j < NUM_CUBE_VERTICES; j++) 
        { ivolv_list[box_corner[j]].in_box = true; }

    }
  }
}


void IVOLDUAL::determine_ivol_vertices_in_isosurface_pseudoboxes
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const std::vector<GRID_CUBE_DATA> & cube_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  const int NUM_CUBE_VERTICES(8);
  VERTEX_INDEX cube_list_index[NUM_CUBE_VERTICES];

  for (IVOL_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) 
    { ivolv_list[ivolv].in_pseudobox = false; }

  for (int i = 0; i < cube_list.size(); i++) {

    if (is_upper_right_grid_vertex_in_isosurface_pseudobox
        (vertex_adjacency_list, cube_list, ivolv_list, i, cube_list_index)) {

      for (int j = 0; j < NUM_CUBE_VERTICES; j++) 
        { set_in_pseudobox(cube_list[cube_list_index[j]], ivolv_list); }
    }
  }
}


void IVOLDUAL::set_in_pseudobox
(const GRID_CUBE_DATA & cube_ivolv,
 DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  for (int i = 0; i < cube_ivolv.num_isov; i++) {
    const IVOL_VERTEX_INDEX ivolv = cube_ivolv.first_isov + i;
    if (ivolv_list[ivolv].NumIncidentIsoQuad() > 0) 
      { ivolv_list[ivolv].in_pseudobox = true; }
  }
}


void IVOLDUAL::determine_ivol_vertices_adjacent_to_upper_and_lower_isosurfaces
(const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  const int NUM_CUBE_VERTICES(8);
  VERTEX_INDEX box_corner[NUM_CUBE_VERTICES];

  for (IVOL_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) {
    ivolv_list[ivolv].flag_adjacent_to_upper_isosurface = false;
    ivolv_list[ivolv].flag_adjacent_to_lower_isosurface = false;
  }


  for (IVOL_VERTEX_INDEX ivolv0 = 0; 
       ivolv0 < vertex_adjacency_list.NumVertices(); ivolv0++) {

    ivolv_list[ivolv0].flag_adjacent_to_lower_isosurface = false;
    ivolv_list[ivolv0].flag_adjacent_to_upper_isosurface = false;

    for (int j = 0; j < vertex_adjacency_list.NumAdjacent(ivolv0); j++) {
      IVOL_VERTEX_INDEX ivolv1 = 
        vertex_adjacency_list.AdjacentVertex(ivolv0, j);

      if (ivolv_list[ivolv1].flag_lower_isosurface) 
        { ivolv_list[ivolv0].flag_adjacent_to_lower_isosurface = true; }

      if (ivolv_list[ivolv1].flag_upper_isosurface) 
        { ivolv_list[ivolv0].flag_adjacent_to_upper_isosurface = true; }
    }
  }
}


void IVOLDUAL::determine_thin_regions
(const DUALISO_GRID & grid,
 const std::vector<GRID_CUBE_DATA> & cube_list,
 const VERTEX_INDEX index_to_cube_list[],
 const IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list)
{
  typedef DUAL_IVOLVERT::DIR_BITS_TYPE DIR_BITS_TYPE;

  const int DIM3(3);
  DIR_BITS_TYPE mask[DIM3] = 
    { DIR_BITS_TYPE(1), DIR_BITS_TYPE(2), DIR_BITS_TYPE(4) };
  IJK::PROCEDURE_ERROR error("determine_thin_regions");

  for (IVOL_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) 
    { ivolv_list[ivolv].in_thin_region = false; }

  for (IVOL_VERTEX_INDEX ivolv0 = 0; 
       ivolv0 < vertex_adjacency_list.NumVertices(); ivolv0++) {

    const VERTEX_INDEX cube_index0 = ivolv_list[ivolv0].cube_index;

    // Skip vertices which are not on isosurface or not bent degree 4.
    if (ivolv_list[ivolv0].NumIncidentIsoQuad() != 4) { continue; }
    if (ivolv_list[ivolv0].separation_vertex >= grid.NumVertices()) 
      { continue; }

    if (ivolv_list[ivolv0].in_thin_region) {
      // Already flagged as in thin region
      continue;
    }

    if (ivolv_list[ivolv0].flag_missing_ivol_hexahedra) { 
      // Separation edge is on grid boundary.
      continue; 
    }

    const int edge_dir = ivolv_list[ivolv0].separation_edge_direction;

    const int d1 = (edge_dir+1)%DIM3;
    const int d2 = (edge_dir+2)%DIM3;
    const DIR_BITS_TYPE connect_dir = ivolv_list[ivolv0].connect_dir;

    VERTEX_INDEX cube_index1;
    if ((mask[d1] & connect_dir) == 0) 
      { cube_index1 = grid.NextVertex(cube_index0, d1); }
    else
      { cube_index1 = grid.PrevVertex(cube_index0, d1); }

    if ((mask[d2] & connect_dir) == 0) 
      { cube_index1 = grid.NextVertex(cube_index1, d2); }
    else
      { cube_index1 = grid.PrevVertex(cube_index1, d2); }

    const int i1 = index_to_cube_list[cube_index1];

    if (i1 == grid.NumVertices()) {
      // Grid cube cube_index1 is not active.  Skip.
      continue;
    }

    IVOL_VERTEX_INDEX ivolv1 = cube_list[i1].first_isov;
    for (int j = 0; j < cube_list[i1].num_isov; j++) {

      // Skip vertices which are not on isosurface or not bent degree 4.
      if (ivolv_list[ivolv1].NumIncidentIsoQuad() != 4) { continue; }
      if (ivolv_list[ivolv1].separation_vertex >= grid.NumVertices()) 
      { continue; }

      if ((ivolv_list[ivolv1].separation_vertex ==
           ivolv_list[ivolv0].separation_vertex) &&
          (ivolv_list[ivolv1].separation_edge_direction ==
           ivolv_list[ivolv0].separation_edge_direction)) {
        ivolv_list[ivolv0].in_thin_region = true;
        ivolv_list[ivolv1].in_thin_region = true;
        break;
      }
    }
  }

}


// **************************************************
// EXTRACT DUAL INTERVAL VOLUME POLYTOPES
// **************************************************

void IVOLDUAL::extract_dual_ivolpoly
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly,
 std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 IVOLDUAL_INFO & dualiso_info)
{
  dualiso_info.time.extract = 0;

  clock_t t0 = clock();

  // Initialize output
  ivolpoly.clear();

  if (encoded_grid.NumCubeVertices() < 1) { return; }

  extract_ivolpoly_dual_to_grid_edges
    (encoded_grid, ivolpoly, poly_vertex, ivolpoly_info, dualiso_info);
  extract_ivolpoly_dual_to_grid_vertices
    (encoded_grid, ivolpoly, poly_vertex, ivolpoly_info, dualiso_info);

  clock_t t1 = clock();
  IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
}

namespace {

  void add_grid_facet_vertices
  (const DUALISO_GRID & grid, const VERTEX_INDEX iv0, const int edge_dir,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly,
   std::vector<POLY_VERTEX_INDEX> & poly_vertex)
  {
    const int num_facet_vertices = grid.NumFacetVertices();

    for (int k = 0; k < num_facet_vertices; k++) {

      // jcube is the index of the cube containing the k'th vertex
      //   of the facet.
      const VERTEX_INDEX jcube = grid.FacetVertex(iv0, edge_dir, k);
      ivolpoly.push_back(jcube);
      poly_vertex.push_back(edge_dir*num_facet_vertices+k);
    }
  }

  // \param[out] flag_reverse_orient True if orientation of dual polytope
  //   should be reversed.
  bool does_grid_edge_have_dual_ivolpoly
  (const IVOLDUAL_ENCODED_GRID & encoded_grid,
   const VERTEX_INDEX iend0, const int edge_dir,
   bool & flag_reverse_orient)
  {
    VERTEX_INDEX iend1 = encoded_grid.NextVertex(iend0, edge_dir);
    GRID_VERTEX_ENCODING s0 = encoded_grid.Scalar(iend0);
    GRID_VERTEX_ENCODING s1 = encoded_grid.Scalar(iend1);

    flag_reverse_orient = false;

    if (s0 > s1) { 
      std::swap(s0,s1); 
      flag_reverse_orient = true;
    }

    if (s0 == 0) {
      if (s1 >= 2) { return(true); }
    }
    else if (s0 == 1) {
      if (s1 > 1)
        { return(true); }
    }

    return(false);
  }

  bool does_grid_vertex_have_dual_ivolpoly
  (const IVOLDUAL_ENCODED_GRID & encoded_grid, const VERTEX_INDEX iv0)
  {
    GRID_VERTEX_ENCODING s0 = encoded_grid.Scalar(iv0);

    if (s0 == 1 || s0 == 2) { return(true); }
    return(false);
  }

}

void IVOLDUAL::extract_ivolpoly_dual_to_grid_edges
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly,
 std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 IVOLDUAL_INFO & dualiso_info)
{
  const int num_facet_vertices = encoded_grid.NumFacetVertices();
  bool flag_reverse_orient;

  if (num_facet_vertices == 0) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_EDGE(iend0, edge_dir, encoded_grid, VERTEX_INDEX) {

    if (does_grid_edge_have_dual_ivolpoly
        (encoded_grid, iend0, edge_dir, flag_reverse_orient)) {

      const VERTEX_INDEX increment = 
        encoded_grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);
      const VERTEX_INDEX iv0 = iend0 - increment;

      for (int j = 0; j < 2; j++) {
        add_grid_facet_vertices
          (encoded_grid, iv0, edge_dir, ivolpoly, poly_vertex);
      }

      IVOLDUAL_POLY_INFO info;
      info.SetDualToEdge(iend0, edge_dir, flag_reverse_orient);
      ivolpoly_info.push_back(info);
    }
  }

}


void IVOLDUAL::extract_ivolpoly_dual_to_grid_vertices
(const IVOLDUAL_ENCODED_GRID & encoded_grid,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly,
 std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 IVOLDUAL_INFO & dualiso_info)
{
  const int num_facet_vertices = encoded_grid.NumFacetVertices();
  const int num_cube_vertices = encoded_grid.NumCubeVertices();

  if (num_facet_vertices == 0) { return; }

  IJK_FOR_EACH_INTERIOR_GRID_VERTEX(iv0, encoded_grid, VERTEX_INDEX) {

    if (does_grid_vertex_have_dual_ivolpoly(encoded_grid, iv0)) {
      const VERTEX_INDEX increment = 
        encoded_grid.CubeVertexIncrement(num_cube_vertices-1);
      const VERTEX_INDEX iv1 = iv0 - increment;


      for (int k = 0; k < num_cube_vertices; k++) {
        ivolpoly.push_back(encoded_grid.CubeVertex(iv1, k));
        poly_vertex.push_back(k);
      }

      IVOLDUAL_POLY_INFO info;
      info.SetDualToVertex(iv0);
      ivolpoly_info.push_back(info);
    }
  }
}


// **************************************************
// SPLIT DUAL INTERVAL VOLUME VERTICES
// **************************************************

namespace {

  void set_dual_ivolv_vertices
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const std::vector<GRID_CUBE_DATA> & cube_list,
   const std::vector<VERTEX_INDEX> & ivolpoly_cube, 
   const std::vector<POLY_VERTEX_INDEX> & poly_vertex,
   const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
   std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert)
  {
    const int dimension = ivoldual_table.Dimension();
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    const int num_cube_vertices = cube.NumVertices();
    const int num_facet_vertices = cube.NumFacetVertices();

    ivolpoly_vert.resize(ivolpoly_cube.size());

    for (ISO_VERTEX_INDEX i = 0; i < ivolpoly_cube.size(); i++) {
      const VERTEX_INDEX k = ivolpoly_cube[i];
      const TABLE_INDEX it = cube_list[k].table_index;
      const int ipoly = i/num_cube_vertices;
      const int poly_vertex_i = poly_vertex[i];
      const int icorner = i%num_cube_vertices;

      // Compute index of vertex opposite to poly_vertex[i]
      if (ivolpoly_info[ipoly].flag_dual_to_edge) {
        const bool on_lower_cube_facet = (icorner < num_facet_vertices);
        const int ifacet = cube.FacetIndex(poly_vertex_i);
        const int j = poly_vertex_i - ifacet*num_facet_vertices;
        int opposite_vertex = (num_facet_vertices-1) - j;
        opposite_vertex += (ifacet*num_facet_vertices);

        if (ivolpoly_info[ipoly].flag_reverse_orient ^
            on_lower_cube_facet) {
          // Vertex on lower facet.
          ivolpoly_vert[i] = cube_list[k].first_isov +
            ivoldual_table.LowerIncident(it, opposite_vertex);
        }
        else {
          // Vertex on upper facet.
          ivolpoly_vert[i] = cube_list[k].first_isov +
            ivoldual_table.UpperIncident(it, opposite_vertex);
        }
      }
      else {
        const int opposite_vertex = cube.OppositeVertex(poly_vertex_i);

        ivolpoly_vert[i] = cube_list[k].first_isov +
          ivoldual_table.IncidentIVolVertex(it, opposite_vertex);
      }
    }
  }


  // Return true if ivol vertices in cube are candidates for splitting
  //   to avoid pair ambiguity.
  bool is_cube_ambig_split_candidate
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const GRID_CUBE_DATA & cube_data)
  {
    const TABLE_INDEX table_index = cube_data.table_index;
    const int numv_in_lower_lifted = 
      ivoldual_table.NumVerticesInLowerLifted(table_index);
    const int numv_in_upper_lifted = 
      ivoldual_table.NumVerticesInUpperLifted(table_index);
    const TABLE_INDEX opposite_table_index =
      ivoldual_table.OppositeTableIndex(table_index);

    if (ivoldual_table.NumAmbiguousFacets(table_index) > 1) 
      { return(false); }
    if (numv_in_lower_lifted > 1) { return(false); }
    if (numv_in_upper_lifted > 1) { return(false); }

    if (ivoldual_table.AmbiguousFacetBitsInLowerLifted(table_index) != 0) {
      // Lower lifted cube has an ambiguous facet.
      if (ivoldual_table.NumVerticesInLowerLifted(opposite_table_index) < 2) {
        // Switching to opposite_table_index does not solve
        //   ambiguity problem.
        return(false);
      }
    }

    if (ivoldual_table.AmbiguousFacetBitsInUpperLifted(table_index) != 0) {
      // Lower lifted cube has an ambiguous facet.
      if (ivoldual_table.NumVerticesInUpperLifted(opposite_table_index) < 2) {
        // Switching to opposite_table_index does not solve
        //   ambiguity problem.
        return(false);
      }
    }

    return(true);
  }


  // *** Probably not necessary
  // Return true if ivol vertices in cube are candidates for splitting
  //   to avoid pair ambiguity.
  // Version which allows more splits.
  // Allow split if numv_in_lower_lifted > 1 but numv_in_upper_lifted = 1
  //   or vice versa.
  bool is_cube_ambig_split_candidateB
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const GRID_CUBE_DATA & cube_data)
  {
    const TABLE_INDEX table_index = cube_data.table_index;
    const int numv_in_lower_lifted = 
      ivoldual_table.NumVerticesInLowerLifted(table_index);
    const int numv_in_upper_lifted = 
      ivoldual_table.NumVerticesInUpperLifted(table_index);
    const TABLE_INDEX opposite_table_index =
      ivoldual_table.OppositeTableIndex(table_index);

    if (ivoldual_table.NumAmbiguousFacets(table_index) > 1) 
      { return(false); }
    if (numv_in_lower_lifted > 1 && numv_in_upper_lifted > 1) 
      { return(false); }

    if (ivoldual_table.AmbiguousFacetBitsInLowerLifted(table_index) != 0) {
      // Lower lifted cube has an ambiguous facet.
      if (ivoldual_table.NumVerticesInLowerLifted(opposite_table_index) < 2) {
        // Switching to opposite_table_index does not solve
        //   ambiguity problem.
        return(false);
      }
    }

    if (ivoldual_table.AmbiguousFacetBitsInUpperLifted(table_index) != 0) {
      // Lower lifted cube has an ambiguous facet.
      if (ivoldual_table.NumVerticesInUpperLifted(opposite_table_index) < 2) {
        // Switching to opposite_table_index does not solve
        //   ambiguity problem.
        return(false);
      }
    }

    return(true);
  }


  // Return true if ivol vertices in cube are candidates for splitting
  //   to avoid pair ambiguity.
  // Version which allows more splits.
  // Allow split if splitting will increase number of vertices
  //   in lower or upper lifted cube.
  //   (Does not necessarily need to do both.)
  bool is_cube_ambig_split_candidateC
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const GRID_CUBE_DATA & cube_data)
  {
    const TABLE_INDEX table_index = cube_data.table_index;
    const int numv_in_lower_lifted = 
      ivoldual_table.NumVerticesInLowerLifted(table_index);
    const int numv_in_upper_lifted = 
      ivoldual_table.NumVerticesInUpperLifted(table_index);
    const TABLE_INDEX opposite_table_index =
      ivoldual_table.OppositeTableIndex(table_index);

    if (ivoldual_table.NumAmbiguousFacets(table_index) > 1) 
      { return(false); }
    if ((numv_in_lower_lifted > 1) && (numv_in_upper_lifted > 1))
      { return(false); }

    return(true);
  }

  // Return true if pair of cubes sharing ambiguous facet
  //   have lower ambiguous facet and only one ivolv in each
  //   or have upper ambiguous facet and only one ivolv in each.
  bool does_cube_pair_have_ambiguity_problem
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const TABLE_INDEX table_indexA,
   const TABLE_INDEX table_indexB)
  {
    const int numv_in_lower_liftedA = 
      ivoldual_table.NumVerticesInLowerLifted(table_indexA);
    const int numv_in_upper_liftedA = 
      ivoldual_table.NumVerticesInUpperLifted(table_indexA);
    const int numv_in_lower_liftedB = 
      ivoldual_table.NumVerticesInLowerLifted(table_indexB);
    const int numv_in_upper_liftedB = 
      ivoldual_table.NumVerticesInUpperLifted(table_indexB);

    if (ivoldual_table.AmbiguousFacetBitsInLowerLifted(table_indexA) != 0) {
      if ((numv_in_lower_liftedA == 1) && (numv_in_lower_liftedB == 1)) 
        { return(true); }
    }

    if (ivoldual_table.AmbiguousFacetBitsInUpperLifted(table_indexA) != 0) {
      if ((numv_in_upper_liftedA == 1) && (numv_in_upper_liftedB == 1)) 
        { return(true); }
    }

    return(false);
  }


  // Return true if pair of cubes sharing ambiguous facet
  //   have lower ambiguous facet and only one ivolv in each
  //   or have upper ambiguous facet and only one ivolv in each.
  // - Version which allows more splits.
  bool does_cube_pair_have_ambiguity_problemD
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const TABLE_INDEX table_indexA,
   const TABLE_INDEX table_indexB)
  {
    const int numv_in_lower_liftedA = 
      ivoldual_table.NumVerticesInLowerLifted(table_indexA);
    const int numv_in_upper_liftedA = 
      ivoldual_table.NumVerticesInUpperLifted(table_indexA);
    const int numv_in_lower_liftedB = 
      ivoldual_table.NumVerticesInLowerLifted(table_indexB);
    const int numv_in_upper_liftedB = 
      ivoldual_table.NumVerticesInUpperLifted(table_indexB);

    if ((numv_in_lower_liftedA > 1) && (numv_in_lower_liftedB > 1)) 
      { return(false); }

    if ((numv_in_upper_liftedA > 1) && (numv_in_upper_liftedB > 1)) 
      { return(false); }

    if (ivoldual_table.AmbiguousFacetBitsInLowerLifted(table_indexA) != 0) {
      if ((numv_in_lower_liftedA == 1) && (numv_in_lower_liftedB == 1)) 
        { return(true); }
    }

    if (ivoldual_table.AmbiguousFacetBitsInUpperLifted(table_indexA) != 0) {
      if ((numv_in_upper_liftedA == 1) && (numv_in_upper_liftedB == 1)) 
        { return(true); }
    }

    return(false);
  }


  // Return true if splitting ivol vertices separates all 
  //   to avoid pair ambiguity.
  // Version which checks split pairs.
  bool is_cube_pair_ambig_split_candidate
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const GRID_CUBE_DATA & cube_dataA,
   const GRID_CUBE_DATA & cube_dataB)
  {
    const TABLE_INDEX table_indexA = cube_dataA.table_index;
    const TABLE_INDEX table_indexB = cube_dataB.table_index;
    const TABLE_INDEX opposite_table_indexA =
      ivoldual_table.OppositeTableIndex(table_indexA);
    const TABLE_INDEX opposite_table_indexB =
      ivoldual_table.OppositeTableIndex(table_indexB);

    if (ivoldual_table.NumAmbiguousFacets(table_indexA) > 1) 
      { return(false); }
    if (ivoldual_table.NumAmbiguousFacets(table_indexB) > 1) 
      { return(false); }

    if (does_cube_pair_have_ambiguity_problem
        (ivoldual_table, table_indexA, table_indexB)) {

      if (!does_cube_pair_have_ambiguity_problem
          (ivoldual_table, opposite_table_indexA, opposite_table_indexB)) 
        { return(true); }
    }

    return(false);
  }


  // Return true if splitting ivol vertices separates all 
  //   to avoid pair ambiguity.
  // - Version which checks split pairs.
  // - Version which allows more splits.
  bool is_cube_pair_ambig_split_candidateD
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const GRID_CUBE_DATA & cube_dataA,
   const GRID_CUBE_DATA & cube_dataB)
  {
    const TABLE_INDEX table_indexA = cube_dataA.table_index;
    const TABLE_INDEX table_indexB = cube_dataB.table_index;
    const TABLE_INDEX opposite_table_indexA =
      ivoldual_table.OppositeTableIndex(table_indexA);
    const TABLE_INDEX opposite_table_indexB =
      ivoldual_table.OppositeTableIndex(table_indexB);

    if (ivoldual_table.NumAmbiguousFacets(table_indexA) > 1) 
      { return(false); }
    if (ivoldual_table.NumAmbiguousFacets(table_indexB) > 1) 
      { return(false); }

    if (does_cube_pair_have_ambiguity_problemD
        (ivoldual_table, table_indexA, table_indexB)) {

      if (!does_cube_pair_have_ambiguity_problemD
          (ivoldual_table, opposite_table_indexA, opposite_table_indexB)) 
        { return(true); }
    }

    return(false);
  }

}

void IVOLDUAL::split_dual_ivolvert
(const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const std::vector<VERTEX_INDEX> & ivolpoly_cube, 
 const std::vector<POLY_VERTEX_INDEX> & poly_vertex,
 const IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 std::vector<GRID_CUBE_DATA> & cube_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list,
 std::vector<ISO_VERTEX_INDEX> & ivolpoly_vert,
 int & num_split)
{
  IJK::construct_dual_isovert_list(ivoldual_table, cube_list, ivolv_list);

  set_dual_ivolv_vertices
    (ivoldual_table, cube_list, ivolpoly_cube, poly_vertex, ivolpoly_info, 
     ivolpoly_vert);

  compute_num_splitB(ivoldual_table, cube_list, num_split);
}

namespace {

  // *** Local version of is_ambiguous_cube_facet_on_grid_boundary.
  // *** Required because IVOLDUAL_CUBE_TABLE_DOUBLE does not currently
  // *** support GetAmbiguousFacet() and does not set
  // *** ISODUAL_TABLE_ENTRY::index_of_some_ambiguous_facet.
  //
  // Return true if ambiguous cube facet is on grid boundary.
  // @pre Cube has one ambiguous facet.
  // @param[out] orth_dir Direction orthogonal to ambiguous facet.
  // @param[out] side  Side (0 or 1) of cube containing ambiguous facet.
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename CUBE_INDEX_TYPE, typename TABLE_INDEX_TYPE,
            typename DIR_TYPE, typename SIDE_TYPE>
  bool is_ambiguous_ivoldual_cube_facet_on_grid_boundary
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const CUBE_INDEX_TYPE cube_index,
   const TABLE_INDEX_TYPE table_index,
   DIR_TYPE & orth_dir,
   SIDE_TYPE & side)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename ISODUAL_TABLE::FACET_BITS_TYPE FACET_BITS_TYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const DTYPE dimension = isodual_table.Dimension();
    const NUM_TYPE num_cube_facets = isodual_table.NumPolyFacets();

    const FACET_BITS_TYPE facet_set = 
      isodual_table.AmbiguousFacetBits(table_index);
    const NUM_TYPE kf = get_first_one_bit(facet_set, num_cube_facets);
    orth_dir = IJK::cube_facet_orth_dir(dimension, kf);
    side = IJK::cube_facet_side(dimension, kf);

    return(grid.IsCubeFacetOnGridBoundary(cube_index, orth_dir, side));
  }

}

void IVOLDUAL::split_non_manifold_ivolv_pairs_ambig
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const VERTEX_INDEX index_to_cube_list[],
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  typedef typename DUALISO_GRID::NUMBER_TYPE NUM_TYPE;

  const int dimension = grid.Dimension();
  const NUM_TYPE num_cube_facets = IJK::compute_num_cube_facets(dimension);

  num_split = 0;
  for (NUM_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

    const TABLE_INDEX it0 = cube_list[i0].table_index;

    if (ivoldual_table.NumAmbiguousFacets(it0) == 1 &&
        is_cube_ambig_split_candidate(ivoldual_table, cube_list[i0])) {
      int orth_dir;
      int side;
      const VERTEX_INDEX cube_index0 = cube_list[i0].cube_index;
      
      if (is_ambiguous_ivoldual_cube_facet_on_grid_boundary
          (grid, ivoldual_table, cube_index0, it0, orth_dir, side)) {

        // Change table index of cube i0 opposite table
        //   indices to split interval volume vertex.
        cube_list[i0].table_index = 
          ivoldual_table.OppositeTableIndex(it0);
        num_split++;
      }
      else {
        const VERTEX_INDEX cube_index1 =
          grid.AdjacentVertex(cube_index0, orth_dir, side);
        const VERTEX_INDEX i1 = index_to_cube_list[cube_index1];
        if (is_cube_ambig_split_candidate(ivoldual_table, cube_list[i1])) {
          const TABLE_INDEX it1 = cube_list[i1].table_index;

          // Change table indices of cubes i0 and i1 to opposite table
          //   indices to split interval volume vertices.
          cube_list[i0].table_index = 
            ivoldual_table.OppositeTableIndex(it0);
          cube_list[i1].table_index = 
            ivoldual_table.OppositeTableIndex(it1);
          num_split += 2;
        }
      }
    }
  }
}


void IVOLDUAL::split_non_manifold_ivolv_pairs_ambig
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  IJK::ARRAY<VERTEX_INDEX> index_to_cube_list(grid.NumVertices());

  IJK::set_index_to_cube_list(cube_list, index_to_cube_list);

  split_non_manifold_ivolv_pairs_ambig
    (grid, ivoldual_table, index_to_cube_list.PtrConst(), 
     cube_list, num_split);
}


// *** Probably not necessary
void IVOLDUAL::split_non_manifold_ivolv_pairs_ambigB
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const VERTEX_INDEX index_to_cube_list[],
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  typedef typename DUALISO_GRID::NUMBER_TYPE NUM_TYPE;

  const int dimension = grid.Dimension();
  const NUM_TYPE num_cube_facets = IJK::compute_num_cube_facets(dimension);

  num_split = 0;
  for (NUM_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

    const TABLE_INDEX it0 = cube_list[i0].table_index;

    if (ivoldual_table.NumAmbiguousFacets(it0) == 1 &&
        is_cube_ambig_split_candidateB(ivoldual_table, cube_list[i0])) {
      int orth_dir;
      int side;
      const VERTEX_INDEX cube_index0 = cube_list[i0].cube_index;
      
      if (is_ambiguous_ivoldual_cube_facet_on_grid_boundary
          (grid, ivoldual_table, cube_index0, it0, orth_dir, side)) {

        // Change table index of cube i0 opposite table
        //   indices to split interval volume vertex.
        cube_list[i0].table_index = 
          ivoldual_table.OppositeTableIndex(it0);
        num_split++;
      }
      else {
        const VERTEX_INDEX cube_index1 =
          grid.AdjacentVertex(cube_index0, orth_dir, side);
        const VERTEX_INDEX i1 = index_to_cube_list[cube_index1];
        if (is_cube_ambig_split_candidateB(ivoldual_table, cube_list[i1])) {
          const TABLE_INDEX it1 = cube_list[i1].table_index;

          // Change table indices of cubes i0 and i1 to opposite table
          //   indices to split interval volume vertices.
          cube_list[i0].table_index = 
            ivoldual_table.OppositeTableIndex(it0);
          cube_list[i1].table_index = 
            ivoldual_table.OppositeTableIndex(it1);
          num_split += 2;
        }
      }
    }
  }
}


// *** Probably not necessary
void IVOLDUAL::split_non_manifold_ivolv_pairs_ambigB
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  IJK::ARRAY<VERTEX_INDEX> index_to_cube_list(grid.NumVertices());

  IJK::set_index_to_cube_list(cube_list, index_to_cube_list);

  split_non_manifold_ivolv_pairs_ambigB
    (grid, ivoldual_table, index_to_cube_list.PtrConst(), 
     cube_list, num_split);
}


void IVOLDUAL::split_non_manifold_ivolv_pairs_ambigC
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const VERTEX_INDEX index_to_cube_list[],
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  typedef typename DUALISO_GRID::NUMBER_TYPE NUM_TYPE;

  const int dimension = grid.Dimension();
  const NUM_TYPE num_cube_facets = IJK::compute_num_cube_facets(dimension);

  num_split = 0;
  for (NUM_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

    const TABLE_INDEX it0 = cube_list[i0].table_index;

    if (ivoldual_table.NumAmbiguousFacets(it0) == 1) {
      int orth_dir;
      int side;
      const VERTEX_INDEX cube_index0 = cube_list[i0].cube_index;

      const bool flag_on_boundary = 
        is_ambiguous_ivoldual_cube_facet_on_grid_boundary
        (grid, ivoldual_table, cube_index0, it0, orth_dir, side);
        
      if (flag_on_boundary) {
        if (is_cube_ambig_split_candidateC(ivoldual_table, cube_list[i0])) {

          // Change table index of cube i0 opposite table
          //   indices to split interval volume vertex.

          cube_list[i0].table_index = 
            ivoldual_table.OppositeTableIndex(it0);
          num_split++;
        }
      }
      else {
        const VERTEX_INDEX cube_index1 =
          grid.AdjacentVertex(cube_index0, orth_dir, side);
        const VERTEX_INDEX i1 = index_to_cube_list[cube_index1];

        if (is_cube_pair_ambig_split_candidate
            (ivoldual_table, cube_list[i0], cube_list[i1])) {
          const TABLE_INDEX it1 = cube_list[i1].table_index;

          // Change table indices of cubes i0 and i1 to opposite table
          //   indices to split interval volume vertices.
          cube_list[i0].table_index = 
            ivoldual_table.OppositeTableIndex(it0);
          cube_list[i1].table_index = 
            ivoldual_table.OppositeTableIndex(it1);
          num_split += 2;
        }
      }
    }
  }
}


void IVOLDUAL::split_non_manifold_ivolv_pairs_ambigC
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  IJK::ARRAY<VERTEX_INDEX> index_to_cube_list(grid.NumVertices());

  IJK::set_index_to_cube_list(cube_list, index_to_cube_list);

  split_non_manifold_ivolv_pairs_ambigC
    (grid, ivoldual_table, index_to_cube_list.PtrConst(), 
     cube_list, num_split);
}


// - Version which allows more splits than ambigC.
void IVOLDUAL::split_non_manifold_ivolv_pairs_ambigD
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const VERTEX_INDEX index_to_cube_list[],
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  typedef typename DUALISO_GRID::NUMBER_TYPE NUM_TYPE;

  const int dimension = grid.Dimension();
  const NUM_TYPE num_cube_facets = IJK::compute_num_cube_facets(dimension);

  num_split = 0;
  for (NUM_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

    const TABLE_INDEX it0 = cube_list[i0].table_index;

    if (ivoldual_table.NumAmbiguousFacets(it0) == 1) {
      int orth_dir;
      int side;
      const VERTEX_INDEX cube_index0 = cube_list[i0].cube_index;

      const bool flag_on_boundary = 
        is_ambiguous_ivoldual_cube_facet_on_grid_boundary
        (grid, ivoldual_table, cube_index0, it0, orth_dir, side);
        
      if (flag_on_boundary) {
        if (is_cube_ambig_split_candidateC(ivoldual_table, cube_list[i0])) {

          // Change table index of cube i0 opposite table
          //   indices to split interval volume vertex.

          cube_list[i0].table_index = 
            ivoldual_table.OppositeTableIndex(it0);
          num_split++;
        }
      }
      else {
        const VERTEX_INDEX cube_index1 =
          grid.AdjacentVertex(cube_index0, orth_dir, side);
        const VERTEX_INDEX i1 = index_to_cube_list[cube_index1];

        if (is_cube_pair_ambig_split_candidateD
            (ivoldual_table, cube_list[i0], cube_list[i1])) {
          const TABLE_INDEX it1 = cube_list[i1].table_index;

          // Change table indices of cubes i0 and i1 to opposite table
          //   indices to split interval volume vertices.
          cube_list[i0].table_index = 
            ivoldual_table.OppositeTableIndex(it0);
          cube_list[i1].table_index = 
            ivoldual_table.OppositeTableIndex(it1);
          num_split += 2;
        }
      }
    }
  }
}


void IVOLDUAL::split_non_manifold_ivolv_pairs_ambigD
(const DUALISO_GRID & grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 std::vector<GRID_CUBE_DATA> & cube_list,
 int & num_split)
{
  IJK::ARRAY<VERTEX_INDEX> index_to_cube_list(grid.NumVertices());

  IJK::set_index_to_cube_list(cube_list, index_to_cube_list);

  split_non_manifold_ivolv_pairs_ambigD
    (grid, ivoldual_table, index_to_cube_list.PtrConst(), 
     cube_list, num_split);
}


// **************************************************
// POSITION INTERVAL VOLUME VERTICES
// **************************************************

void IVOLDUAL::position_all_dual_ivol_vertices
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_TYPE * vertex_coord)
{
  const int dimension = scalar_grid.Dimension();

  if (dimension < 1) { return; }

  IJK::ARRAY<COORD_TYPE> coord0(dimension);
  IJK::ARRAY<COORD_TYPE> coord1(dimension);
  IJK::ARRAY<COORD_TYPE> coord2(dimension);
  CUBE_FACE_INFO cube(dimension);

  for (ISO_VERTEX_INDEX ivolv = 0; ivolv < ivolv_list.size(); ivolv++) {
    position_dual_ivolv_centroid_multi
      (scalar_grid, ivoldual_table, isovalue0, isovalue1, 
       ivolv_list[ivolv], cube, 
       vertex_coord+ivolv*dimension, 
       coord0.Ptr(), coord1.Ptr(), coord2.Ptr());
  }

}

void IVOLDUAL::position_all_dual_ivol_vertices
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const IJKDUAL::ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord)
{
  const int dimension = scalar_grid.Dimension();

  vertex_coord.resize(ivolv_list.size()*dimension);
  position_all_dual_ivol_vertices
    (scalar_grid, ivoldual_table, isodual_table, isovalue0, isovalue1, 
     ivolv_list, &(vertex_coord.front()));
}


void IVOLDUAL::position_dual_ivolv_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const DUAL_IVOLVERT & ivolv_info,
 const CUBE_FACE_INFO & cube,
 COORD_TYPE * vcoord, 
 COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, COORD_TYPE * temp_coord2)
{
  const int dimension = scalar_grid.Dimension();
  const int ivolv = ivolv_info.patch_index;
  const TABLE_INDEX it = ivolv_info.table_index;

  if (ivoldual_table.OnLowerIsosurface(it, ivolv)) {
    position_dual_ivolv_on_lower_isosurface_centroid_multi
      (scalar_grid, ivoldual_table, isovalue0, ivolv_info, cube,
       vcoord, temp_coord0, temp_coord1, temp_coord2);
  }
  else if (ivoldual_table.OnUpperIsosurface(it, ivolv)) {
    position_dual_ivolv_on_upper_isosurface_centroid_multi
      (scalar_grid, ivoldual_table, isovalue1, ivolv_info, cube,
       vcoord, temp_coord0, temp_coord1, temp_coord2);
  }
  else {
    position_dual_ivolv_in_interval_volume_centroid_multi
      (scalar_grid, ivoldual_table, isovalue0, isovalue1, ivolv_info, cube,
       vcoord, temp_coord0, temp_coord1, temp_coord2);
  }

}


void IVOLDUAL::position_dual_ivolv_on_lower_isosurface_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const SCALAR_TYPE isovalue,
 const DUAL_IVOLVERT & ivolv_info,
 const CUBE_FACE_INFO & cube,
 COORD_TYPE * vcoord, 
 COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, COORD_TYPE * temp_coord2)
{
  const int dimension = scalar_grid.Dimension();
  const VERTEX_INDEX icube = ivolv_info.cube_index;
  const int ivolv = ivolv_info.patch_index;
  const TABLE_INDEX table_index = ivolv_info.table_index;
  int num_intersected_edges = 0;
  IJK::set_coord(dimension, 0.0, vcoord);

  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    int k0 = cube.EdgeEndpoint(ie, 0);
    int k1 = cube.EdgeEndpoint(ie, 1);

    if (ivoldual_table.IsBelowIntervalVolume(table_index, k0) &&
        ivoldual_table.IsBelowIntervalVolume(table_index, k1)) {
      // Edge does not intersect lower isosurface.
      continue;
    }

    if (!ivoldual_table.IsBelowIntervalVolume(table_index, k0) &&
        !ivoldual_table.IsBelowIntervalVolume(table_index, k1)) {
      // Edge does not intersect lower isosurface.
      continue;
    }

    if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) {
      if (ivoldual_table.LowerIncident(table_index, ie) != ivolv) {
        // Dual polytope is not incident on vertex ivolv.
        continue;
      }
    }
    else {

      if (ivoldual_table.IsInIntervalVolume(table_index, k0)) {
        if (ivoldual_table.IncidentIVolVertex(table_index, k0) != ivolv) {
          // Dual polytope is not incident on vertex ivolv.
          continue;
        }
      }

      if (ivoldual_table.IsInIntervalVolume(table_index, k1)) {
        if (ivoldual_table.IncidentIVolVertex(table_index, k1) != ivolv) {
          // Dual polytope is not incident on vertex ivolv.
          continue;
        }
      }
    }

    VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
    VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);

    SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
    SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

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

  if (num_intersected_edges > 0) {
    IJK::multiply_coord
      (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
  }
  else {
    scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
  }

}


void IVOLDUAL::position_dual_ivolv_on_upper_isosurface_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const SCALAR_TYPE isovalue,
 const DUAL_IVOLVERT & ivolv_info,
 const CUBE_FACE_INFO & cube,
 COORD_TYPE * vcoord, 
 COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, COORD_TYPE * temp_coord2)
{
  const int dimension = scalar_grid.Dimension();
  const VERTEX_INDEX icube = ivolv_info.cube_index;
  const int ivolv = ivolv_info.patch_index;
  const TABLE_INDEX table_index = ivolv_info.table_index;
  int num_intersected_edges = 0;
  IJK::set_coord(dimension, 0.0, vcoord);

  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    int k0 = cube.EdgeEndpoint(ie, 0);
    int k1 = cube.EdgeEndpoint(ie, 1);

    if (ivoldual_table.IsAboveIntervalVolume(table_index, k0) &&
        ivoldual_table.IsAboveIntervalVolume(table_index, k1)) {
      // Edge does not intersect upper isosurface.
      continue;
    }

    if (!ivoldual_table.IsAboveIntervalVolume(table_index, k0) &&
        !ivoldual_table.IsAboveIntervalVolume(table_index, k1)) {
      // Edge does not intersect upper isosurface.
      continue;
    }

    if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) {
      if (ivoldual_table.UpperIncident(table_index, ie) != ivolv) {
        // Dual polytope is not incident on vertex ivolv.
        continue;
      }
    }
    else {

      if (ivoldual_table.IsInIntervalVolume(table_index, k0)) {
        if (ivoldual_table.IncidentIVolVertex(table_index, k0) != ivolv) {
          // Dual polytope is not incident on vertex ivolv.
          continue;
        }
      }

      if (ivoldual_table.IsInIntervalVolume(table_index, k1)) {
        if (ivoldual_table.IncidentIVolVertex(table_index, k1) != ivolv) {
          // Dual polytope is not incident on vertex ivolv.
          continue;
        }
      }
    }

    VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
    VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);

    SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
    SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

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

  if (num_intersected_edges > 0) {
    IJK::multiply_coord
      (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
  }
  else {
    scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
  }

}



namespace {

  /// Return true if interval volume polytope dual to iv0 and
  ///   incident on ivolv intersects grid edge (iv0,iv1) 
  ///   and some vertex of (iv0,iv1) is not in the interval volume.
  /// @param isovalueX Isovalue of intersection point.
  bool does_poly_dual_to_grid_vertex_intersect_grid_edge
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const TABLE_INDEX table_index,
   const int ivolv,
   const int ie,
   const int iv0,
   const int iv1,
   const SCALAR_TYPE s0,
   const SCALAR_TYPE s1,
   SCALAR_TYPE & isovalueX)
  {
    if (ivoldual_table.IsInIntervalVolume(table_index, iv0)) {
      if (ivoldual_table.IncidentIVolVertex(table_index, iv0) == ivolv) {
        if (ivoldual_table.IsBelowIntervalVolume(table_index, iv1)) {
          if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) 
            { isovalueX = (s0+isovalue0)/2.0; }
          else 
            { isovalueX = isovalue0; }
          return(true);
        }
        else if (ivoldual_table.IsAboveIntervalVolume(table_index, iv1)) {
          if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) 
            { isovalueX = (s0+isovalue1)/2.0; }
          else 
            { isovalueX = isovalue1; }
          return(true);
        }
        else {
          // iv0 and iv1 are in the interval volume.
          isovalueX = (s0+s1)/2.0;
          return(true);
        }
      }
    }

    return(false);
  }
   

  /// Return true if grid edge ie intersects interval volume polytope incident 
  /// on ivolv and some vertex of ie is not in the interval volume.
  /// @param isovalueX Isovalue of intersection point.
  bool does_grid_edge_intersect_incident_poly
  (const IVOLDUAL_CUBE_TABLE & ivoldual_table,
   const SCALAR_TYPE isovalue0,
   const SCALAR_TYPE isovalue1,
   const TABLE_INDEX table_index,
   const int ivolv,
   const int ie,
   const int iv0,
   const int iv1,
   const SCALAR_TYPE s0,
   const SCALAR_TYPE s1,
   SCALAR_TYPE & isovalueX)
  {

    if (does_poly_dual_to_grid_vertex_intersect_grid_edge
        (ivoldual_table, isovalue0, isovalue1, table_index, ivolv,
         ie, iv0, iv1, s0, s1, isovalueX))
      { return(true); }
    else if (does_poly_dual_to_grid_vertex_intersect_grid_edge
        (ivoldual_table, isovalue0, isovalue1, table_index, ivolv,
         ie, iv1, iv0, s1, s0, isovalueX))
      { return(true); }
    else if (ivoldual_table.EdgeHasDualIVolPoly(table_index, ie)) {
      if (ivoldual_table.UpperIncident(table_index, ie) == ivolv) {
        isovalueX = isovalue1;
        return(true);
      }
      else if (ivoldual_table.LowerIncident(table_index, ie) == ivolv) {
        isovalueX = isovalue0;
        return(true);
      }
    }

    return(false);
  }

}

void IVOLDUAL::position_dual_ivolv_in_interval_volume_centroid_multi
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 const SCALAR_TYPE isovalue0,
 const SCALAR_TYPE isovalue1,
 const DUAL_IVOLVERT & ivolv_info,
 const CUBE_FACE_INFO & cube,
 COORD_TYPE * vcoord, 
 COORD_TYPE * temp_coord0, COORD_TYPE * temp_coord1, COORD_TYPE * temp_coord2)
{
  const int dimension = scalar_grid.Dimension();
  const VERTEX_INDEX icube = ivolv_info.cube_index;
  const int ivolv = ivolv_info.patch_index;
  const TABLE_INDEX table_index = ivolv_info.table_index;
  SCALAR_TYPE isovalueX;

  int num_intersected_edges = 0;
  IJK::set_coord(dimension, 0.0, vcoord);

  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    const int k0 = cube.EdgeEndpoint(ie, 0);
    const int k1 = cube.EdgeEndpoint(ie, 1);
    const VERTEX_INDEX iend0 = scalar_grid.CubeVertex(icube, k0);
    const VERTEX_INDEX iend1 = scalar_grid.CubeVertex(icube, k1);
    const SCALAR_TYPE s0 = scalar_grid.Scalar(iend0);
    const SCALAR_TYPE s1 = scalar_grid.Scalar(iend1);

    if (does_grid_edge_intersect_incident_poly
        (ivoldual_table, isovalue0, isovalue1, table_index, ivolv,
         ie, k0, k1, s0, s1, isovalueX)) {

      scalar_grid.ComputeCoord(iend0, temp_coord0);
      scalar_grid.ComputeCoord(iend1, temp_coord1);

      if ((s0 < isovalueX && s1 < isovalueX) || 
          (s0 > isovalueX && s1 > isovalueX)) {
        // Use edge midpoint.
        IJK::linear_interpolate_coord
          (dimension, 0.5, temp_coord0, temp_coord1, temp_coord2);
      }
      else {
        IJK::linear_interpolate_coord
          (dimension, s0, temp_coord0, s1, temp_coord1, isovalueX, temp_coord2);
      }

      IJK::add_coord(dimension, vcoord, temp_coord2, vcoord);
      num_intersected_edges++;
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
