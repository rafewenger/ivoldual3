/// @file ijkdual_datastruct.cxx
/// @brief ijkdual data structures.
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

#include <assert.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "ijkcoord.tpp"

#include "ijkdual_types.h"
#include "ijkdual_datastruct.h"

using namespace IJK;
using namespace IJKDUAL;


// **************************************************
// CLASS DUALISO_DATA_FLAGS
// **************************************************

// Initialize DUALISO_DATA_FLAGS
void DUALISO_DATA_FLAGS::Init()
{
  flag_interval_volume = false;
  interpolation_type = LINEAR_INTERPOLATION;
  vertex_position_method = CENTROID_EDGE_ISO;
  allow_multiple_iso_vertices = true;
  flag_split_non_manifold = true;
  flag_separate_neg = true;
  flag_select_split = false;
  flag_connect_ambiguous = false;
  flag_iso_quad_dual_to_grid_edges = true;
  flag_dual_collapse = false;
  max_small_magnitude = 0.0001;
  use_triangle_mesh = false;
  quad_tri_method = UNDEFINED_TRI;
  quad_edge_intersection_method = QEI_INTERPOLATE_SCALAR;

  flag_tri4_quad = false;
  flag_tri5_pentagon = false;
  tri4_position_method = TRI4_ON_GRID_EDGE;
  min_distance_use_tri4 = 0.4;
  min_distance_allow_tri4 = 0.1;
}


// Set
void DUALISO_DATA_FLAGS::Set(const DUALISO_DATA_FLAGS & data_flags)
{
  *this = data_flags;
}

// Return false if flag_iso_quad_dual_to_grid_edges is false.
bool DUALISO_DATA_FLAGS::CheckIsoQuadDualToGridEdges
(const char * proc_name, IJK::ERROR & error) const
{
  if (flag_iso_quad_dual_to_grid_edges) {
    return (true);
  }
  else {
    error.AddMessage
      ("Programming error.  Isosurface quadrilaterals may not be dual to grid edges.");
    error.AddMessage
      ("  Procedure ", proc_name, " requires ");
    error.AddMessage
      ("  that all isosurface quadrilaterals are dual to grid edges.");
    return(false);
  }
}

// **************************************************
// DUALISO TIME
// **************************************************

IJKDUAL::DUALISO_TIME::DUALISO_TIME()
// constructor
{
  Clear();
}

void IJKDUAL::DUALISO_TIME::Clear()
{
  preprocessing = 0.0;
  extract = 0.0;
  merge = 0.0;
  position = 0.0;
  collapse = 0.0;
  total = 0.0;
}

void IJKDUAL::DUALISO_TIME::Add(const DUALISO_TIME & dualiso_time)
{
  preprocessing += dualiso_time.preprocessing;
  extract += dualiso_time.extract;
  merge += dualiso_time.merge;
  position += dualiso_time.position;
  collapse += dualiso_time.collapse;
  total += dualiso_time.total;
}

// **************************************************
// INFO CLASSES
// **************************************************

IJKDUAL::GRID_INFO::GRID_INFO()
{
  Clear();
}

void IJKDUAL::GRID_INFO::Clear()
{
  num_cubes = 0;
}

void IJKDUAL::SCALAR_INFO::Init(const int dimension)
{
  this->dimension = 0;
  SetDimension(dimension);

  Clear();
}

void IJKDUAL::SCALAR_INFO::FreeAll()
{
  dimension = 0;
}

void IJKDUAL::SCALAR_INFO::Clear()
{
  num_non_empty_cubes = 0;
  num_bipolar_edges = 0;
}

void IJKDUAL::SCALAR_INFO::SetDimension(const int dimension)
{
  FreeAll();

  this->dimension = dimension;

  Clear();
}

void IJKDUAL::SCALAR_INFO::Copy(const SCALAR_INFO & info)
{
  Init(info.Dimension());
  num_non_empty_cubes = info.num_non_empty_cubes;
  num_bipolar_edges = info.num_bipolar_edges;
}

/// Copy assignment.
const SCALAR_INFO &  IJKDUAL::SCALAR_INFO::operator =
(const SCALAR_INFO & right)
{
  if (&right != this) {
    FreeAll();
    Copy(right);
  }

  return *this;
}

IJKDUAL::SCALAR_INFO::~SCALAR_INFO()
{
  dimension = 0;
  Clear();
}

IJKDUAL::MULTI_ISOV_INFO::MULTI_ISOV_INFO()
{
  Clear();
}

void IJKDUAL::MULTI_ISOV_INFO::Clear()
{
  num_cubes_single_isov = 0;
  num_cubes_multi_isov = 0;
  num_non_manifold_split = 0;
  num_1_2_changed = 0;
  num_connect_changed = 0;
  num_ambig_ridge_cubes_changed = 0;
  num_non_ambig_ridge_cubes_changed = 0;
}

IJKDUAL::DUALISO_INFO::DUALISO_INFO()
{
  Clear();
}

IJKDUAL::DUALISO_INFO::DUALISO_INFO(const int dimension):scalar(dimension)
{
  Clear();
}

void IJKDUAL::DUALISO_INFO::Clear()
{
  grid.Clear();
  scalar.Clear();
  time.Clear();
  multi_isov.Clear();
}


// **************************************************
// MERGE DATA
// **************************************************

void IJKDUAL::MERGE_DATA::Init
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MERGE_INDEX num_obj_per_vertex, const MERGE_INDEX num_obj_per_edge)
{
  this->num_obj_per_vertex = num_obj_per_vertex;
  this->num_obj_per_edge = num_obj_per_edge;
  this->num_obj_per_grid_vertex = 
    dimension*num_obj_per_edge + num_obj_per_vertex;
  compute_num_grid_vertices(dimension, axis_size, num_vertices);
  num_edges = dimension*num_vertices;
  vertex_id0 = num_obj_per_edge*num_edges;
  MERGE_INDEX num_obj = 
    num_obj_per_vertex*num_vertices + num_obj_per_edge*num_edges;
  INTEGER_LIST<MERGE_INDEX,MERGE_INDEX>::Init(num_obj);
}

bool IJKDUAL::MERGE_DATA::Check(ERROR & error) const
{
  if (MaxNumInt() < 
      NumObjPerVertex()*NumVertices() + NumObjPerEdge()*NumEdges()) {
    error.AddMessage("Not enough allocated memory.");
    return(false);
  };

  return(true);
}

