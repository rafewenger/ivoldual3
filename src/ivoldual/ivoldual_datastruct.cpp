/// @file ivoldual_datastruct.cpp
/// @brief ivoldual data structures.
/// - Version 0.3.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2026 Rephael Wenger

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

#include "ivoldual_datastruct.h"
#include "ijkgrid_macros.h"
#include "ivoldual_ivolpoly.tpp"
#include "ivoldualtable.h"

// **************************************************
// CLASS IVOLDUAL_DATA_FLAGS
// **************************************************

// Initialize
void IVOLDUAL::IVOLDUAL_DATA_FLAGS::Init()
{
  const COORD_TYPE ONE_THIRD(1.0/3.0);

  hex_tri_method = UNDEFINED_HEX_TRI;
  flag_split_ambig_pairs = false;
  flag_split_ambig_pairsB = false;
  flag_split_ambig_pairsC = false;
  flag_split_ambig_pairsD = false;
  flag_rm_non_manifold = false;
  flag_add_isov_dual_to_hexahedra = false;
  flag_orient_in = false;
  flag_lsmooth_elength = false;
  flag_lsmooth_jacobian = false;
  flag_gsmooth_jacobian = false;
  flag_split_hex = false;
  flag_collapse_hex = false;
  flag_set_interior_code_from_scalar = false;
  default_interior_code = 2;
  lsmooth_elength_iter = 1;
  lsmooth_jacobian_iter = 1;
  gsmooth_jacobian_iter = 1;
  elength_threshold = 0.1;
  jacobian_threshold = 0.0;
  split_hex_threshold = 0.0;
  collapse_hex_threshold = 0.0;

  flag_expand_thin_regions = false;
  thin_separation_distance = ONE_THIRD;
}

// **************************************************
// CLASS IVOLDUAL_POLY_INFO
// **************************************************

// Initialize
void IVOLDUAL::IVOLDUAL_POLY_INFO::Init()
{
  flag_subdivide_hex = false;
}

// **************************************************
// CLASS DUALISO INFO MEMBER FUNCTIONS
// **************************************************
IVOLDUAL::IVOLDUAL_INFO::IVOLDUAL_INFO() 
{
  Clear();
}

IVOLDUAL::IVOLDUAL_INFO::IVOLDUAL_INFO(const int dimension) 
{
  Clear();
}

void IVOLDUAL::IVOLDUAL_INFO::Clear() 
{
  num_non_manifold_changes = 0;
  num_vertices_moved_in_expand_thin = 0;
}

// **************************************************
// CLASS IVOLDUAL_DATA MEMBER FUNCTIONS
// **************************************************

// Copy, subsample, supersample or subdivide scalar grid.
void IVOLDUAL::IVOLDUAL_DATA::SetScalarGrid
(const DUALISO_SCALAR_GRID_BASE & scalar_grid2, 
 const bool flag_subsample, const int subsample_resolution,
 const bool flag_supersample, const int supersample_resolution, 
 const bool flag_subdivide, const bool flag_rm_diag_ambig,
 const bool flag_add_outer_layer, 
 const GRID_VERTEX_ENCODING interior_code,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1)
{
  IJK::PROCEDURE_ERROR error("IVOLDUAL_DATA::SetScalarGrid");

  if (flag_subsample + flag_supersample + flag_subdivide > 1) {
    error.AddMessage
      ("At most one of -subsample, -supersample or -subdivide may be used.");
    throw error;
  }

  // Set interior code
  default_interior_code = interior_code;
  
  if (flag_subsample) {
    // subsample grid
    SubsampleScalarGrid(scalar_grid2, subsample_resolution);
  }
  else if (flag_supersample) {
    // supersample grid
    SupersampleScalarGrid(scalar_grid2, supersample_resolution);
  }
  else if (flag_subdivide) {
    // subdivide grid
    int subdivide_resolution(2);
    SupersampleScalarGrid(scalar_grid2, subdivide_resolution);
    SubdivideScalarGrid(isovalue0, isovalue1);
  }
  else {
    CopyScalarGrid(scalar_grid2);
  };
  // Eliminate non-manifold of diagonal '++' corner.
  if (flag_rm_diag_ambig) {
    int recur = 0;
    while (RmDiagonalAmbig(isovalue0, isovalue1, 0)) {
      SubdivideScalarGrid(isovalue0, isovalue1);
      if (++recur == 10) {
        printf("Too many iterations in -rm_diag_ambig, break!\n");
        break;
      }
    }
    RmDiagonalAmbig(isovalue0, isovalue1, 1);
  }
  // Add outer layer to the scalar grid.
  if (flag_add_outer_layer) {
    scalar_grid.AddOuterLayer();
  }
}

void IVOLDUAL::IVOLDUAL_DATA::SubdivideScalarGrid
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1) 
{

  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  int dx = axis_size[0], dy = axis_size[1], dz = axis_size[2];
  int dxy = dx * dy;

  for (int i = 0; i < scalar_grid.NumVertices(); i++) {
    int x = i % dx;
    int y = (i / dx) % dy;
    int z = i / dx / dy;

    // Facet center in X-Z plane
    if (x % 2 == 1 && y % 2 == 0 && z % 2 == 1) {
      int corner[] = {i-dxy-1, i-dxy+1, i+dxy+1, i+dxy-1};
      int edge[] = {i-dxy, i+1, i+dxy, i-1};
      EvaluateSubdivideCenter(corner, edge, i, isovalue0, isovalue1);
    }
    // Facet center in Y-Z plane
    else if (x % 2 == 0 && y % 2 == 1 && z % 2 == 1) {
      int corner[] = {i-dxy-dx, i-dxy+dx, i+dxy+dx, i+dxy-dx};
      int edge[] = {i-dxy, i+dx, i+dxy, i-dx};
      EvaluateSubdivideCenter(corner, edge, i, isovalue0, isovalue1);
    }
    // Facet center in X-Y plane
    else if (x % 2 == 1 && y % 2 == 1 && z % 2 == 0) {
      int corner[] = {i-dx-1, i-dx+1, i+dx+1, i+dx-1};
      int edge[] = {i-dx, i+1, i+dx, i-1};
      EvaluateSubdivideCenter(corner, edge, i, isovalue0, isovalue1);
    }
  }
  EvaluateCubeCenter(isovalue0, isovalue1);
}

void IVOLDUAL::IVOLDUAL_DATA::EvaluateCubeCenter
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1) 
{
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  int dx = axis_size[0], dy = axis_size[1], dz = axis_size[2];
  int dxy = dx * dy;

  for (int i = 0; i < scalar_grid.NumVertices(); i++) {
    int x = i % dx;
    int y = (i / dx) % dy;
    int z = i / dx / dy;

    std::vector<SCALAR_TYPE> non_manifold_list;
    // Vertex at cube center
    if (x % 2 == 1 && y % 2 == 1 && z % 2 == 1) {
      SCALAR_TYPE minus_val = 0.0, equal_val = 0.0, plus_val = 0.0;

      int cornerXY[] = {i-dx-1, i-dx+1, i+dx+1, i+dx-1};
      int edgeXY[] = {i-dx, i+1, i+dx, i-1};      

      int cornerXZ[] = {i-dxy-1, i-dxy+1, i+dxy+1, i+dxy-1};
      int edgeXZ[] = {i-dxy, i+1, i+dxy, i-1};

      int cornerYZ[] = {i-dxy-dx, i-dxy+dx, i+dxy+dx, i+dxy-dx};
      int edgeYZ[] = {i-dxy, i+dx, i+dxy, i-dx};

      // Evaluate cube center in X-Y plane
      EvaluateSubdivideCenter(cornerXY, edgeXY, i, isovalue0, isovalue1);
      if (CheckManifold(cornerXY, edgeXY, i, isovalue0, isovalue1) &&
          CheckManifold(cornerXZ, edgeXZ, i, isovalue0, isovalue1) &&
          CheckManifold(cornerYZ, edgeYZ, i, isovalue0, isovalue1)) {
        continue;
      }
      else {
        if (symbol(i, isovalue0, isovalue1) == 0) 
          { minus_val = scalar_grid.Scalar(i); }
        else if (symbol(i, isovalue0, isovalue1) == 2)
          { equal_val = scalar_grid.Scalar(i); }
        else if (symbol(i, isovalue0, isovalue1) == 3)
          { plus_val = scalar_grid.Scalar(i); }
      }

      // Evaluate cube center in X-Z plane
      EvaluateSubdivideCenter(cornerXZ, edgeXZ, i, isovalue0, isovalue1);
      if (CheckManifold(cornerXY, edgeXY, i, isovalue0, isovalue1) &&
          CheckManifold(cornerXZ, edgeXZ, i, isovalue0, isovalue1) &&
          CheckManifold(cornerYZ, edgeYZ, i, isovalue0, isovalue1)) {
        continue;
      }
      else {
        if (symbol(i, isovalue0, isovalue1) == 0) 
          { minus_val = scalar_grid.Scalar(i); }
        else if (symbol(i, isovalue0, isovalue1) == 2)
          { equal_val = scalar_grid.Scalar(i); }
        else if (symbol(i, isovalue0, isovalue1) == 3)
          { plus_val = scalar_grid.Scalar(i); }
      }

      EvaluateSubdivideCenter(cornerYZ, edgeYZ, i, isovalue0, isovalue1);
      if (CheckManifold(cornerXY, edgeXY, i, isovalue0, isovalue1) &&
          CheckManifold(cornerXZ, edgeXZ, i, isovalue0, isovalue1) &&
          CheckManifold(cornerYZ, edgeYZ, i, isovalue0, isovalue1)) {
        continue;
      }
      else {
        if (symbol(i, isovalue0, isovalue1) == 0) 
          { minus_val = scalar_grid.Scalar(i); }
        else if (symbol(i, isovalue0, isovalue1) == 2)
          { equal_val = scalar_grid.Scalar(i); }
        else if (symbol(i, isovalue0, isovalue1) == 3)
          { plus_val = scalar_grid.Scalar(i); }
      }

      // Set the cube center value for rule confliction cases.
      if (equal_val > 0 && plus_val > 0) {
        scalar_grid.Set(i, plus_val);
      }
      else if (equal_val > 0 && minus_val > 0) {
        scalar_grid.Set(i, minus_val);
      }
    }
  }
}

bool IVOLDUAL::IVOLDUAL_DATA::EvaluateSubdivideCenter
(int corner[], int edge[], int icenter,
 const SCALAR_TYPE v0, const SCALAR_TYPE v1)
{
  const int NUM_SIDES(4);

  // Subdivide Rule 4
  for (int i = 0; i < NUM_SIDES/2; i++) {
    int cur = edge[i], oppo = edge[i+2];
    if (symbol(cur, v0, v1) == symbol(oppo, v0, v1)) {
      SCALAR_TYPE val = 
          0.5*(scalar_grid.Scalar(cur) + scalar_grid.Scalar(oppo));
      scalar_grid.Set(icenter, val);
      return true;
    }
  }

  // Count number of symbols +/=/-
  int num_plus = 0, num_minus = 0, num_equal = 0;
  for (int i = 0; i < NUM_SIDES; i++) {
    int cur = corner[i];
    if (scalar_grid.Scalar(cur) < v0) num_minus++;
    else if (scalar_grid.Scalar(cur) > v1) num_plus++;
    else num_equal++;
  }

  // Subdivide Rule 2
  if (num_plus == 3 || num_equal == 3 || num_minus == 3) {
    for (int i = 0; i < NUM_SIDES; i++) {
      int cur = corner[i];
      int left = corner[(i+NUM_SIDES-1)%NUM_SIDES];
      int right = corner[(i+1)%NUM_SIDES];

      if (symbol(cur, v0, v1) == symbol(left, v0, v1) &&
          symbol(cur, v0, v1) == symbol(right, v0, v1))
      {
        SCALAR_TYPE val = 
          0.5*(scalar_grid.Scalar(left) + scalar_grid.Scalar(right));
        scalar_grid.Set(icenter, val);
        return true;
      }
    }
  }
  // Subdivide Rule 3 and Rule 5
  else if (num_plus == 2 || num_equal == 2 || num_minus == 2) {
    bool flag_rule_five = false;
    // Subdivide Rule 3
    for (int i = 0; i < NUM_SIDES; i++) {
      int cur = corner[i], oppo = corner[(i+2)%4];

      if (symbol(cur, v0, v1) == symbol(oppo, v0, v1)) {
        // Scalar grid index of vertex on edge
        int e1 = edge[i], e2 = edge[(i+1)%NUM_SIDES], 
            e3 = edge[(i+2)%NUM_SIDES], e4 = edge[(i+3)%NUM_SIDES];

        if ((symbol(cur, v0, v1) == symbol(e1, v0, v1) &&
             symbol(cur, v0, v1) == symbol(e2, v0, v1))
            ||
            (symbol(cur, v0, v1) == symbol(e3, v0, v1) &&
             symbol(cur, v0, v1) == symbol(e4, v0, v1)))
        {
          
          SCALAR_TYPE val = 
            0.5*(scalar_grid.Scalar(cur) + scalar_grid.Scalar(oppo));
          scalar_grid.Set(icenter, val);  
          return true;       
        }  
        flag_rule_five = true;
      }
    }
    // Subdivide Rule 5
    if (flag_rule_five) {
      SCALAR_TYPE val = 0.5*(v0 + v1);
      scalar_grid.Set(icenter, val);
      return true; 
    }
  }
  return false;
}

int IVOLDUAL::IVOLDUAL_DATA::symbol
(const int cur, const SCALAR_TYPE v0, const SCALAR_TYPE v1)
{
  int mark;
  if (scalar_grid.Scalar(cur) < v0) mark = 0;
  else if (scalar_grid.Scalar(cur) > v1) mark = 3;
  else mark = 2;
  return mark;
}

bool IVOLDUAL::IVOLDUAL_DATA::CheckManifold
(int corner[], int edge[], int icenter,
 const SCALAR_TYPE v0, const SCALAR_TYPE v1)
{
  const int NUM_SIDES(4);
  for (int i = 0; i < NUM_SIDES; i++) {
    int cur = corner[i];
    int left = edge[(i+NUM_SIDES-1)%NUM_SIDES];
    int right = edge[i];

    if (symbol(icenter, v0, v1) == symbol(cur, v0, v1) &&
      symbol(left, v0, v1) != symbol(cur, v0, v1) &&
      symbol(right, v0, v1) != symbol(cur, v0, v1)) 
    {
      if (symbol(left, v0, v1) == symbol(right, v0, v1) || 
          symbol(cur, v0, v1) != 2)
        { return false; }
    }

    if (symbol(left, v0, v1) == symbol(right, v0, v1) &&
      symbol(left, v0, v1) != symbol(cur, v0, v1) &&
      symbol(left, v0, v1) != symbol(icenter, v0, v1) &&
      symbol(left, v0, v1) != 2) 
    {
      return false;
    }
  }
  return true;
}

void IVOLDUAL::IVOLDUAL_DATA::EliminateAmbigFacets
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1, 
 VERTEX_INDEX & changes_of_ambiguity)
{
  int last_changes = 1;
  while (last_changes) {
    // Eliminate non-manifold cases caused by opposite diagonal plus vertices.
    last_changes = RmDiagonalAmbig(isovalue0, isovalue1, 0);
    // Eliminate non-manifold cases caused by ambiguous facets.
    last_changes += scalar_grid.EliminateAmbiguity(isovalue0, isovalue1);
    changes_of_ambiguity += last_changes;
  }
  scalar_grid.AddOuterLayer();
}

int IVOLDUAL::IVOLDUAL_DATA::RmDiagonalAmbig
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1, 
 int caseID)
{
  int changes_of_cube = 0;
  const int DIM3(3);
  IVOLDUAL_CUBE_TABLE ivoldual_table(DIM3, true);
  TABLE_INDEX table_index;

  IJK_FOR_EACH_GRID_CUBE(icube, scalar_grid, VERTEX_INDEX) {

    int idx_p1 = -1, idx_p2 = -1, num_plus = 0;
    int idx_m1 = -1, idx_m2 = -1, num_minus = 0;

    for (int k = 0; k < scalar_grid.NumCubeVertices(); k++) {
      int idx = scalar_grid.CubeVertex(icube, k);
      if (scalar_grid.Scalar(idx) > isovalue1) {
        num_plus++;
        if (num_plus == 1) idx_p1 = k;
        else if (num_plus == 2) idx_p2 = k;
      }
      else if (scalar_grid.Scalar(idx) < isovalue0) {
        num_minus++;
        if (num_minus == 1) idx_m1 = k;
        else if (num_minus == 2) idx_m2 = k;
      }
    }

    if (num_plus == 2 && idx_p1 + idx_p2 == 7 && default_interior_code == 2)  {
      int num_change = 0;
      scalar_grid.EliminateDiagonalPlus
        (isovalue0, isovalue1, icube, caseID, num_change);
      changes_of_cube += num_change;
    }
    else if (num_minus == 2 && idx_m1 + idx_m2 == 7 && default_interior_code == 1) {
      int num_change = 0;
      scalar_grid.EliminateDiagonalMinus
        (isovalue0, isovalue1, icube, caseID, num_change);
      changes_of_cube += num_change;
    }
  }
  return changes_of_cube;
}

/// Eliminate ambiguity facets in supersample grid.
int IVOLDUAL::IVOLDUAL_SCALAR_GRID::EliminateAmbiguity
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1) 
{
  const DTYPE dim = this->dimension;
  int changes_of_ambiguity = 0;

  ATYPE sizex = this->AxisSize(0);
  ATYPE sizey = this->AxisSize(1);
  ATYPE sizez = this->AxisSize(2);


  for (VTYPE i = 0; i < this->NumVertices(); i++) {
    
    for (DTYPE d = 0; d < dim; d++) {

      VTYPE v01_idx = this->NextVertex(i, d);
      VTYPE v10_idx = this->NextVertex(i, (d+1)%dim);
      VTYPE v11_idx = this->NextVertex(this->NextVertex(i, d), (d+1)%dim);

      // If current vertex is at the end of a direction, then continue.
      if ((d == 0 && v01_idx%sizex == 0) ||
          (d == 2 && v10_idx%sizex == 0) ||
          (d == 1 && (v01_idx/sizex)%sizey == 0) ||
          (d == 0 && (v10_idx/sizex)%sizey == 0) ||
          (d == 2 && (v01_idx/sizex/sizey)%sizez == 0) ||
          (d == 1 && (v10_idx/sizex/sizey)%sizez == 0) ) 
        { continue; } 

      // Get scalar values of a square.
      STYPE *s00 = &this->scalar[i];
      STYPE *s01 = &this->scalar[v01_idx];
      STYPE *s10 = &this->scalar[v10_idx];
      STYPE *s11 = &this->scalar[v11_idx];


      // In a unit cube before subdivide:
      // Original vertex is level 0
      // Vertex at unit cube edge center is level 1
      // Vertex at unit cube facet center is level 2
      // Vertex at unit cube body center is level 3
      if (*s00 < isovalue0 && *s01 > isovalue0 && *s10 > isovalue0 && *s11 < isovalue0) {
        
        VTYPE z_level = v11_idx/sizex/sizey;
        VTYPE y_level = (v11_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v11_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s00 = 0.5 * (isovalue0 + isovalue1); 
        else *s11 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
      else if (*s00 > isovalue0 && *s01 < isovalue0 && *s10 < isovalue0 && *s11 > isovalue0) {
        
        VTYPE z_level = v01_idx/sizex/sizey;
        VTYPE y_level = (v01_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v01_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s10 = 0.5 * (isovalue0 + isovalue1); 
        else *s01 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
      else if (*s00 > isovalue1 && *s01 < isovalue1 && *s10 < isovalue1 && *s11 > isovalue1) {
        
        VTYPE z_level = v11_idx/sizex/sizey;
        VTYPE y_level = (v11_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v11_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s00 = 0.5 * (isovalue0 + isovalue1); 
        else *s11 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
      else if (*s00 < isovalue1 && *s01 > isovalue1 && *s10 > isovalue1 && *s11 < isovalue1) {
        
        VTYPE z_level = v01_idx/sizex/sizey;
        VTYPE y_level = (v01_idx%(sizex*sizey))/sizex;
        VTYPE x_level = v01_idx%sizex;
        VTYPE level = z_level%2 + y_level%2 + x_level%2;

        if (level == 0) *s10 = 0.5 * (isovalue0 + isovalue1); 
        else *s01 = 0.5 * (isovalue0 + isovalue1); 
        changes_of_ambiguity++; 
      }
    }  
  }

  return changes_of_ambiguity;
}

/// Add outer layer to subdivide grid
void IVOLDUAL::IVOLDUAL_SCALAR_GRID::AddOuterLayer() 
{

  const DTYPE dim = this->Dimension();
  
  int dim_x = this->AxisSize(0);
  int dim_y = this->AxisSize(1);
  int dim_z = this->AxisSize(2);

  AXIS_SIZE_TYPE axis_size[3] = {dim_x, dim_y, dim_z};

  IVOLDUAL_SCALAR_GRID scalar_grid3;
  scalar_grid3.SetSize(3, axis_size);
  scalar_grid3.CopyScalar(*this);

  ATYPE dim_x2 = dim_x + 2;
  ATYPE dim_y2 = dim_y + 2;
  ATYPE dim_z2 = dim_z + 2;

  AXIS_SIZE_TYPE axis_size2[3] = {dim_x2, dim_y2, dim_z2};
  this->SetSize(dim, axis_size2);

  VERTEX_INDEX iv = 0;
  VERTEX_INDEX iw = 0;

  for(int k = 0; k < dim_z; k++) {
    for(int j = 0; j < dim_y; j++) {
      for(int i = 0; i < dim_x; i++) {
        this->Set(iw++, scalar_grid3.Scalar(iv));
        if (i == 0 || i+1 == dim_x)
          this->Set(iw++, scalar_grid3.Scalar(iv));
        iv++;
      }
      if (j == 0 || j+1 == dim_y) {
        for (int i = 0; i < dim_x2; i++) {
          this->Set(iw, this->Scalar(iw-dim_x2));
          iw++;
        }
      }
    }

    if (k == 0 || k+1 == dim_z) {
      for (int j = 0; j < dim_y2; j++) {
        for (int i = 0; i < dim_x2; i++) {
          this->Set(iw, this->Scalar(iw-dim_x2*dim_y2));
          iw++;
        }
      }
    }
  }
}

void IVOLDUAL::IVOLDUAL_SCALAR_GRID::EliminateDiagonalPlus
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 VTYPE icube, int caseID, int& num_change) 
{
  int dx = this->axis_size[0];
  int dy = this->axis_size[1];
  int dz = this->axis_size[2];

  for (int k = 0; k < this->NumCubeVertices(); k++) {
    int idx = this->CubeVertex(icube, k);
    int x = idx % dx;
    int y = (idx / dx) % dy;
    int z = idx / dx / dy;

    // Vertex is at cube edge center
    if (caseID == 0 && this->Scalar(idx) > isovalue1 && 
        (x+y+z)%2 == 1 && (x%2)+(y%2)+(z%2) < 3) {    
      this->Set(idx, 0.5*(isovalue0 + isovalue1));
      num_change = 1;
      break;
    }
    // Vertex is at the cube center
    else if (caseID == 1 && this->Scalar(idx) > isovalue1 && 
        (x+y+z)%2 == 1 && (x%2)+(y%2)+(z%2) == 3) {      
      this->Set(idx, 0.5*(isovalue0 + isovalue1));
      num_change = 1;
      break;
    }
  }
}

void IVOLDUAL::IVOLDUAL_SCALAR_GRID::EliminateDiagonalMinus
(const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 VTYPE icube, int caseID, int& num_change) 
{
  int dx = this->axis_size[0];
  int dy = this->axis_size[1];
  int dz = this->axis_size[2];

  for (int k = 0; k < this->NumCubeVertices(); k++) {
    int idx = this->CubeVertex(icube, k);
    int x = idx % dx;
    int y = (idx / dx) % dy;
    int z = idx / dx / dy;

    // Vertex is at cube edge center
    if (caseID == 0 && this->Scalar(idx) < isovalue0 && 
        (x+y+z)%2 == 1 && (x%2)+(y%2)+(z%2) < 3) {    
      this->Set(idx, 0.5*(isovalue0 + isovalue1));
      num_change = 1;
      break;
    }
    // Vertex is at the cube center
    else if (caseID == 1 && this->Scalar(idx) < isovalue0 && 
        (x+y+z)%2 == 1 && (x%2)+(y%2)+(z%2) == 3) {      
      this->Set(idx, 0.5*(isovalue0 + isovalue1));
      num_change = 1;
      break;
    }
  }
}


// **************************************************
// DUAL_IVOLVERT MEMBER FUNCTIONS
// **************************************************

void IVOLDUAL::DUAL_IVOLVERT::Init()
{
  // Data from interval volume lookup table.
  num_incident_hex = 0;
  num_incident_iso_quad = 0;
  flag_lower_isosurface = false;
  flag_upper_isosurface = false;

  separation_vertex = 0;
  separation_edge_direction = 0;
  flag_missing_ivol_hexahedra = false;
  is_doubly_connected = false;
  in_loop = false;
  in_box = false;
  in_pseudobox = false;
  flag_adjacent_to_lower_isosurface = false;
  flag_adjacent_to_upper_isosurface = false;
}
