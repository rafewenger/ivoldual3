/// @file ivoldual_reposition.cxx
/// @brief Reposition vertices for mesh quality optimization.
/// - Version 0.3.0

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

#include "ivoldual_compute.h"
#include "ivoldual_reposition.h"
#include "ivoldual_divide_hex.h"
#include "ijktriangulate.txx"

using namespace IJK;
using namespace IVOLDUAL;

void IVOLDUAL::eliminate_non_manifold_grid
(IVOLDUAL_DATA & ivoldual_data, 
 const SCALAR_TYPE isovalue0, 
 const SCALAR_TYPE isovalue1, 
 IVOLDUAL_INFO & dualiso_info) 
{
  int num_changes = 0;
  ivoldual_data.EliminateAmbigFacets(isovalue0, isovalue1, num_changes);
  dualiso_info.num_non_manifold_changes = num_changes;
}

void IVOLDUAL::laplacian_smooth_elength
(const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord,
 float laplacian_smooth_limit, 
 int iteration)
{
  const int d = 3;
  float dist;
  COORD_TYPE * vcoord = &(vertex_coord.front());

  for (int it = 0; it < 2*iteration+1; it++) {

    bool skipSurfaceVert = (it % 2 == 0);

    // Loop over all vertices
    for (int cur = 0; cur < vertex_adjacency_list.NumVertices(); cur++) {
      
      // Current node coordinates.
      COORD_TYPE *cur_coord = vcoord + cur*d;

      // Check if current node is on isosurface.
      const int ivolv_cur = ivolv_list[cur].patch_index;
      const TABLE_INDEX table_cur = ivolv_list[cur].table_index;
      bool curOnLower = ivoldual_table.OnLowerIsosurface(table_cur, ivolv_cur);
      bool curOnUpper = ivoldual_table.OnUpperIsosurface(table_cur, ivolv_cur);
      bool isOnSurface = curOnLower || curOnUpper;

      if (isOnSurface == skipSurfaceVert) continue;

      // Store sum of neighbor coordinates.
      COORD_TYPE neigh_sum[d]; 
      IJK::set_coord(d, 0.0, neigh_sum);
      bool flag_moving = false;
      int adj_count = 0;

      // Loop over adjacent vertices of the current vertex
      for (int  k = 0; k < vertex_adjacency_list.NumAdjacent(cur); k++) {

        // Neighbor node coordinates
        int adj = vertex_adjacency_list.AdjacentVertex(cur, k);
        COORD_TYPE *neigh_coord = vcoord + adj*d;

        // Check if neighbor node is on isosurface.
        const int ivolv_adj = ivolv_list[adj].patch_index;
        const TABLE_INDEX table_adj = ivolv_list[adj].table_index;
        bool adjOnLower = ivoldual_table.OnLowerIsosurface(table_adj, ivolv_adj);
        bool adjOnUpper = ivoldual_table.OnUpperIsosurface(table_adj, ivolv_adj);

        // Skip if a vertex and its adjacent vertex are not on the same surface.
        if ((curOnLower && !adjOnLower) || (curOnUpper && !adjOnUpper))
            continue;

        IJK::add_coord(d, neigh_sum, neigh_coord, neigh_sum);
        IJK::compute_distance(d, cur_coord, neigh_coord, dist);

        // Check if minimum distance is valid.
        if (dist < laplacian_smooth_limit) {
          flag_moving = true;
        }

        adj_count++;
      }

      // Update current node coordinate.
      if (flag_moving) {
        IJK::divide_coord(d, adj_count, neigh_sum, neigh_sum);
        IJK::copy_coord(d, neigh_sum, cur_coord);
      }
    }
  }
}

void IVOLDUAL::laplacian_smooth_jacobian
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord, 
 float jacobian_limit, 
 int iteration)
 {
  for (int it = 0; it < iteration; it++) {
    std::vector<int> neg_jacob_list;
    // Find all vertices with negative Jacobian.
    for (int ihex = 0; ihex < ivolpoly_vert.size()/8; ihex++) {
      for (int i = 0; i < 8; i++) {
        // Compute Jacobian at current vertex
        COORD_TYPE jacob;        
        compute_hexahedron_normalized_Jacobian_determinant
          (ivolpoly_vert, ihex, vertex_coord, i, jacob);

        if (jacob < jacobian_limit) {
          neg_jacob_list.push_back(ivolpoly_vert[ihex * 8 + i]);
        }
      }
    }
    laplacian_smooth_jacobian
      (ivolpoly_vert, ivoldual_table, vertex_adjacency_list,  
       vertex_poly_incidence, ivolv_list, vertex_coord, neg_jacob_list);
  }
}


void IVOLDUAL::laplacian_smooth_jacobian
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord, 
 const std::vector<int> & neg_jacobian_list)
{
	const int DIM3(3);
  COORD_TYPE * vcoord = &(vertex_coord.front());

	for (int cur : neg_jacobian_list) {

    // Current node coordinates.
    COORD_TYPE *cur_coord = vcoord + cur * DIM3;

    // Check if current node is on isosurface.
    const int ivolv_cur = ivolv_list[cur].patch_index;
    const TABLE_INDEX table_cur = ivolv_list[cur].table_index;
    const VERTEX_INDEX cube_cur = ivolv_list[cur].cube_index;
    bool curOnLower = ivoldual_table.OnLowerIsosurface(table_cur, ivolv_cur);
    bool curOnUpper = ivoldual_table.OnUpperIsosurface(table_cur, ivolv_cur);

    // Loop over adjacent vertices of the current vertex
    for (int j = 0; j < vertex_adjacency_list.NumAdjacent(cur); j++) {

      // Neighbor node coordinates
      int adj = vertex_adjacency_list.AdjacentVertex(cur, j);
      COORD_TYPE *neigh_coord = vcoord + adj * DIM3;
      
      // Check if neighbor node is on isosurface.
      const int ivolv_adj = ivolv_list[adj].patch_index;
      const TABLE_INDEX table_adj = ivolv_list[adj].table_index;
      const VERTEX_INDEX cube_adj = ivolv_list[adj].cube_index;
      bool adjOnLower = ivoldual_table.OnLowerIsosurface(table_adj, ivolv_adj);
      bool adjOnUpper = ivoldual_table.OnUpperIsosurface(table_adj, ivolv_adj);

      if (cube_cur == cube_adj) {
      	gradient_move_vertex
      	(ivolpoly_vert, ivoldual_table, vertex_adjacency_list, vertex_poly_incidence, 
         ivolv_list, vertex_coord,neigh_coord, adj, adjOnLower, adjOnUpper);
	      gradient_move_vertex
      	(ivolpoly_vert, ivoldual_table, vertex_adjacency_list, vertex_poly_incidence, 
         ivolv_list, vertex_coord,cur_coord, cur, curOnLower, curOnUpper);
      }
    }
  }
}

void IVOLDUAL::gradient_smooth_jacobian
(std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 DUAL_IVOLVERT_ARRAY & ivolv_list,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 COORD_ARRAY & vertex_coord, 
 float jacobian_limit, 
 int iteration)
{
  const int NUM_VERT_PER_HEX(8); 
  COORD_TYPE * vcoord = &(vertex_coord.front());

  for (int it = 0; it < iteration; it++) {
    std::vector<int> neg_jacobian_list;
    std::vector<std::vector<int>> facet_list, edge_list;
    std::unordered_map<int, COORD_TYPE> neg_jacob_value;

    for (int ihex = 0; ihex < ivolpoly_vert.size()/8; ihex++) {

      std::vector<int> internal_vert;

      bool small_jacob_hex = false;

      // Check is current hex has Jacobian below threshold
      for (int i = 0; i < NUM_VERT_PER_HEX; i++) {
        // Compute Jacobian at current vertex
        COORD_TYPE jacob;        
        compute_hexahedron_normalized_Jacobian_determinant
          (ivolpoly_vert, ihex, vertex_coord, i, jacob);

        if (jacob < jacobian_limit) {
          small_jacob_hex = true;
          break;
        }
      }

      if (small_jacob_hex == true) {
        for (int i = 0; i < NUM_VERT_PER_HEX; i++) {

          int ivert = ivolpoly_vert[ihex * 8 + i];
          
          // Check if current node is on isosurface.
          const int ivolv_cur = ivolv_list[ivert].patch_index;
          const TABLE_INDEX table_cur = ivolv_list[ivert].table_index;
          bool curOnLower = ivoldual_table.OnLowerIsosurface(table_cur, ivolv_cur);
          bool curOnUpper = ivoldual_table.OnUpperIsosurface(table_cur, ivolv_cur);
          if (curOnLower || curOnUpper) continue;
          
          internal_vert.push_back(ivert);

          COORD_TYPE jacob;        
          compute_hexahedron_normalized_Jacobian_determinant
            (ivolpoly_vert, ihex, vertex_coord, i, jacob);

          if (jacob < jacobian_limit) {
            // A vertex is not in the negative jacobian list
            if (neg_jacob_value.count(ivert) == 0) {
              neg_jacobian_list.push_back(ivert);
              neg_jacob_value[ivert] = jacob;
            }
            else if (jacob < neg_jacob_value[ivert]) {
              neg_jacob_value[ivert] = jacob;
            } 
          }   
        }

        facet_list.push_back(internal_vert);

        for (int i = 0; i < internal_vert.size() - 1; i++) {
          for (int j = i + 1; j < internal_vert.size(); j++) {
            int iv0 = internal_vert[i], iv1 = internal_vert[j];
            if (vertex_adjacency_list.IsAdjacent(iv0, iv1)) {
              edge_list.push_back({iv0, iv1});
            }
          }
        }
      }
    }

    // Smoothing vertices.
    gradient_smooth_jacobian
    (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, vertex_poly_incidence, 
     ivolv_list, vertex_coord, neg_jacobian_list, neg_jacob_value, it);
    
    if (it == iteration - 1) {
      // Smoothing flat hex facet
      expand_flat_hex
      (ivolpoly_vert, ivoldual_table, vertex_adjacency_list,  
       vertex_poly_incidence, ivolv_list, vertex_coord, facet_list); 
      // Smoothing flat hex edge
      expand_flat_hex
      (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, 
       vertex_poly_incidence, ivolv_list, vertex_coord, edge_list); 
    }
  }

  // Split hex to fix indented small Jacobian. Only works for indented vertex.
  // split_hex
  // (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, ivolv_list,
  //  ivolpoly_info, vertex_coord, 0.1);
}

void IVOLDUAL::gradient_smooth_jacobian
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord, 
 const std::vector<int> & neg_jacobian_list,
 std::unordered_map<int, COORD_TYPE> & neg_jacob_value,
 int iter)
{
  for (int i = 0; i < neg_jacobian_list.size(); i++) {
    int ivert = neg_jacobian_list[i];
    COORD_TYPE cur_min_jacob = neg_jacob_value[ivert];

    move_vertex_all_direction
    (ivolpoly_vert, vertex_adjacency_list, vertex_poly_incidence, 
     vertex_coord, cur_min_jacob, ivert, iter);
  } 
}

void IVOLDUAL::expand_flat_hex
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord, 
 const std::vector<std::vector<int>> & flat_hex)
{
  const int DIM3(3);
  COORD_TYPE * vcoord = &(vertex_coord.front());

  // Loop over hexahedra facets with small Jacobian.
  for (int ifacet = 0; ifacet < flat_hex.size(); ifacet++) {

    float move_dist = 0.0;
    std::vector<int> dir(3);
    float pre_min_at_facet = 1.0, pre_min_around_facet = 1.0;

    // Find min Jacobian at/around a facet/edge
    for (int j = 0; j < flat_hex[ifacet].size(); j++) {
      int ivert = flat_hex[ifacet][j];
      float min_jacob_at_cur, min_jacob_around_cur;

      min_jacob_around_vertex
      (ivolpoly_vert, vertex_poly_incidence, vertex_coord,
       ivert, min_jacob_at_cur, min_jacob_around_cur);

      pre_min_at_facet = std::min(pre_min_at_facet, min_jacob_at_cur);
      pre_min_around_facet = std::min(pre_min_around_facet, min_jacob_around_cur);
    } 

    find_optimal_jacobian_point
    (ivolpoly_vert, vertex_poly_incidence, vertex_coord, flat_hex, 
     pre_min_at_facet, pre_min_around_facet, ifacet, move_dist, dir);

    // Move to optiminal position.
    for (int j = 0; j < flat_hex[ifacet].size(); j++) {
      int ivert = flat_hex[ifacet][j];
      COORD_TYPE *cur_coord = vcoord + ivert * DIM3;
      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] += move_dist * dir[d];
      }
    }
  } 
}

void IVOLDUAL::find_optimal_jacobian_point
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 COORD_ARRAY & vertex_coord, 
 const std::vector<std::vector<int>> & flat_hex, 
 float pre_min_at_facet, float pre_min_around_facet, 
 int ifacet, 
 float & move_dist, std::vector<int> & dir)
{
  const int DIM3(3);
  const float step_base(0.05);
  float jacob_limit = 0.2;
  COORD_TYPE * vcoord = &(vertex_coord.front());

  for (int k = 1; k * step_base < 0.4; k++) {
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
        for (int iz = -1; iz <= 1; iz++) {

          if (ix == 0 && iy == 0 && iz == 0) continue;

          std::vector<int> dir_temp(3);
          dir_temp[0] = ix;
          dir_temp[1] = iy;
          dir_temp[2] = iz;

          // Move in dir_temp
          for (int j = 0; j < flat_hex[ifacet].size(); j++) {
            int ivert = flat_hex[ifacet][j];
            COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

            for (int d = 0; d < DIM3; d++) {
              cur_coord[d] += step_base * k * dir_temp[d];
            }
          }

          // Check Jacobian.
          float min_at_facet = 1.0, min_around_facet = 1.0;
          for (int j = 0; j < flat_hex[ifacet].size(); j++) {
            int ivert = flat_hex[ifacet][j];
            float min_jacob_at_cur, min_jacob_around_cur;

            min_jacob_around_vertex
            (ivolpoly_vert, vertex_poly_incidence, vertex_coord,
             ivert, min_jacob_at_cur, min_jacob_around_cur);

            min_at_facet = std::min(min_at_facet, min_jacob_at_cur);
            min_around_facet = std::min(min_around_facet, min_jacob_around_cur);
          }

          // If Jacobian improved, store the distance and direction.
          if (min_at_facet > std::min(pre_min_at_facet, jacob_limit) && 
              min_around_facet >= std::min(pre_min_around_facet, jacob_limit) ) {
            pre_min_at_facet = min_at_facet;
            pre_min_around_facet = min_around_facet;
            move_dist = k * step_base;
            dir = dir_temp;
          }

          // Move back to original positions
          for (int j = 0; j < flat_hex[ifacet].size(); j++) {
            int ivert = flat_hex[ifacet][j];
            COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

            for (int d = 0; d < DIM3; d++) {
              cur_coord[d] -= step_base * k * dir_temp[d];
            }
          }
        }
      }
    }
  }
}

void IVOLDUAL::expand_flat_hex_normal_direction
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord, 
 const std::vector<int> & flat_hex)
{
  const int DIM3(3);
  const float step_base(0.05);
  const int NUM_VERT_PER_HEXAHEDRON(8); 
  float jacob_limit = 0.2;
  COORD_TYPE * vcoord = &(vertex_coord.front());

  for (int i = 0; i < flat_hex.size(); i++) {
    int ihex = flat_hex[i];

    std::vector<int> surface_vert, internal_vert;

    // Find surface and internal vertices
    for (int i = 0; i < 8; i++) {

      int ivert = ivolpoly_vert[ihex * 8 + i];
      
      // Check if current node is on isosurface.
      const int ivolv_cur = ivolv_list[ivert].patch_index;
      const TABLE_INDEX table_cur = ivolv_list[ivert].table_index;
      bool curOnLower = ivoldual_table.OnLowerIsosurface(table_cur, ivolv_cur);
      bool curOnUpper = ivoldual_table.OnUpperIsosurface(table_cur, ivolv_cur);
      if (curOnLower || curOnUpper) {
        surface_vert.push_back(ivert);
        continue;
      }

      internal_vert.push_back(ivert);  
    }

    if (surface_vert.size() != 4) continue;

    // Find normal direction of surface
    COORD_TYPE normal_dir[3];
    surface_normal_direction
     (vertex_coord, surface_vert, normal_dir);

    move_vertex_normal_direction
     (ivolpoly_vert, vertex_poly_incidence, 
      vertex_coord, internal_vert, normal_dir);

    for (int j = 0; j < internal_vert.size() - 1; j++) {
      std::vector<int> internal_edge;
      for (int k = j + 1; k < internal_vert.size(); k++) {
        int iv0 = internal_vert[j], iv1 = internal_vert[k];
        if (vertex_adjacency_list.IsAdjacent(iv0, iv1)) {
          internal_edge.push_back(iv0);
          internal_edge.push_back(iv1);
          move_vertex_normal_direction
           (ivolpoly_vert, vertex_poly_incidence, 
            vertex_coord, internal_edge, normal_dir);
        }
      }
    }
  }
}

void IVOLDUAL::move_vertex_along_edge
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 COORD_ARRAY & vertex_coord,
 COORD_TYPE cur_min_jacob,
 int ivert)
{
  const int DIM3(3);
  const float step_base(0.1);
  const int NUM_VERT_PER_HEXAHEDRON(8); 

  COORD_TYPE * vcoord = &(vertex_coord.front());
  COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

  COORD_TYPE target[DIM3];
  COORD_TYPE pre_jacobian = cur_min_jacob;

  // Optimal position to be moved to.
  for (int d = 0; d < DIM3; d++) {
    target[d] = cur_coord[d];
  }

  // Backup current vertex coordinates
  COORD_TYPE cur_coord_backup[DIM3];
  for (int d = 0; d < DIM3; d++) {
    cur_coord_backup[d] = cur_coord[d];
  }

  // Loop over all neighbor vertices to find the best moving direction
  for (int k = 0; k < vertex_adjacency_list.NumAdjacent(ivert); k++) {
    int adj = vertex_adjacency_list.AdjacentVertex(ivert, k);
    COORD_TYPE *neigh_coord = vcoord + adj*DIM3;

    // For different step length
    for (int i = -4; step_base * i < 0.7; i++) {
      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] = (1.0-step_base*i)*cur_coord[d] + step_base*i*neigh_coord[d];
      }

      // Find the min Jacobian of all neighbor hex
      COORD_TYPE min_jacobian = 1.0;
      for (int ipoly = 0; ipoly < vertex_poly_incidence.NumIncidentPoly(ivert); ipoly++) {
        const int jhex = vertex_poly_incidence.IncidentPoly(ivert, ipoly);

        for (int i = 0; i < 8; i++) {
          // Compute Jacobian at current vertex
          COORD_TYPE jacob;        
          compute_hexahedron_normalized_Jacobian_determinant
            (ivolpoly_vert, jhex, vertex_coord, i, jacob);

          min_jacobian = std::min(jacob, min_jacobian);
        }
      }

      if (min_jacobian > pre_jacobian) {
        pre_jacobian = min_jacobian;
        for (int d = 0; d < DIM3; d++) {
          target[d] = cur_coord[d];
        }
      }
      // Restore neighbor coordinate to backup coord
      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] = cur_coord_backup[d];
      }
    }
  }

  // Move the coordinate to the optimal position.
  for (int d = 0; d < DIM3; d++) {
    cur_coord[d] = target[d];
  }
}


void IVOLDUAL::move_vertex_all_direction
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 COORD_ARRAY & vertex_coord,
 COORD_TYPE pre_jacob,
 int ivert, int iter)
{
  COORD_TYPE pre_jacob_around = pre_jacob;

  const int DIM3(3);
  float step_base(0.05);
  const int NUM_VERT_PER_HEXAHEDRON(8); 

  COORD_TYPE * vcoord = &(vertex_coord.front());
  COORD_TYPE *cur_coord = vcoord + ivert * DIM3;
  COORD_TYPE target[DIM3];

  // Optimal position to be moved to.
  for (int d = 0; d < DIM3; d++) {
    target[d] = cur_coord[d];
  }

  // Backup current vertex coordinates
  COORD_TYPE cur_coord_backup[DIM3];
  for (int d = 0; d < DIM3; d++) {
    cur_coord_backup[d] = cur_coord[d];
  }

  // For different step length
  for (int i = 1; step_base * i <= 0.5; i++) {

    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
        for (int iz = -1; iz <= 1; iz++) {

          if (ix == 0 && iy == 0 && iz == 0) continue;

          cur_coord[0] += step_base*i*ix;
          cur_coord[1] += step_base*i*iy;
          cur_coord[2] += step_base*i*iz;

          COORD_TYPE min_jacob_at_cur, min_jacob_around_cur;
          min_jacob_around_vertex
            (ivolpoly_vert, vertex_poly_incidence, vertex_coord,
             ivert, min_jacob_at_cur, min_jacob_around_cur);

          if (min_jacob_at_cur > pre_jacob &&
              min_jacob_around_cur >= pre_jacob_around) {
            pre_jacob = min_jacob_at_cur;
            for (int d = 0; d < DIM3; d++) {
              target[d] = cur_coord[d];
            }
          }

          // Restore neighbor coordinate to backup coord
          for (int d = 0; d < DIM3; d++) {
            cur_coord[d] = cur_coord_backup[d];
          }
        }
      }
    }
  }
  

  // Move the coordinate to the optimal position.
  for (int d = 0; d < DIM3; d++) {
    cur_coord[d] = target[d];
  }

}

void IVOLDUAL::gradient_move_vertex
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 const DUAL_IVOLVERT_ARRAY & ivolv_list,
 COORD_ARRAY & vertex_coord,
 COORD_TYPE *ver_coord, int ver_index,
 bool flag_onLower,  bool flag_onUpper)
{
	const int DIM3(3);
  const float step_base(0.1);
  const int NUM_VERT_PER_HEXAHEDRON(8); 
  COORD_TYPE * vcoord = &(vertex_coord.front());

  COORD_TYPE target[DIM3];
  float pre_jacobian = -1.0;

  // Optimal position to be moved to.
  for (int d = 0; d < DIM3; d++) {
    target[d] = ver_coord[d];
  }

  // Find the direction along which Jacobian changes most.
  for (int k = 0; k < vertex_adjacency_list.NumAdjacent(ver_index); k++) {

    int adj = vertex_adjacency_list.AdjacentVertex(ver_index, k);
    COORD_TYPE *neigh_coord = vcoord + adj*DIM3;

    const int ivolv_adj = ivolv_list[adj].patch_index;
    const TABLE_INDEX table_adj = ivolv_list[adj].table_index;
    const VERTEX_INDEX cube_adj = ivolv_list[adj].cube_index;
    bool adjOnLower = ivoldual_table.OnLowerIsosurface(table_adj, ivolv_adj);
    bool adjOnUpper = ivoldual_table.OnUpperIsosurface(table_adj, ivolv_adj);

    // Skip if two vertices are not on same surface.
    if ((flag_onLower && !adjOnLower) || (flag_onUpper && !adjOnUpper))
      continue;

    // copy adj coord to a temp_adj_coord
    COORD_TYPE neigh_temp[DIM3];
    
    // Backup Coord
    for (int d = 0; d < DIM3; d++) {
      neigh_temp[d] = ver_coord[d];
    }

    for (int d = 0; d < DIM3; d++) {
      ver_coord[d] = (1.0-step_base)*ver_coord[d] + step_base*neigh_coord[d];
    }

    COORD_TYPE min_jacobian = 1.0;
    for (int ipoly = 0; ipoly < vertex_poly_incidence.NumIncidentPoly(ver_index); ipoly++) {
      const int ihex = vertex_poly_incidence.IncidentPoly(ver_index, ipoly);

      for (int i = 0; i < 8; i++) {
        // Compute Jacobian at current vertex
        COORD_TYPE jacob;        
        compute_hexahedron_normalized_Jacobian_determinant
          (ivolpoly_vert, ihex, vertex_coord, i, jacob);

        min_jacobian = std::min(jacob, min_jacobian);
      }
    }

    if (min_jacobian > pre_jacobian) {
      pre_jacobian = min_jacobian;
      for (int d = 0; d < DIM3; d++) {
        target[d] = neigh_coord[d];
      }
    }
    // Restore neighbor coordinate to backup coord
    for (int d = 0; d < DIM3; d++) {
      ver_coord[d] = neigh_temp[d];
    }
  }

  // Move the vertex along max gradient direction.
  COORD_TYPE step[DIM3], optimal[DIM3];
  pre_jacobian = -1.0;

  for (int d = 0; d < DIM3; d++) {
  	step[d] = (target[d] - ver_coord[d]) * step_base;
  }
  for (int i = 1; step_base * i <= 0.6; i++) {
		for (int d = 0; d < DIM3; d++) {
			ver_coord[d] += step[d];
		}
		COORD_TYPE min_jacobian = 1.0;
    for (int ipoly = 0; ipoly < vertex_poly_incidence.NumIncidentPoly(ver_index); ipoly++) {
      const int ihex = vertex_poly_incidence.IncidentPoly(ver_index, ipoly);

      for (int i = 0; i < 8; i++) {
        // Compute Jacobian at current vertex
        COORD_TYPE jacob;        
        compute_hexahedron_normalized_Jacobian_determinant
          (ivolpoly_vert, ihex, vertex_coord, i, jacob);
        min_jacobian = std::min(jacob, min_jacobian);
      }
    }
    if (min_jacobian > pre_jacobian) {
      pre_jacobian = min_jacobian;
      for (int d = 0; d < DIM3; d++) {
        optimal[d] = ver_coord[d];
      }
    }
  }

  // Move the coordinate to the optimal position.
  for (int d = 0; d < DIM3; d++) {
    ver_coord[d] = optimal[d];
  }
}

void IVOLDUAL::move_vertex_normal_direction
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 COORD_ARRAY & vertex_coord,
 const std::vector<int> & internal_vert,
 COORD_TYPE normal_dir[])
 {
  const int DIM3(3);
  const float step_base(0.05);
  const int NUM_VERT_PER_HEXAHEDRON(8); 
  float jacob_limit = 0.2;
  COORD_TYPE * vcoord = &(vertex_coord.front());

  float pre_min_at_facet = 1.0, pre_min_around_facet = 1.0;

  for (int j = 0; j < internal_vert.size(); j++) {
    int ivert = internal_vert[j];
    float min_jacob_at_cur, min_jacob_around_cur;

    min_jacob_around_vertex
    (ivolpoly_vert, vertex_poly_incidence, vertex_coord,
     ivert, min_jacob_at_cur, min_jacob_around_cur);

    pre_min_at_facet = std::min(pre_min_at_facet, min_jacob_at_cur);
    pre_min_around_facet = std::min(pre_min_around_facet, min_jacob_around_cur);
  }

  float move_dist = 0.0;

  for (int k = 1; k * step_base < 0.6; k++) {

    // Move in dir_temp
    for (int j = 0; j < internal_vert.size(); j++) {
      int ivert = internal_vert[j];
      COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] += step_base * k * normal_dir[d];
      }
    }

    float min_at_facet = 1.0, min_around_facet = 1.0;
    for (int j = 0; j < internal_vert.size(); j++) {
      int ivert = internal_vert[j];
      float min_jacob_at_cur, min_jacob_around_cur;

      min_jacob_around_vertex
      (ivolpoly_vert, vertex_poly_incidence, vertex_coord,
       ivert, min_jacob_at_cur, min_jacob_around_cur);

      min_at_facet = std::min(min_at_facet, min_jacob_at_cur);
      min_around_facet = std::min(min_around_facet, min_jacob_around_cur);
    }

    if (min_at_facet > std::min(pre_min_at_facet, jacob_limit) && 
        min_around_facet >= std::min(pre_min_around_facet, jacob_limit) ) {
      pre_min_at_facet = min_at_facet;
      pre_min_around_facet = min_around_facet;
      move_dist = k * step_base;
    }

    // Move back to original positions
    for (int j = 0; j < internal_vert.size(); j++) {
      int ivert = internal_vert[j];
      COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] -= step_base * k * normal_dir[d];
      }
    }
  }

  for (int k = -1; k * step_base > -0.6; k--) {

    // Move in dir_temp
    for (int j = 0; j < internal_vert.size(); j++) {
      int ivert = internal_vert[j];
      COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] += step_base * k * normal_dir[d];
      }
    }

    float min_at_facet = 1.0, min_around_facet = 1.0;
    for (int j = 0; j < internal_vert.size(); j++) {
      int ivert = internal_vert[j];
      float min_jacob_at_cur, min_jacob_around_cur;

      min_jacob_around_vertex
      (ivolpoly_vert, vertex_poly_incidence, vertex_coord,
       ivert, min_jacob_at_cur, min_jacob_around_cur);

      min_at_facet = std::min(min_at_facet, min_jacob_at_cur);
      min_around_facet = std::min(min_around_facet, min_jacob_around_cur);
    }

    if (min_at_facet > std::min(pre_min_at_facet, jacob_limit) && 
        min_around_facet >= std::min(pre_min_around_facet, jacob_limit) ) {
      pre_min_at_facet = min_at_facet;
      pre_min_around_facet = min_around_facet;
      move_dist = k * step_base;
    }

    // Move back to original positions
    for (int j = 0; j < internal_vert.size(); j++) {
      int ivert = internal_vert[j];
      COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] -= step_base * k * normal_dir[d];
      }
    }
  }

  // Find an optimization point to move
  if (move_dist > 0.0) {
    for (int j = 0; j < internal_vert.size(); j++) {
      int ivert = internal_vert[j];
      COORD_TYPE *cur_coord = vcoord + ivert * DIM3;

      for (int d = 0; d < DIM3; d++) {
        cur_coord[d] += move_dist * normal_dir[d];
      }
    }
  }
}


void IVOLDUAL::min_jacob_around_vertex
(const std::vector<VERTEX_INDEX> & ivolpoly_vert,
 IJK::VERTEX_POLY_INCIDENCE<int,int> & vertex_poly_incidence,
 COORD_ARRAY & vertex_coord,
 int & ivert,
 COORD_TYPE & min_jacob_at_cur,
 COORD_TYPE & min_jacob_around_cur)
{
  // Find min Jacobian around the vertex
  min_jacob_at_cur = 1.0; 
  min_jacob_around_cur = 1.0;

  for (int ipoly = 0; ipoly < vertex_poly_incidence.NumIncidentPoly(ivert); ipoly++) {
    const int ihex = vertex_poly_incidence.IncidentPoly(ivert, ipoly);

    for (int j = 0; j < 8; j++) {
      // Compute Jacobian at current vertex
      COORD_TYPE jacob;        
      compute_hexahedron_normalized_Jacobian_determinant
        (ivolpoly_vert, ihex, vertex_coord, j, jacob);

      min_jacob_around_cur = std::min(jacob, min_jacob_around_cur);
      if (ivert == ivolpoly_vert[ihex * 8 + j]) {
        min_jacob_at_cur = std::min(jacob, min_jacob_at_cur);
      }
    }
  }
}

void IVOLDUAL::surface_normal_direction
(COORD_ARRAY & vertex_coord,
 const std::vector<int> & surface_vert,
 COORD_TYPE normal_dir[])
{
  const int DIM3(3);
  COORD_TYPE * vcoord = &(vertex_coord.front());

  // Three vertices on the surface
  COORD_TYPE * v0 = vcoord + surface_vert[0]*DIM3;
  COORD_TYPE * v1 = vcoord + surface_vert[1]*DIM3;
  COORD_TYPE * v2 = vcoord + surface_vert[2]*DIM3;

  // The edges formed by the three vertices
  COORD_TYPE e1[3], e2[3], dir[3];

  IJK::subtract_coord(DIM3, v1, v0, e1);
  IJK::subtract_coord(DIM3, v2, v0, e2);

  // Compute normal direction
  IJK::compute_cross_product_3D(v0, v1, dir);

  // Normalize
  float magnitude;
  bool flag_zero;
  normalize_vector
   (DIM3, dir, 0.0, normal_dir, magnitude, flag_zero);
}
