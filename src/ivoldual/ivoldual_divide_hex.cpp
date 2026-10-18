/// @file ivoldual_divide_hex.cpp
/// @brief Reposition vertices for mesh quality optimization.
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

#include <unordered_map>
#include "ivoldual_compute.h"
#include "ivoldual_divide_hex.h"

#include "ijktriangulate.txx"

using namespace IJK;
using namespace IVOLDUAL;

void IVOLDUAL::collapse_hex
(std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 COORD_ARRAY & vertex_coord, 
 COORD_TYPE jacobian_limit)
 {
  const int NUM_VERT_PER_HEX(8);
  std::vector<std::pair<VERTEX_INDEX, VERTEX_INDEX>> vertex_subdivide_list;

  // Loop over polytopeS to find vertex with negative Jacobian and indentation.
  for (int ihex = 0; ihex < ivolpoly_vert.size() / NUM_VERT_PER_HEX; ihex++) {
    for (int icorner = 0; icorner < NUM_VERT_PER_HEX; icorner++) {
    	
    	int cur = ivolpoly_vert[ihex*8+icorner];
      // Compute Jacobian at current vertex
      COORD_TYPE jacob;        
      compute_hexahedron_normalized_Jacobian_determinant
        (ivolpoly_vert, ihex, vertex_coord, icorner, jacob);
      if (jacob < jacobian_limit && ivolv_list[cur].num_incident_hex == 4 && 
      	  ivolv_list[cur].num_incident_iso_quad == 0)
      {
      	// Push back the opposite vertex in the cube to the negative Jacobian vertex.
      	int index_indent = ivolpoly_vert[ihex*8+icorner];
      	int index_indent_oppo = ivolpoly_vert[ihex*8 + 7 - icorner];

      	if (ivolv_list[index_indent_oppo].num_incident_iso_quad != 0 && 
      		  ivolv_list[index_indent_oppo].num_incident_hex == 1) {
	      	ivolv_list[index_indent].map_to = index_indent_oppo;
	      	ivolpoly_info[ihex].flag_subdivide_hex = true;
	      }
      }
    }
  }

  // Remove deleted hexahedra.
  std::vector<VERTEX_INDEX> ivolpoly_vert_new;
  for (int i = 0; i < ivolpoly_info.size(); i++) {
  	if (!ivolpoly_info[i].flag_subdivide_hex) {
			for (int j = 0; j < 8; j++) {
				int cur = ivolpoly_vert[i*8 + j];
				ivolpoly_vert_new.push_back(ivolv_list[cur].map_to);
			} 		
  	}
  }
  ivolpoly_vert = move(ivolpoly_vert_new);
}

void IVOLDUAL::split_hex
(std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 COORD_ARRAY & vertex_coord, 
 COORD_TYPE jacobian_limit)
 {
  const int NUM_VERT_PER_HEX(8);
  std::vector<std::pair<VERTEX_INDEX, VERTEX_INDEX>> vertex_subdivide_list;

  // Loop over every polytope to find vertex with negative Jacobian and indentation feature.
  for (int ihex = 0; ihex < ivolpoly_vert.size() / NUM_VERT_PER_HEX; ihex++) {
    for (int icorner = 0; icorner < NUM_VERT_PER_HEX; icorner++) {
    	int cur = ivolpoly_vert[ihex*8+icorner];

      // Compute Jacobian at current vertex
      COORD_TYPE jacob;        
      compute_hexahedron_normalized_Jacobian_determinant
        (ivolpoly_vert, ihex, vertex_coord, icorner, jacob);

      if (jacob < jacobian_limit)
      {
      	// Push back the opposite vertex in the cube to the negative Jacobian vertex.
      	vertex_subdivide_list.push_back
      		(std::make_pair(ivolpoly_vert[ihex*8+(7-icorner)], ivolpoly_vert[ihex*8+icorner]));
      }
    }
  }

  // Subdivide each polytope into four hexahedra.
  subdivide_hex_to_four
    (ivolpoly_vert, ivoldual_table, vertex_adjacency_list, ivolv_list, 
     ivolpoly_info, vertex_subdivide_list, vertex_coord);

  // Remove deleted hexahedra.
  std::vector<VERTEX_INDEX> ivolpoly_vert_new;
  for (int i = 0; i < ivolpoly_info.size(); i++) {
  	if (!ivolpoly_info[i].flag_subdivide_hex) {
  		for (int j = 0; j < 8; j++) {
  			ivolpoly_vert_new.push_back(ivolpoly_vert[i*8 + j]);
  		} 		
  	}
  }
  ivolpoly_vert = move(ivolpoly_vert_new);
}

void IVOLDUAL::subdivide_hex_to_four
(std::vector<VERTEX_INDEX> & ivolpoly_vert,
 const IVOLDUAL_CUBE_TABLE & ivoldual_table,
 IVOL_VERTEX_ADJACENCY_LIST & vertex_adjacency_list,
 DUAL_IVOLVERT_ARRAY & ivolv_list,
 IVOLDUAL_POLY_INFO_ARRAY & ivolpoly_info,
 std::vector<std::pair<VERTEX_INDEX, VERTEX_INDEX>> &
 	vertex_subdivide_list,
 COORD_ARRAY & vertex_coord)
{
	VERTEX_INDEX corner[8];
  VERTEX_INDEX iv[8];
  VERTEX_INDEX iw[8];
  std::unordered_map<int, int> new_vertex;
  const int DIM3(3);
  const int NUM_VERT_PER_HEX(8);
  IJK::CUBE_FACE_INFO<int, int, int> cube(DIM3);
  std::unordered_map<VERTEX_INDEX, int> forbiden;

  // Polytopes dual to vertex.
  IJK::POLYMESH_DATA<VERTEX_INDEX,int, 
    IJK::HEX_TRIANGULATION_INFO<char,char>> hex_data;
  hex_data.AddPolytopes(ivolpoly_vert, NUM_VERT_PER_HEX);
  IJK::VERTEX_POLY_INCIDENCE<int,int> vertex_poly_incidence(hex_data);

  // Loop over all vertices around which the hexhedra needs subdivid.
  for (auto ipair : vertex_subdivide_list) {
  	VERTEX_INDEX ivertex = ipair.first, indent_vertex = ipair.second;
  	// Skip vertex in polytopes which are already subdivided
  	if (forbiden[ivertex] == 1) continue;
  	// Loop over all polytopes incident on the vertex.
  	for (int ipoly = 0; ipoly < vertex_poly_incidence.NumIncidentPoly(ivertex); ipoly++) {
			const int ihex = vertex_poly_incidence.IncidentPoly(ivertex, ipoly);
			// Flag subdivided hex.
			ivolpoly_info[ihex].flag_subdivide_hex = true;
			ivolpoly_info.resize(ivolpoly_info.size() + 4);

			// Forbid subdividing other polytopes around this vertex
			for (int i = 0; i < 8; i++) {
				forbiden[ivolpoly_vert[ihex*8+i]] = 1;
			}

			for (int i = 0; i < 8; i++) {
				if (ivolpoly_vert[ihex*8+i] == ivertex) {
					corner[0] = i;
					corner[1] = cube.VertexNeighbor(corner[0], 0);
					break;
				}
			}

			// Get all 8 corner indices.			
			for (int i = 0; i < 2; i++) corner[i+2] = cube.VertexNeighbor(corner[i], 1);
			for (int i = 0; i < 4; i++) corner[i+4] = cube.VertexNeighbor(corner[i], 2);
			// Get all 8 vertex indices.
			for (int i = 0; i < 8; i++) {
				iv[i] = ivolpoly_vert[ihex*8+corner[i]];
			}


			// Generate new vertices for subdivision.
			for (int i = 1; i < 8; i++) {
				// Create a new vertex
				if (new_vertex.count(iv[i]) == 0) {					
					iw[i] = vertex_coord.size()/DIM3;
					vertex_coord.resize(vertex_coord.size() + DIM3);
					for (int d = 0; d < DIM3; d++) {
						vertex_coord[iw[i]*DIM3+d] = 0.5 * (vertex_coord[iv[0] * DIM3 + d] + 
															vertex_coord[iv[i] * DIM3 + d]);
					}
					new_vertex[iv[i]] = iw[i];
					ivolv_list.push_back(ivolv_list.back());
				}
				// Get a existing vertex.
				else {
					iw[i] = new_vertex[iv[i]];
				}
			}

			// Special scale ratio for indentation hexahedra.
			if (iv[7] == indent_vertex) {
				float dist0, dist1;
				COORD_TYPE * vcoord0 = &(vertex_coord.front()) + iv[0] * DIM3;
				COORD_TYPE * vcoord1 = &(vertex_coord.front()) + iv[7] * DIM3;
				COORD_TYPE * vcoord_iw7 = &(vertex_coord.front()) + iw[7] * DIM3;

				for (int d = 0; d < DIM3; d++) {
					vertex_coord[iw[7]*DIM3+d] = 0.25*vcoord0[d] + 0.75*vcoord1[d];
				}


				IJK::compute_distance(DIM3, vcoord0, vcoord_iw7, dist0);
				for (int i = 1; i < 7; i++) {
					vcoord1 = &(vertex_coord.front()) + iv[i] * DIM3;
					IJK::compute_distance(DIM3, vcoord0, vcoord1, dist1);
					for (int d = 0; d < DIM3; d++) {
						vertex_coord[iw[i]*DIM3+d] = vcoord0[d] + dist0/dist1*(vcoord1[d]-vcoord0[d]);
					}
				}
			}
			// Subdivide hexahedron.
			generate_children_polytope(ivolpoly_vert, iv, iw, corner);
    }
  }
}

void IVOLDUAL::generate_children_polytope
(std::vector<VERTEX_INDEX> & ivolpoly_vert,
 int iv[], int iw[], int corner[])
{
	std::vector<VERTEX_INDEX> new_ivolpoly_vert;
	bool flag_reverse[8] = {false, true, true, false, 
		                      true, false, false, true};
	
	new_ivolpoly_vert.push_back(iv[0]);
	new_ivolpoly_vert.push_back(iw[1]);
	new_ivolpoly_vert.push_back(iw[2]);
	new_ivolpoly_vert.push_back(iw[3]);
	new_ivolpoly_vert.push_back(iw[4]);
	new_ivolpoly_vert.push_back(iw[5]);
	new_ivolpoly_vert.push_back(iw[6]);
	new_ivolpoly_vert.push_back(iw[7]);

	new_ivolpoly_vert.push_back(iw[1]);
	new_ivolpoly_vert.push_back(iv[1]);
	new_ivolpoly_vert.push_back(iw[3]);
	new_ivolpoly_vert.push_back(iv[3]);
	new_ivolpoly_vert.push_back(iw[5]);
	new_ivolpoly_vert.push_back(iv[5]);
	new_ivolpoly_vert.push_back(iw[7]);
	new_ivolpoly_vert.push_back(iv[7]);

	new_ivolpoly_vert.push_back(iw[2]);
	new_ivolpoly_vert.push_back(iw[3]);
	new_ivolpoly_vert.push_back(iv[2]);
	new_ivolpoly_vert.push_back(iv[3]);
	new_ivolpoly_vert.push_back(iw[6]);
	new_ivolpoly_vert.push_back(iw[7]);
	new_ivolpoly_vert.push_back(iv[6]);
	new_ivolpoly_vert.push_back(iv[7]);

	new_ivolpoly_vert.push_back(iw[4]);
	new_ivolpoly_vert.push_back(iw[5]);
	new_ivolpoly_vert.push_back(iw[6]);
	new_ivolpoly_vert.push_back(iw[7]);
	new_ivolpoly_vert.push_back(iv[4]);
	new_ivolpoly_vert.push_back(iv[5]);
	new_ivolpoly_vert.push_back(iv[6]);
	new_ivolpoly_vert.push_back(iv[7]);

	if (flag_reverse[corner[0]]) {
		reverse(new_ivolpoly_vert.begin(), new_ivolpoly_vert.end());
	}
	ivolpoly_vert.insert(ivolpoly_vert.end(), 
		new_ivolpoly_vert.begin(), new_ivolpoly_vert.end());
}
