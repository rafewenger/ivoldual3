/// @file ivoldualtable.cxx
/// @brief ivoldual lookup table.
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

#include <vector>
#include "ivoldualtable.h"

// **************************************************
// CLASS IVOLDUAL_CUBE_TABLE MEMBER FUNCTIONS
// **************************************************

void IVOLDUAL::IVOLDUAL_CUBE_TABLE::Init()
{
  vertex_info.Set(*this);
  compute_ivoldual_table_num_incident_poly(*this, vertex_info);
  compute_ivoldual_table_num_incident_isosurface_poly(*this, vertex_info);
  determine_ivol_vertex_connection_directions(*this, vertex_info);
  determine_ivol_vertex_iso_connection_directions(*this, vertex_info);
  determine_separation_vertices(*this, vertex_info);
  determine_separation_edges(*this, vertex_info);
  determine_doubly_connected_ivol3D_vertices(*this, vertex_info);
}
