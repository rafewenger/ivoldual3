/// @file ivoldual_main.cxx
/// @brief Generate interval volume using dual contouring algorithm.
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


#include <iostream>

#include "isodual.h"
#include "ivoldualIO.h"
#include "ivoldual.h"

#include "ivoldual_triangulate.h"
#include "ivoldual_reposition.h"

using namespace IJK;
using namespace IJKDUAL;
using namespace ISODUAL;
using namespace IVOLDUAL;

using namespace std;

// local subroutines
void memory_exhaustion();
void eliminate_non_manifold
(const IO_INFO & io_info, IVOLDUAL_DATA & ivoldual_data, 
 IVOLDUAL_INFO & dualiso_info);
void construct_interval_volume
(const IO_INFO & io_info, const IVOLDUAL_DATA & ivoldual_data,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time, IVOLDUAL_INFO & dualiso_info);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);

  DUALISO_TIME dualiso_time;
  IO_TIME io_time = {0.0, 0.0};
  IO_INFO io_info;
  IJK::ERROR error;

  try {

    std::set_new_handler(memory_exhaustion);

    // Set to interval volume.
    io_info.flag_interval_volume = true;

    parse_command_line(argc, argv, io_info);

    DUALISO_SCALAR_GRID full_scalar_grid, scalar_grid_4D;
    NRRD_HEADER nrrd_header;
    read_nrrd_file
      (io_info.input_filename, full_scalar_grid,  nrrd_header, io_time);

    if (!check_input(io_info, full_scalar_grid, error)) 
      { throw(error); };

    nrrd_header.GetSpacing(io_info.grid_spacing);

    // set DUAL datastructures and flags
    IVOLDUAL_DATA ivoldual_data;

    // Note: ivoldual_data.SetScalarGrid must be called before set_mesh_data.
    // subsample and supersample parameters are hard-coded here.
    ivoldual_data.SetScalarGrid
      (full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution, 
       io_info.flag_supersample, io_info.supersample_resolution, 
       io_info.flag_subdivide, io_info.flag_rm_diag_ambig,
       io_info.flag_add_outer_layer,
       io_info.default_interior_code,
       io_info.isovalue[0], io_info.isovalue[1]);

    // set ivoldual info
    int dimension = ivoldual_data.ScalarGrid().Dimension();
    IVOLDUAL_INFO dualiso_info(dimension);

    // Eliminate non-manifold facets and cubes
    if (io_info.flag_rm_non_manifold) {
      eliminate_non_manifold
        (io_info, ivoldual_data, dualiso_info);
    }

    ivoldual_data.Set(io_info);
    warn_non_manifold(io_info);
    report_num_cubes(full_scalar_grid, io_info, ivoldual_data);

    if (io_info.flag_write_scalar) {
      write_nrrd_file
        (io_info.write_scalar_filename, ivoldual_data.ScalarGrid(), false);
    }

    construct_interval_volume
      (io_info, ivoldual_data, dualiso_time, io_time, dualiso_info);

    /* OBSOLETE.  MOVED TO report_ivol_info.
    // print out total number of changes for eliminating non-manifold
    if (io_info.flag_rm_non_manifold) {
      report_non_manifold_changes(dualiso_info);
    }
    */
    
    if (io_info.flag_report_time) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cout << endl;
      report_time(io_info, io_time, dualiso_time, total_elapsed_time);
    };

  } 
  catch (ERROR & error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

void eliminate_non_manifold
(const IO_INFO & io_info, IVOLDUAL_DATA & ivoldual_data, 
 IVOLDUAL_INFO & dualiso_info)
{
  eliminate_non_manifold_grid
    (ivoldual_data, io_info.isovalue[0], io_info.isovalue[1], dualiso_info);
}

void construct_interval_volume
(const IO_INFO & io_info, const IVOLDUAL_DATA & ivoldual_data,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time, IVOLDUAL_INFO & dualiso_info)
{
  int dimension = ivoldual_data.ScalarGrid().Dimension();
  const int num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
  const int num_cubes = ivoldual_data.ScalarGrid().ComputeNumCubes();

  io_time.write_time = 0;

  for (unsigned int i = 0; i+1 < io_info.isovalue.size(); i++) {

    const SCALAR_TYPE isovalue0 = io_info.isovalue[i];
    const SCALAR_TYPE isovalue1 = io_info.isovalue[i+1];

    dualiso_info.grid.num_cubes = num_cubes;

    // Dual contouring.  
    DUAL_INTERVAL_VOLUME interval_volume(dimension, num_cube_vertices);

    dual_contouring_interval_volume
      (ivoldual_data, isovalue0, isovalue1, interval_volume,
       dualiso_info);

    // Time info
    dualiso_time.Add(dualiso_info.time);

    // Rescale vertex coordinates. 
    rescale_vertex_coord
      (dimension, ivoldual_data.ScalarGrid().SpacingPtrConst(),
       interval_volume.vertex_coord);

    if (ivoldual_data.UseTriangleMesh()) {
      triangulate_interval_volume(ivoldual_data, interval_volume);
    }


    OUTPUT_INFO output_info;
    output_info.SetDimension(dimension, num_cube_vertices);
    set_output_info(io_info, i, output_info);

    output_dual_interval_volume
      (output_info, ivoldual_data, interval_volume, dualiso_info, io_time);

    if (output_info.flag_report_all_isov) {
      report_all_ivol_vert
        (output_info, ivoldual_data.ScalarGrid(), interval_volume);
    }

    if (output_info.flag_report_all_ivol_poly) {
      report_all_ivol_hex
        (output_info, ivoldual_data.ScalarGrid(), interval_volume);
    }
  }
}


void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

