/// @file ijkdualIO.txx
/// @brief IO templates for ijkdual.
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

#ifndef _IJKDUALIO_TXX_
#define _IJKDUALIO_TXX_

#include <vector>

#include "ijkdual_types.h"

namespace IJKDUAL {

  // **************************************************
  // OUTPUT ISOSURFACE
  // **************************************************

  /// Output dual isosurface.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (output_info.use_triangle_mesh) {
      output_dual_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.tri_vert, dualiso_info, io_time);
    }
    else if (output_info.flag_dual_collapse) {
      output_dual_quad_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dual_isosurface.tri_vert, 
         dualiso_info, io_time);
    }
    else {
      output_dual_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dualiso_info, io_time);
    }
  }


  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & slist,
   const DUALISO_INFO_TYPE & dualiso_info, 
   IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info(output_info, dualiso_data, 
                      vertex_coord, slist, dualiso_info);
    }

    if (!output_info.flag_nowrite) 
      { write_dual_mesh(output_info, vertex_coord, slist, io_time); }
  }


  /// Output isosurface of triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_tri_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info(output_info, dualiso_data,
                      vertex_coord, tri_vert, dualiso_info);
    }

    if (!output_info.flag_nowrite) {
      write_dual_tri_mesh(output_info, vertex_coord, tri_vert, io_time);
    }
  }


  /// Output isosurface of quadrilaterals and triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_quad_tri_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_quad_tri_iso_info
        (output_info, dualiso_data, vertex_coord, quad_vert, tri_vert, 
         dualiso_info);
    }

    if (!output_info.flag_nowrite) {
      write_dual_quad_tri_mesh
        (output_info, vertex_coord, quad_vert, tri_vert, io_time);
    }
  }


  /// Output isosurface with colored facets.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename COLOR_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_isosurface_color
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & slist,
   const COLOR_TYPE * front_color, 
   const COLOR_TYPE * back_color,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info
        (output_info, dualiso_data, vertex_coord, slist, dualiso_info);
    }
  
    if (!output_info.flag_nowrite) {
      write_dual_mesh_color
        (output_info, vertex_coord, slist, front_color, back_color, io_time);
    }
  }


  /// Output isosurface with colored facets.
  /// - Version with input DUAL_ISOSURFACE_TYPE.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename COLOR_TYPE, typename IO_TIME_TYPE>
  void output_dual_isosurface_color
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    output_dual_isosurface_color
      (output_info, dualiso_data, 
       dual_isosurface.vertex_coord, dual_isosurface.isopoly_vert,
       front_color, back_color, dualiso_info, io_time);
  }


  /// Output isosurface.  Color facets in alternating colors.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_isosurface_color_alternating
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    const VERTEX_INDEX num_poly = dual_isosurface.NumIsoPoly();

    typedef typename OUTPUT_INFO_TYPE::COLOR_TYPE COLOR_TYPE;

    IJK::ARRAY<COLOR_TYPE> front_color(4*num_poly);
    IJK::ARRAY<COLOR_TYPE> back_color(4*num_poly);
    set_color_alternating
      (dualiso_data.ScalarGrid(), dual_isosurface.cube_containing_isopoly,
       front_color.Ptr());
    set_color_alternating
      (dualiso_data.ScalarGrid(), dual_isosurface.cube_containing_isopoly,
       back_color.Ptr());

    output_dual_isosurface_color
      (output_info, dualiso_data, dual_isosurface.vertex_coord,
       dual_isosurface.isopoly_vert, 
       front_color.PtrConst(), back_color.PtrConst(),
       dualiso_info, io_time);
  }


  // **************************************************
  // REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // **************************************************

  /// Report information about cubes with multiple vertices.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_INFO_TYPE>  
  void report_multi_isov_info
  (const OUTPUT_INFO_TYPE & output_info, const DUALISO_INFO_TYPE & dualiso_info)
  {
    const char * indent4 = "    ";
    const char * indent6 = "      ";
    const char * indent8 = "        ";

    using namespace std;

    if (output_info.allow_multiple_iso_vertices) {
      cout << indent4 << "# active (non-empty) cubes: "
           << dualiso_info.scalar.num_non_empty_cubes << endl;
      cout << indent4 << "# cubes with single isosurface vertex: "
           << dualiso_info.multi_isov.num_cubes_single_isov << endl;
      cout << indent4 << "# cubes with multiple isosurface vertices: "
           << dualiso_info.multi_isov.num_cubes_multi_isov << endl;

      if (output_info.flag_split_non_manifold) {
        cout << indent4 << "# cubes changed to 2 iso vertices to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_non_manifold_split << endl;
        cout << indent4 
             << "# ambiguous ridge cubes changed to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_ambig_ridge_cubes_changed << endl;
        cout << indent4 
             << "# non-ambiguous ridge cubes changed to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_non_ambig_ridge_cubes_changed << endl;
      }

      if (output_info.flag_select_split) {
        cout << indent4 << "# cubes changed in selecting isosurface patch splits: "
             << dualiso_info.multi_isov.num_1_2_changed << endl;
      }

      if (output_info.flag_connect_ambiguous) {
        cout << indent4 << "# ambiguous cubes changed to connect isosurface: "
             << dualiso_info.multi_isov.num_connect_changed << endl;
      }
    }
  }


  /// Report isosurface information.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUALISO_INFO_TYPE>
  void report_iso_info
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & plist, 
   const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int dimension = output_info.dimension;
    const int numv_per_simplex = output_info.num_vertices_per_isopoly;

    const char * indent4 = "    ";

    using namespace std;

    VERTEX_INDEX numv = (vertex_coord.size())/dimension;
    VERTEX_INDEX num_poly = (plist.size())/numv_per_simplex;
    VERTEX_INDEX num_grid_cubes = dualiso_info.grid.num_cubes;

    if (output_info.flag_interval_volume) {
      cout << "  Interval volume [" 
           << output_info.isovalue[0] << ":"
           << output_info.isovalue[1] << "].  "
           << numv << " ivol vertices.  "
           << num_poly << " ivol polytopes." << endl;
    }
    else {
      cout << "  Isovalue " << output_info.isovalue[0] << ".  "
           << numv << " isosurface vertices.  "
           << num_poly << " isosurface polytopes." << endl;
    }

    report_multi_isov_info(output_info, dualiso_info);
  }


  /// Report information about isosurface quadrilaterals and triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUALISO_INFO_TYPE>
  void report_quad_tri_iso_info
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & quad_list, 
   const std::vector<VERTEX_INDEX> & tri_list, 
   const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int dimension = output_info.dimension;
    const int numv_per_simplex = output_info.num_vertices_per_isopoly;
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_TRI(3);
    const char * indent4 = "    ";

    using namespace std;

    VERTEX_INDEX numv = (vertex_coord.size())/dimension;
    VERTEX_INDEX num_quad = (quad_list.size())/NUM_VERT_PER_QUAD;
    VERTEX_INDEX num_tri = (tri_list.size())/NUM_VERT_PER_TRI;
    VERTEX_INDEX num_grid_cubes = dualiso_info.grid.num_cubes;

    cout << "  Isovalue " << output_info.isovalue[0] << ".  " 
         << numv << " isosurface vertices.  "
         << num_quad+num_tri << " isosurface polygons." << endl;
    cout << indent4 << num_quad << " quadrilaterals.  "
         << num_tri << " triangles." << endl;

    report_multi_isov_info(output_info, dualiso_info);
  }

}

#endif
