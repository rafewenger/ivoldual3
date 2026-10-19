/// @file ivoldualIO.h
/// @brief IO classes and routines for ivoldual.
/// - Version 0.3.0

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

#ifndef _IVOLDUALIO_
#define _IVOLDUALIO_

#include <ctime>
#include <string>

#include "ijk.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkIO.tpp"
#include "ijkisoIO.tpp"
#include "ijkstring.tpp"

#include "ijkdualIO.tpp"

#include "ivoldual_types.h"
#include "ivoldual_datastruct.h"


namespace IVOLDUAL {

  // **************************************************
  // TYPE DEFINITIONS
  // **************************************************

  typedef float COLOR_TYPE;           //!< Color type.

  //! Nrrd header.
  typedef IJK::NRRD_DATA<int, AXIS_SIZE_TYPE> NRRD_HEADER; 

  typedef enum { OFF, PLY, VTK } OUTPUT_FORMAT;    //!< Output format.


  // **************************************************
  // IO INFORMATION
  // **************************************************

  //! IO information
  class IO_INFO:public IVOLDUAL_DATA_FLAGS {

  protected:
    void Init();

  public:
    SCALAR_ARRAY isovalue;        ///< List of isovalues.
    std::vector<std::string> isovalue_string;
    COORD_ARRAY grid_spacing;
    std::string input_filename;
    std::string output_filename;
    std::string output_filename_prefix;
    bool label_with_isovalue;
    std::string output_off_filename;
    std::string output_ply_filename;
    std::string output_vtk_filename;
    std::string output_iv_filename;
    bool are_output_filenames_set;
    std::string isotable_directory;
    bool flag_output_off;    ///< Output Geomview .off file.
    bool flag_output_ply;    ///< Output PLY file.
    bool flag_output_iv;     ///< Output OpenInventor file.
    bool flag_output_vtk;    ///< Output vtk file.
    bool flag_report_time;
    bool flag_report_info;
    bool flag_use_stdout;
    bool flag_nowrite;
    bool flag_silent;
    bool flag_no_warn;
    bool flag_subsample;
    bool flag_report_all_isov;
    std::string report_isov_filename;
    bool flag_report_all_ivol_poly;
    std::string report_ivol_poly_filename;
    bool flag_write_scalar;
    std::string write_scalar_filename;
    int subsample_resolution;
    bool flag_supersample;
    int supersample_resolution;
    bool flag_subdivide;
    bool flag_rm_diag_ambig;
    bool flag_add_outer_layer;
    bool flag_color_alternating;  ///< Color simplices in alternating cubes
    int region_length;

    /// Color isosurface boundary vertices.
    bool flag_color_vert;         

    /// List of high resolution arguments,
    ///   e.g., "-highres {coord list}".
    std::vector<std::string> high_resolution_option;

    /// Flag for options set.
    bool is_qei_method_set;
    bool is_tri4_position_method_set;
    bool is_file_format_set;
    bool is_flag_orient_in_set;

    /// Return number of selected output formats.
    int NumOutputFormats() const;

    /// Copy io_info.
    void Set(const IO_INFO & io_info);

    /// Set output format to output_format.
    /// - Note: More than one output format may be set.
    void SetOutputFormat(const OUTPUT_FORMAT output_format);

    /// Set output filename for output format which is flagged true.
    /// @pre At most one output format is flagged true.
    void SetOutputFilename(const char * output_filename);

    /// Set output filename for output format which is flagged true.
    /// @pre At most one output format is flagged true.
    void SetOutputFilename(const std::string & output_filename)
    { SetOutputFilename(output_filename.c_str()); }

    /// Set output filename for a given output format.
    void SetOutputFilename
    (const OUTPUT_FORMAT output_format, const char * output_filename);

    /// Set output filename for a given output format.
    void SetOutputFilename
    (const OUTPUT_FORMAT output_format, const std::string & output_filename)
    { SetOutputFilename(output_format, output_filename.c_str()); }

    /// Construct output filenames.
    /// @param i Construct filenames for isovalue i.
    void ConstructOutputFilenames(const int i);


  public:
    IO_INFO() { Init(); };
    ~IO_INFO() { Init(); };
  };


  // **************************************************
  // OUTPUT INFORMATION
  // **************************************************

  /// Interval volume output data.
  class IVOL_OUTPUT_DATA:public IVOLDUAL_DATA_FLAGS {

  public:
    /* OBSOLETE
    bool flag_interval_volume;
    */
    int num_vertices_per_isopoly;
    bool flag_report_all_ivol_poly;
    bool use_triangle_mesh;
    bool is_flag_orient_in_set;
    /* OBSOLETE
    bool flag_orient_in;
    bool flag_rm_non_manifold;
    */
    std::string report_isov_filename;
    std::string report_ivol_poly_filename;

    IVOL_OUTPUT_DATA():
      /* OBSOLETE
      flag_interval_volume(true),
      */
      num_vertices_per_isopoly(8),
      is_flag_orient_in_set(false)
      /* OBSOLETE
      flag_report_all_ivol_poly(false),
      is_flag_orient_in_set(false),
      flag_rm_non_manifold(false)
      */
    {}
      
  };
    
  /// Output information.
  class OUTPUT_INFO:
    public IJK::ISO_OUTPUT_PARAM_BASE
  <SCALAR_TYPE,COORD_TYPE,IVOL_OUTPUT_DATA>
  {
  public:
    OUTPUT_INFO() {};
    ~OUTPUT_INFO() {};

    void SetDimension(const int d, const int numv_per_isopoly)
    {
      this->dimension = d;
      this->num_vertices_per_isopoly = numv_per_isopoly;
    }

    template <typename IO_INFO_TYPE>
    void Set(const IO_INFO_TYPE & io_info)
    {
      // TO BE CONTINUED...
    }

    void ConstructOutputFilenames(const int i)
    { this->IJK::ISO_OUTPUT_PARAM_BASE
        <SCALAR_TYPE,COORD_TYPE,IVOL_OUTPUT_DATA>
        ::ConstructOutputFilenames(i, true); }
  };
    
  /* OBSOLETE
  /// Output information.
  class OUTPUT_INFO:public IO_INFO {

  protected:
    void Init();

  public:
    typedef IVOLDUAL::COLOR_TYPE COLOR_TYPE;
    
  public:
    int dimension;
    int num_vertices_per_isopoly;
    SCALAR_TYPE isovalue[2];
    COORD_ARRAY grid_spacing;
    int grow_factor;
    int shrink_factor;

    OUTPUT_INFO() { Init(); };
    ~OUTPUT_INFO() { Init(); };

    void SetDimension(const int d, const int numv_per_isopoly);
  };
  */
 

  // **************************************************
  // TIMING FUNCTIONS/CLASSES
  // **************************************************

  /// Elapsed CPU time.
  class ELAPSED_CPU_TIME {

  protected:
    clock_t t;

  public:
    ELAPSED_CPU_TIME() { t = clock(); };

    clock_t getElapsed() {
      clock_t old_t = t;
      t = clock();
      return(t - old_t);
    };
  };

  /// Elapsed wall time.
  class ELAPSED_TIME {

  protected:
    time_t t;

  public:
    ELAPSED_TIME() { time(&t);  };

    double getElapsed() {
      time_t old_t = t;
      time(&t);
      return(difftime(t,old_t));
    };
  };

  /// IO time.
  struct IO_TIME {
    double read_nrrd_time;  ///< Wall time to read nrrd file.
    double write_time;      ///< Wall time to write output.
  };


  // **************************************************
  // PARSE COMMAND LINE
  // **************************************************

  /// Parse the command line.
  void parse_command_line(int argc, char **argv, IO_INFO & io_info);

  /// Check input information in io_info
  bool check_input
    (const IO_INFO & io_info, 
     const DUALISO_SCALAR_GRID_BASE & scalar_grid,
     IJK::ERROR & error);


  // **************************************************
  // OUTPUT DUAL INTERVAL VOLUME
  // **************************************************

  /// Output dual interval volume
  void output_dual_interval_volume
  (const OUTPUT_INFO & output_info, 
   const IVOLDUAL_DATA & ivoldual_data,
   const DUAL_INTERVAL_VOLUME & interval_volume,
   const IVOLDUAL_INFO & ivoldual_info, IO_TIME & io_time);

  /// Output dual interval volume simplices.
  void output_dual_interval_volume_simplices
  (const OUTPUT_INFO & output_info, 
   const IVOLDUAL_DATA & ivoldual_data,
   const COORD_ARRAY & vertex_coord,
   const VERTEX_INDEX_ARRAY & simplex_vert,
   const IVOLDUAL_INFO & ivoldual_info, IO_TIME & io_time);

  /// Output dual interval volume hexahedra.
  void output_dual_interval_volume
  (const OUTPUT_INFO & output_info, 
   const IVOLDUAL_DATA & ivoldual_data,
   const COORD_ARRAY & vertex_coord,
   const VERTEX_INDEX_ARRAY & hex_vert,
   const IVOLDUAL_INFO & ivoldual_info, IO_TIME & io_time);


  // **************************************************
  // READ NEARLY RAW RASTER DATA (nrrd) FILE
  // **************************************************

  /// Read a nearly raw raster data (nrrd) file.
  void read_nrrd_file
    (const char * input_filename, DUALISO_SCALAR_GRID & scalar_grid, 
     NRRD_HEADER & nrrd_header, IO_TIME & io_time);

  /// Read a nearly raw raster data (nrrd) file.
  /// - Pass C++ string instead of char * input filename.
  void read_nrrd_file
  (const std::string & input_filename, DUALISO_SCALAR_GRID & scalar_grid, 
   NRRD_HEADER & nrrd_header, IO_TIME & io_time);


  // **************************************************
  // WRITE NEARLY RAW RASTER DATA (nrrd) FILE
  // **************************************************

  /// Write scalar grid to nrrd file output_filename.
  void write_nrrd_file
  (const char * output_filename, const DUALISO_SCALAR_GRID_BASE & grid, 
   const bool flag_gzip);

  /// Write scalar grid to nrrd file output_filename.
  /// - Pass C++ string instead of char * output filename.
  void write_nrrd_file
  (const std::string & output_filename, const DUALISO_SCALAR_GRID_BASE & grid, 
   const bool flag_gzip);


  // **************************************************
  // RESCALE ROUTINES
  // **************************************************

  /// Rescale vertex coordinates by grid_spacing.
  void rescale_vertex_coord
    (const int dimension, const COORD_TYPE * grid_spacing,
     std::vector<COORD_TYPE> & vertex_coord);

  /// Rescale vertex coordinates by grid_spacing.
  /// @pre grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord(const std::vector<COORD_TYPE> & grid_spacing,
                            std::vector<COORD_TYPE> & vertex_coord);

  /// Rescale subsampled/supersampled vertex coordinates.
  void rescale_vertex_coord
  (const int grow_factor, const int shrink_factor, COORD_ARRAY & vertex_coord);

  /// Rescale vertex coordinates by grow and shrink factor and by grid_spacing.
  /// @pre grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord
    (const int grow_factor, const int shrink_factor,
     const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord);


  // **************************************************
  // WRITE_DUAL_MESH
  // **************************************************

  /// Write dual mesh with output format output_format.
  /// @param output_format Output format.
  void write_dual_mesh
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & plist);

  /// Write dual mesh.
  void write_dual_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist);

  /// Write dual mesh and record output time.
  void write_dual_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     IO_TIME & io_time);

  /// Write dual mesh and color facets with output format output_format.
  void write_dual_mesh_color
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & plist,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// Write dual mesh and color facets.
  void write_dual_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// Write dual mesh, color facets, and record output time.
  void write_dual_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     IO_TIME & io_time);

  /// Write dual isosurface triangular mesh
  /// @param output_info Output information.
  /// @param output_format Output format.
  /// @param vertex_coord List of vertex coordinates.
  /// @param tri_vert[] List of triangle vertices.
  ///        tri_vert[3*i+k] is k'th vertex of triangle i.
  void write_dual_tri_mesh
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert);

  /// Write dual isosurface triangular mesh.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  /// @param tri_vert[] List of triangle vertices.
  ///        tri_vert[3*i+k] is k'th vertex of triangle i.
  void write_dual_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert);

  /// Write dual isosurface triangular mesh.
  /// Record write time.
  void write_dual_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   IO_TIME & io_time);

  /// Write dual isosurface mesh of quad and triangles.
  /// @param output_info Output information.
  /// @param output_format Output format.
  /// @param vertex_coord List of vertex coordinates.
  void write_dual_quad_tri_mesh
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert);

  /// Write dual isosurface mesh of quad and triangles.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  void write_dual_quad_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert);

  /// Write dual isosurface mesh of quad and triangles and record output time.
  void write_dual_quad_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   IO_TIME & io_time);

  /// Write dual isosurface triangular mesh and color vertices
  ///   with output format output_format.
  /// @param output_info Output information.
  /// @param output_format Output format.
  /// @param vertex_coord List of vertex coordinates.
  /// @param tri_vert[] List of triangle vertices.
  ///        tri_vert[3*i+k] is k'th vertex of triangle i.
  void write_dual_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// Write dual isosurface triangular mesh, color vertices.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  /// @param tri_vert[] List of triangle vertices.
  ///        tri_vert[3*i+k] is k'th vertex of triangle i.
  void write_dual_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// Write dual isosurface triangular mesh, color vertices 
  ///   and record output time.
  void write_dual_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
   IO_TIME & io_time);

  /// Write dual isosurface mesh of quad and triangles and color vertices
  ///   with output format output_format.
  /// @param output_info Output information.
  /// @param output_format Output format.
  /// @param vertex_coord List of vertex coordinates.
  void write_dual_quad_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// Write dual isosurface mesh of quad and triangles and color vertices.
  /// @param output_info Output information.
  /// @param vertex_coord List of vertex coordinates.
  void write_dual_quad_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// Write dual isosurface mesh of quad and triangles, color vertices,
  ///   and report output time.
  void write_dual_quad_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
   IO_TIME & io_time);


  // **************************************************
  // SET ROUTINES
  // **************************************************

  /// Set output_info based on isotable, io_info and isovalue index i.
  void set_output_info
  (const IO_INFO & io_info, 
   const int i, OUTPUT_INFO & output_info);

  /// Set simplices in alternating cubes to have different colors.
  void set_color_alternating
  (const DUALISO_GRID & grid, const std::vector<VERTEX_INDEX> & cube_list, 
   COLOR_TYPE * color);


  // **************************************************
  // REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // **************************************************

  /// Report scalar field or interval volume information.
  void report_ivol_info
  (const OUTPUT_INFO & output_info, 
   const IVOLDUAL_DATA & ivoldual_data,
   const COORD_ARRAY & vertex_coord, 
   const VERTEX_INDEX_ARRAY & plist, 
   const IVOLDUAL_INFO & ivoldual_info);

  void report_num_cubes
    (const DUALISO_GRID & full_grid, const IO_INFO & io_info, 
     const IVOLDUAL_DATA & ivoldual_data);

  void report_num_cubes
    (const DUALISO_GRID & full_grid, const IO_INFO & io_info, 
     const DUALISO_GRID & dualiso_data_grid);

  /// Report number of changes for eliminating non-manifold
  void report_non_manifold_changes(const IVOLDUAL_INFO & dualiso_info);

  void warn_non_manifold(const IO_INFO & io_info);


  // **************************************************
  // REPORT INTERVAL VOLUME VERTICES AND POLYTOPES
  // **************************************************

  /// Report all interval volume vertices.
  void report_all_ivol_vert
  (std::ostream & out,
   const DUALISO_GRID & grid, 
   const DUAL_INTERVAL_VOLUME & interval_volume);

  /// Report all interval volume vertices.
  /// - Version which opens output stream from output info.
  void report_all_ivol_vert
  (const OUTPUT_INFO & output_info,
   const DUALISO_GRID & grid, const DUAL_INTERVAL_VOLUME & interval_volume);

  /// Report all interval volume polytopes.
  void report_all_ivol_poly
  (std::ostream & out, const DUALISO_GRID & grid,  
   const DUAL_INTERVAL_VOLUME & interval_volume);

  /// Report all interval volume hexahedra.
  /// - Version which opens output stream from output info.
  void report_all_ivol_hex
  (const OUTPUT_INFO & output_info,
   const DUALISO_GRID & grid, const DUAL_INTERVAL_VOLUME & interval_volume);

  /// Report all interval volume hexahedra.
  void report_all_ivol_hex
  (std::ostream & out, const DUALISO_GRID & grid,  
   const DUAL_INTERVAL_VOLUME & interval_volume);


  // **************************************************
  // REPORT TIMING INFORMATION
  // **************************************************

  void report_dualiso_time
    (const IO_INFO & io_info, const DUALISO_TIME & dualiso_time, 
     const char * mesh_type_string);

  void report_time
    (const IO_INFO & io_info, const IO_TIME & io_time, 
     const DUALISO_TIME & dualiso_time, const double total_elapsed_time);


  // **************************************************
  // USAGE/HELP MESSAGES
  // **************************************************

  void usage_error();
  void usage(std::ostream & out, const int return_code);
  void usage_all(std::ostream & out, const int return_code);
  void help(), help_all();

}

#endif
