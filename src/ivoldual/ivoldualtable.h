/// @file ivoldualtable.h
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

#ifndef _IVOLDUALTABLE_
#define _IVOLDUALTABLE_

#include "ijkcube.tpp"
#include "ijkdualtableX.txx"

#include "ivoldual_types.h"


namespace IVOLDUAL {

  // **************************************************
  // INTERVAL VOLUME LOOKUP TABLE ENTRY
  // **************************************************

  class IVOLDUAL_TABLE_ENTRY:public 
  IJKDUALTABLE::IVOLDUAL_TABLE_ENTRY
  <int,ISO_VERTEX_INDEX,FACET_BITS_TYPE,TABLE_INDEX> 
  {
  public:
    // ADDITIONAL DATA ABOUT THE IVOLDUAL CONFIGURATION.
  };

  // **************************************************
  // INTERVAL VOLUME VERTEX INFORMATION
  // **************************************************

  class IVOLDUAL_TABLE_VERTEX_INFO {

  public:
    typedef unsigned char CUBE_VERTEX_TYPE;
    typedef unsigned char CUBE_EDGE_TYPE;
    typedef unsigned char DIR_BITS_TYPE;

    const static CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX = 255;
    const static CUBE_EDGE_TYPE UNDEFINED_CUBE_EDGE = 255;

    int num_incident_poly;
    int num_incident_isopoly;

    /// Separation vertex.
    CUBE_VERTEX_TYPE separation_vertex;

    /// Separation edge.
    CUBE_EDGE_TYPE separation_edge;

    /// If ivol vertex is doubly connected, then ivol vertex is
    ///   doubly connected across facet doubly_connected_facet.
    FACET_INDEX doubly_connected_facet;

    bool is_doubly_connected;        ///< True, if vertex is doubly connected.

    /// Bit flag for connection direction.
    /// If bit i is 1, then vertex connects across facet i.
    DIR_BITS_TYPE connect_dir;

    /// Bit flag for connection direction on isosurface.
    /// If bit i is 1, then vertex connects across facet i
    ///   and edge across facet i is on the isosurface.
    DIR_BITS_TYPE iso_connect_dir;

    /// Number of interval volume polytopes incident on vertex.
    int NumIncidentPoly() const
    { return(num_incident_poly); }

    /// Number of isosurface polytopes incident on vertex.
    int NumIncidentIsoPoly() const
    { return(num_incident_isopoly); }

  };


  // **************************************************
  // INTERVAL VOLUME LOOKUP TABLE
  // **************************************************

  class IVOLDUAL_CUBE_TABLE:
    public IJKDUALTABLE::IVOLDUAL_CUBE_DOUBLE_TABLE
  <4,int,int,TABLE_INDEX,IVOLDUAL_TABLE_ENTRY>
  {
  
  protected:
    void Init();

  public:
    IJKDUALTABLE::DUAL_TABLE_VERTEX_INFO<int,IVOLDUAL_TABLE_VERTEX_INFO>
    vertex_info;

  public:
    IVOLDUAL_CUBE_TABLE(const int dimension, const bool flag_separate_neg):
      IVOLDUAL_CUBE_DOUBLE_TABLE(dimension, flag_separate_neg)
    { Init(); }; 

    const IVOLDUAL_TABLE_VERTEX_INFO & VertexInfo
    (const TABLE_INDEX ientry, const int ivolv) const
    { return(vertex_info.VertexInfo(ientry, ivolv)); }
  };

}

#endif
