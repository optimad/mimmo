/*---------------------------------------------------------------------------*\
 * 
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

#ifndef PORT_DEFINITIONS_HPP
#define PORT_DEFINITIONS_HPP

/*!
 * \ingroup Ports
 * \{
 */

//PORTS DEFINITION AS CONSTANTS
#define M_GEOM            "M_GEOM"              /**< Port dedicated to communication of pointers to a MimmoObject object*/
#define M_GEOM2           "M_GEOM2"             /**< Port dedicated to communication of pointers to a MimmoObject object*/
#define M_GEOM3           "M_GEOM3"             /**< Port dedicated to communication of pointers to a MimmoObject object*/
#define M_GEOM4           "M_GEOM4"             /**< Port dedicated to communication of pointers to a MimmoObject object*/
#define M_GEOMOFOAM       "M_GEOMOFOAM"         /**< Port dedicated to communication of pointers to a MimmoObject object used as I/O between OFOAM blocks*/
#define M_GEOMOFOAM2      "M_GEOMOFOAM2"        /**< Port dedicated to communication of pointers to a MimmoObject object used as I/O between OFOAM blocks*/
#define M_VECGEOM         "M_VECGEOM"           /**< Port dedicated to communication of list of pointers to a MimmoObject object [ std::vector< MimmoObject* > ] */
#define M_COORDS          "M_COORDS"            /**< Port dedicated to communication of coordinates of points [vector < array < double,3 > > ] */
#define M_GLOBAL          "M_GLOBAL"            /**< Port dedicated to communication of coordinates of points in a global reference system [ vector < array < double,3 > > ] */
#define M_LOCAL           "M_LOCAL"             /**< Port dedicated to communication of coordinates of points in a local reference system [ vector < array < double,3 > > ] */
#define M_DISPLS          "M_DISPLS"            /**< Port dedicated to communication of displacements of 3D points [ vector < array < double,3 > > ] */
#define M_GDISPLS         "M_GDISPLS"           /**< Port dedicated to communication of displacements relative to geometry vertices [ vector < array < double,3 > > ] */
#define M_GDISPLS2        "M_GDISPLS2"           /**< Port dedicated to communication of displacements relative to geometry vertices [ vector < array < double,3 > > ] */
#define M_FILTER          "M_FILTER"            /**< Port dedicated to communication of a scalar field used as a filter function [ vector < double > ] */
#define M_FILTER2         "M_FILTER2"            /**< Port dedicated to communication of a scalar field used as a filter function [ vector < double > ] */
#define M_DATAFIELD       "M_DATAFIELD"         /**< Port dedicated to communication of a scalar field used as a filter function [ vector < double > ] */
#define M_VECTORSI        "M_VECTORSI"          /**< Port dedicated to communication of a generic list of short integers [ vector < short int > ] */
#define M_VECTORLI        "M_VECTORLI"          /**< Port dedicated to communication of a generic list of long integers [ vector < long int > ] */
#define M_SCALARFIELD     "M_SCALARFIELD"       /**< Port dedicated to communication of a generic scalar field [ MimmoPiercedvector < double > ] */
#define M_SCALARFIELD2    "M_SCALARFIELD2"      /**< Port dedicated to communication of a generic scalar field [ MimmoPiercedvector < double > ] */
#define M_VECTORFIELD     "M_VECTORFIELD"       /**< Port dedicated to communication of a generic vector field [ MimmoPiercedvector < array< double,3> > ] */
#define M_VECTORFIELD2    "M_VECTORFIELD2"      /**< Port dedicated to communication of a generic vector field [ MimmoPiercedvector < array< double,3> > ] */
#define M_VECSFIELDS      "M_VECSFIELDS"        /**< Port dedicated to communication of a std::vector of generic scalar fields [ std::vector< MimmoPiercedvector < double > >] */
#define M_VECVFIELDS      "M_VECVFIELDS"        /**< Port dedicated to communication of a std::vector of generic vector fields [ std::vector< MimmoPiercedvector < array< double,3> > >] */
#define M_POINT           "M_POINT"             /**< Port dedicated to communication of 3D point coordinates [ array < double,3 > ] */
#define M_AXIS            "M_AXIS"              /**< Port dedicated to communication of an axis direction [ array < double,3 > ] */
#define M_AXES            "M_AXES"              /**< Port dedicated to communication of a 3 axis reference system [ array < array< double,3 >, 3>] */
#define M_SPAN            "M_SPAN"              /**< Port dedicated to communication of the span(dimensions) of an elemental object [ array < double,3 > ] */
#define M_DIMENSION       "M_DIMENSION"         /**< Port dedicated to communication of the lattice mesh dimension [ array < int,3 > ] */
#define M_INFLIMITS       "M_INFLIMITS"         /**< Port dedicated to communication of the nferior limits for building an elemental shape [ array < double,3 > ]*/
#define M_SHAPE           "M_SHAPE"             /**< Port dedicated to communication of an elemental shape type [ mimmo::ShapeType enum ] */
#define M_COPYSHAPE       "M_COPYSHAPE"         /**< Port dedicated to communication of an elemental shape itself, passed by pointer, instantiated elsewhere [BasicShape *]*/
#define M_SHAPEI          "M_SHAPEI"            /**< Port dedicated to communication of an elemental shape type [int -> casted from mimmo::ShapeType enum] */
#define M_PLANE           "M_PLANE"             /**< Port dedicated to communicate planes in a*x+b*y+c*z +d = 0 implicit form [ array < double,4 > ]*/  
#define M_VALUED          "M_VALUED"            /**< Port dedicated to communication of a single scalar value [double] */
#define M_VALUEI          "M_VALUEI"            /**< Port dedicated to communication of a single scalar value [int] */
#define M_VALUEB          "M_VALUEB"            /**< Port dedicated to communication of a single scalar value [bool] */
#define M_BMATRIX         "M_BMATRIX"           /**< Port dedicated to communication of a 3x3 matrix of float values [array< array< double,3>,3>] */
#define M_BCOEFFS         "M_BCOEFFS"           /**< Port dedicated to communication of a 3x3 matrix of float vector values [array < array< vector< double>,3>,3>] */
#define M_VALUESI         "M_VALUESI"           /**< Port dedicated to communication of a single scalar value [short int] */
#define M_DEG             "M_DEG"               /**< Port dedicated to communication of degrees of freedom numbers in each of 3D space dimensions [array< int,3 >] */
#define M_NURBSCOORDTYPE  "M_NURBSCOORDTYPE"    /**< Port dedicated to communicate NURBS design boundary condition on FFDLattice [array < mimmo::CoordType,3> ] */
#define M_NURBSWEIGHTS    "M_NURBSWEIGHTS"      /**< Port dedicated to communicate NURBS weights on FFDLattice nodes [ vector < double > ] */
#define M_BCCGNS          "M_BCCGNS"            /**< Port dedicated to communicate a pointer to a BCCGNS object (Boundary Conditions Info for CGNS export) */
#define M_POINT2          "M_POINT2"            /**< Port dedicated to communication of 3D point coordinates [ array < double,3 > ] */
#define M_VALUED2         "M_VALUED2"           /**< Port dedicated to communication of a single scalar value [double] */
#define M_VALUEB2         "M_VALUEB2"           /**< Port dedicated to communication of a single scalar value [bool] */
#define M_VALUEI2         "M_VALUEI2"           /**< Port dedicated to communication of a single scalar value [int] */
#define M_POLYDATA_       "M_POLYDATA_"         /**< Port dedicated to communicate a pointer to a vtk polydata mesh [vtkPolyData *] */
#define M_UMAPIDS      	  "M_UMAPIDS"           /**< Port dedicated to communicate a map between IDs [std::unordered_map<long, long> *] */

/*!
 * \}
 */

/*!
 * \ingroup PortContainers
 * \{
 */

#define  MC_SCALAR      "MC_SCALAR"         /**< Single value container identifier */
#define  MC_VECTOR      "MC_VECTOR"         /**< std::vector< . > container identifier */
#define  MC_ARRAY3      "MC_ARRAY3"         /**< std::array< . , 3 > container identifier */
#define  MC_ARRAY4      "MC_ARRAY4"         /**< std::array< . , 4 > container identifier */
#define  MC_VECARR3     "MC_VECARR3"        /**< std::vector< std::array< . , 3 > > container identifier */
#define  MC_ARR3ARR3    "MC_ARR3ARR3"       /**< std::array< std::array< . , 3 > , 3 > container identifier */
#define  MC_ARR3ARR3VEC "MC_ARR3ARR3VEC"    /**< std::array< std::array< std::array< . , 3 > , 3 > , 3 > container identifier */
#define  MC_MPVECTOR    "MC_MPVECTOR"       /**< MimmoPiercedVector< . > container identifier */
#define  MC_MPVECARR3   "MC_MPVECARR3"      /**< MimmoPiercedVector< std::array< . , 3 > > container identifier */
#define  MC_UMAP        "MC_UMAP"    		/**< std::unordered_map< . , . > container identifier */

/*!
 * \}
 */

/*!
 * \ingroup PortData
 * \{
 */
#define  MD_MIMMO_                  "MD_MIMMO_"                  /**< mimmo::MimmoObject pointer data identifier*/
#define  MD_INT                     "MD_INT"                     /**< integer data identifier*/
#define  MD_SHORT                   "MD_SHORT"                   /**< short integer data identifier*/
#define  MD_LONG                    "MD_LONG"                    /**< long integer data identifier*/
#define  MD_FLOAT                   "MD_FLOAT"                   /**< float/double data identifier*/
#define  MD_BOOL                    "MD_BOOL"                    /**< boolean data identifier*/
#define  MD_SHAPET                  "MD_SHAPET"                  /**< mimmo::ShapeType data identifier*/
#define  MD_SHAPE_                  "MD_SHAPE_"                  /**< mimmo::BasicShape pointer data identifier*/
#define  MD_COORDT                  "MD_COORDT"                  /**< mimmo::CoordType data identifier*/
#define  MD_POLYDATA_               "MD_POLYDATA_"               /**< VTK::vtkPolyData* data identifier*/
#define  MD_BCCGNS_                 "MD_BCCGNS_"                 /**< mimmo::BCCGNS (Boundary Conditions Info for CGNS export class) pointer data identifier */
#define  MD_MPVECFLOAT              "MD_MPVECFLOAT"              /**< MimmoPiercedVector< double > data identifier */
#define  MD_MPVECARR3FLOAT          "MD_MPVECARR3FLOAT"          /**< MimmoPiercedVector< std::array< double,3>  > data identifier */
#define  MD_PAIRLONGLONG            "MD_PAIRLONGLONG"            /**< std::pair< long, long > data identifier */

/*!
 * \}
 */

#endif

