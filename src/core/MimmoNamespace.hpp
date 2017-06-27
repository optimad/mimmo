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
\*---------------------------------------------------------------------------*/
#ifndef __MIMMONAMESPACE_HPP__
#define __MIMMONAMESPACE_HPP__

#include "enum.hpp"
#include "bitpit.hpp"

BETTER_ENUM(PortType, int,	
            M_COORDS 			= 0, 
            M_GLOBAL 			= 1,
            M_LOCAL 			= 2,
            M_DISPLS 			= 10,
            M_GDISPLS 			= 11,
            M_FILTER            = 12,
            M_FILTER_           = 13,
            M_DATAFIELD         = 14,
            M_VECTORSI			= 16,
            M_VECTORLI			= 17,
            M_SCALARFIELD       = 18,
            M_VECTORFIELD       = 19,
            M_POINT				= 20,
            M_AXIS				= 21,
            M_AXES				= 22,
            M_SPAN				= 23,
            M_DIMENSION			= 24,
            M_INFLIMITS			= 25,
            M_SHAPE	 			= 26,
            M_COPYSHAPE			= 27,
            M_SHAPEI	 		= 28,
            M_PLANE				= 29,
            M_VALUED			= 30,
            M_VALUEI			= 31,
            M_VALUEB			= 32,
            M_BMATRIX			= 33,
            M_BCOEFFS			= 34,
            M_VALUESI			= 35,
            M_DEG				= 40,
            M_RANGE				= 41,
            M_BOOLS3			= 42,
            M_NURBSCOORDTYPE	= 43,
            M_NURBSWEIGHTS 		= 44,
            M_VECSFIELDS        = 60,
            M_VECVFIELDS        = 61,
            M_GEOM2				= 98,
            M_GEOM 				= 99,
            M_FINFO				= 102,
            M_FINFO2			= 103,
            M_BCCGNS            = 108,
            M_POINT2			= 120,
            M_VALUED2			= 130,
            M_VALUED3			= 131,
            M_VALUEB2			= 140,
            M_VALUEB3			= 141,
            M_VALUEB4			= 142,
            M_VALUEB5			= 143,
            M_VALUEI2			= 150,
            M_POLYDATA_         = 1100
        );


namespace mimmo{

class BaseManipulation;

template<typename value_t, typename id_t>
class MimmoPiercedVector;

/*!
 * \class FileDataInfo
 * \brief FileDataInfo is a struct to stock data relative to names of external files. 
 * \ingroup core
 * 
 *
 * FileDataInfo has three fields: an integer relative to the type of file stored, 
 * a string reporting the path to the file and a string reporting the name of the file
 * 
 */
struct FileDataInfo{
    int ftype;          /**< file type identifier*/
    std::string fname;  /**< file name*/
    std::string fdir;   /**< file directory*/
    
    FileDataInfo();
    virtual ~FileDataInfo();
    FileDataInfo(const FileDataInfo & other);
    FileDataInfo & operator=(const FileDataInfo & other);
};


/*!
* 
*  \brief Utilities to create port connections between executable blocks.
*  \ingroup core
* Here are collected enums and methods to create connections between blocks throughout established ports.
* All the available port types are collected in a special enum <B>PortType</B>, istantiated through 
* BETTER_ENUM macro, and each type reports a specific TAG in order to be retracked in each moment 
* for compatibility checking purposes.
* 
* A type of data is related to each TAG. Same type of data can be related to
* multiple TAG of ports but with different meaning.
* A port with a TAG can communicate with other TAG in function of its pre-coded
* compatibility.
* The ports defined in the base library mimmo must have prefix M_ before their name, and can be recall as 
* an enum writing PortType::M_xxxxx.
* 
* Ports actually available in <B>PortType</B> special enum are:
* 
* - <B>M_COORDS        </B>= 0    Port dedicated to communicate coordinates of points [vector<array<double,3>>].,
* - <B>M_GLOBAL        </B>= 1    Port dedicated to communicate coordinates of points in a global reference system [vector<array<double,3>>].,
* - <B>M_LOCAL         </B>= 2    Port dedicated to communicate coordinates of points in a local reference system. [vector<array<double,3>>].,
* - <B>M_DISPLS        </B>= 10   Port dedicated to communicate displacements of points [vector<array<double,3>>].,
* - <B>M_GDISPLS       </B>= 11   Port dedicated to communicate displacements (or a generic field of array<double,3>) of geometry vertex [MimmoPiercedVector<array<double,3>>].,
* - <B>M_FILTER        </B>= 12   Port dedicated to communicate a scalar field used as filter function [MimmoPiercedVector<double>].,
* - <B>M_FILTER_       </B>= 13   Port dedicated to communicate a pointer to a scalar field used as filter function [MimmoPiercedVector<double>*].,
* - <B>M_DATAFIELD     </B>= 14   Port dedicated to communicate a generic scalar field [std::vector<double>].,
* - <B>M_VECTORSI      </B>= 16   Port dedicated to communicate a generic list of short integers [vector<short int>].,
* - <B>M_VECTORLI      </B>= 17   Port dedicated to communicate a generic list of long integers [vector<long int>].,
* - <B>M_SCALARFIELD   </B>= 18   Port dedicated to communicate a generic scalar field [MimmoPiercedVector<double>].,
* - <B>M_VECTORFIELD   </B>= 19   Port dedicated to communicate a generic vector field [MimmoPiercedVector<array<double,3>>].,
* - <B>M_POINT         </B>= 20   Port dedicated to communicate the coordinate of a point field [array<double,3>].,
* - <B>M_AXIS          </B>= 21   Port dedicated to communicate the direction of an axis [array<double,3>].,
* - <B>M_AXES          </B>= 22   Port dedicated to communicate a reference system [array<array<double,3>,3>].,
* - <B>M_SPAN          </B>= 23   Port dedicated to communicate the span of elemental object [array<double,3>].,
* - <B>M_DIMENSION     </B>= 24   Port dedicated to communicate the dimension of lattice mesh [array<int,3>].,
* - <B>M_INFLIMITS     </B>= 25   Port dedicated to communicate a inferior limits for building elemental shape [array<double,3>].,
* - <B>M_SHAPE         </B>= 26   Port dedicated to communicate a type of an elemental shape itself [mimmo::ShapeType].,
* - <B>M_COPYSHAPE     </B>= 27   Port dedicated to communicate a an elemental shape itself, passed by pointer, instantiated elsewhere [BasicShape *].,
* - <B>M_SHAPEI        </B>= 28   Port dedicated to communicate a type of an elemental shape itself [int].,
* - <B>M_PLANE         </B>= 29   Port dedicated to communicate planes in a*x+b*y+c*z +d = 0 implicit form [array<double,4>].,	
* - <B>M_VALUED        </B>= 30   Port dedicated to communicate a scalar value [double].,
* - <B>M_VALUEI        </B>= 31   Port dedicated to communicate a scalar value [int].,
* - <B>M_VALUEB        </B>= 32   Port dedicated to communicate a condition value [bool].,
* - <B>M_BMATRIX       </B>= 33   Port dedicated to communicate a matrix of values [array<array<double,3>,3>].,
* - <B>M_BCOEFFS       </B>= 34   Port dedicated to communicate a matrix of vector of values [array<array<vector<double>,3>,3>].,
* - <B>M_VALUESI       </B>= 35   Port dedicated to communicate a scalar short int value [short int].,
* - <B>M_DEG           </B>= 40   Port dedicated to communicate a set of number of degrees of freedom [array<int,3>] .,
* - <B>M_RANGE         </B>= 41   Port dedicated to communicate a set of span interval [array<array<double,2>,3>] .,
* - <B>M_BOOLS3        </B>= 42   Port dedicated to communicate a set of conditions [array<bool,3>].,
* - <B>M_NURBSCOORDTYPE</B>= 43   Port dedicated to communicate condition to design NURBS on FFDLattice [array<mimmo::CoordType,3>] .,
* - <B>M_NURBSWEIGHTS  </B>= 44   Port dedicated to communicate condition to exchange NURBS weights on FFDLattice [vector<double>],
* - <B>M_VECSFIELDS    </B>= 60   Port dedicated to communicate a vector of scalar fields [vector<mimm::MimmoPiercedVector<double>].,
* - <B>M_VECVFIELDS    </B>= 61   Port dedicated to communicate a vector of vector (displacements) fields related to a geometry [vector<mimm::MimmoPiercedVector<double>].,
* - <B>M_VIOLATION     </B>= 82   Port dedicated to communicate a double value, with a BaseManipulation object which refers to as std::pair<BaseManipulation *, double>,
* - <B>M_GEOM2         </B>= 98   Port dedicated to communicate a pointer to a geometry [MimmoObject *],
* - <B>M_GEOM          </B>= 99   Port dedicated to communicate a pointer to a geometry [MimmoObject *],
* - <B>M_FINFO         </B>= 102  Port dedicated to communicate a list of FileInfoData classes std::vector<FileInfoData>,
* - <B>M_FINFO2        </B>= 103  Port dedicated to communicate a list of FileInfoData classes std::vector<FileInfoData>.,	 
* - <B>M_BCCGNS        </B>= 108  Port dedicated to communicate a pointer to a BCCGNS object (Boundary Conditions Info for CGNS export),
* - <B>M_POINT2        </B>= 120  Port dedicated to communicate the coordinate of a point field [array<double,3>].,
* - <B>M_VALUED2       </B>= 130  Port dedicated to communicate a scalar value [double].,
* - <B>M_VALUEB2       </B>= 140  Port dedicated to communicate a scalar value [bool].,
* - <B>M_VALUEB3       </B>= 141  Port dedicated to communicate a scalar value [bool].,
* - <B>M_VALUEB4       </B>= 142  Port dedicated to communicate a scalar value [bool].,
* - <B>M_VALUEB5       </B>= 143  Port dedicated to communicate a scalar value [bool].,
* - <B>M_VALUEI2       </B>= 150  Port dedicated to communicate a scalar value [int].,
* - <B>M_POLYDATA_     </B>= 1100 Port dedicated to communicate a pointer to a vtk polydata mesh [vtkPolyData *].
*/
namespace pin{


/*!
*\enum ConnectionType
*\brief Type of allowed connections of the object: bidirectional, only input or only output.
*/
enum class ConnectionType{
    BOTH 		/**<Bidirectional object. It allows both input and output connections.*/,
    BACKWARD 	/**<Uni-directional backward object. It allows only input connections.*/,
    FORWARD 	/**<Uni-directional forwadr object. It allows only output connections.*/
};

typedef	short int	PortID; /**< mimmo custom definition */


/*! 
 *  \enum containerTAG
 *  \brief Specification of possible containers.
 * 
 * Container TAG are used to define wich kind of data a port is able to
 * communicate/receive. A data type is formed by a containerTAG and a dataTAG.
 * The data type is used to check the compatibility between ports.
 *
 */
enum class containerTAG{

    SCALAR			/**< TAG related to single value container.*/,
    VECTOR			/**< TAG related to std::vector< . > container.*/,
    ARRAY3			/**< TAG related to std::array< . , 3 > container.*/,
    ARRAY4			/**< TAG related to std::array< . , 4 > container.*/,
    VECVEC			/**< TAG related to std::vector< std::vector< . > > container.*/,
    VECARR3         /**< TAG related to std::vector< std::array< . , 3 > > container.*/,
    VECARR2			/**< TAG related to std::vector< std::array< . , 2 > > container */,	
    VECVECARR2		/**< TAG related to std::vector< std::vector< std::array< . , 2 > > > container */,	
    ARR3ARR3		/**< TAG related to std::array< std::array< . , 3 > , 3 > container.*/,
    ARR3ARR3VEC		/**< TAG related to std::array< std::array< std::array< . , 3 > , 3 > , 3 > container.*/,
    MAP				/**< TAG related to std::map< . , . > container.*/,
    MPVECTOR        /**< TAG related to mimmo::MimmoPiercedVector< . > container.*/,
    MPVECARR3       /**< TAG related to mimmo::MimmoPiercedVector< std::array< . , 3 > > container.*/

};

/*! 
 *  \enum dataTAG
 *  \brief Specification of possible data allowed.
 * 
 * Data TAG are used to define wich kind of data a port is able to
 * communicate/receive. A data type is formed by a containerTAG and a dataTAG.
 * The data type is used to check the compatibility between ports.
 *
 */
enum class dataTAG{

    MIMMO_						/**< TAG related to a mimmo::MimmoObject* data.*/,
    FILEINFODATA    			/**< TAG related to mimmo::FileInfoData data.*/,
    INT							/**< TAG related to a int data.*/,
    SHORT						/**< TAG related to a short (int) data.*/,
    LONG						/**< TAG related to a long (int) data.*/,
    FLOAT						/**< TAG related to a double data.*/,
    BOOL						/**< TAG related to a bool data.*/,
    STRING						/**< TAG related to a string data.*/,
    STRING_						/**< TAG related to a string pointer data.*/,
    MPVECFLOAT                  /**< TAG related to a std::vector of bitpit::PiercedVector<double> data.*/,
    MPVECARR3FLOAT              /**< TAG related to a std::vector of bitpit::PiercedVector<std::array<double,3> > data.*/,
    PVECFLOAT_                  /**< TAG related to a bitpit::PiercedVector<double>* data pointer.*/,
    LONGSHORT					/**< TAG related to a couple (normally used in pair or map containers) of a long and a short element data*/,
    SHAPET						/**< TAG related to a mimmo::ShapeType data.*/,
    SHAPE_						/**< TAG related to a mimmo::BasicShape* data.*/,
    COORDT						/**< TAG related to a mimmo::CoordType data.*/,
    POLYDATA_					/**< TAG related to a VTK::vtkPolyData* data.*/,
    TRACKINGPTR_                /**< TAG related to a generic object derived from mimmo::TrackingPointer class */,
    BCCGNS_                     /**< TAG related to a mimmo::BCCGNS* object (class with Boundary Conditions Info for CGNS export) */
};


bool addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced = false);

void removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec);

void removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

bool checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);


} //end of pin namespace

//log variables
extern std::string MIMMO_LOG_FILE; /**<Name of logger file.*/

void    setLogger(std::string log);

void    warningXML(bitpit::Logger* log, std::string name);

//expert variable
extern bool MIMMO_EXPERT; /**<Flag that defines expert mode (true) or safe mode (false).
                                In case of expert mode active the mandatory ports are not checked. */

void setExpertMode(bool flag = true);

//miscellanea
double  maxvalmp(const MimmoPiercedVector<double, long int> & field);

}//end namespace mimmo

#endif
