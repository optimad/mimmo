/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __MIMMONAMESPACE_HPP__
#define __MIMMONAMESPACE_HPP__

#include "enum.hpp"

BETTER_ENUM(PortType, int,	
			M_COORDS 			= 0, 
			M_GLOBAL 			= 1,
			M_LOCAL 			= 2,
			M_DISPLS 			= 10,
			M_GDISPLS 			= 11,
			M_FILTER			= 12,
			M_VECTORSI			= 17,
			M_VECTORLI			= 18,
			M_SCALARFIELD		= 19,
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
			M_FILENAME			= 50,
			M_DIR				= 51,
			M_FILENAMEPTR		= 52,
			M_DIRPTR			= 53,
			M_PAIRVECFIELD 		= 80,
			M_PAIRSCAFIELD 		= 81,
			M_VIOLATION			= 82,
			M_GEOM2				= 98,
			M_GEOM 				= 99,
			M_VECGEOM			= 100,
			M_MAPGEOM           = 101,
			M_FINFO				= 102,
			M_FINFO2			= 103,
			M_MAPDCELL			= 104,
			M_MAPDVERT			= 105,
			M_UMGEOSFD			= 106,
			M_UMGEOVFD			= 107,
			M_POINT2			= 120,
			M_VALUED2			= 130,
			M_VALUED3			= 131,
			M_VALUEB2			= 140,
			M_VALUEB3			= 141,
			M_VALUEB4			= 142,
			M_VALUEB5			= 143,
			M_VALUEI2			= 150,
			M_VECPAIRSF			= 200,
			M_VECPAIRVF			= 201
   		);



namespace mimmo{

class BaseManipulation;

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief FileDataInfo is a struct to stock data relative to names of external files.
 *
 *	FileDataInfo has three fields:
 *	an integer relative to the type of file stored 
 *  a string reporting the path of the file
 *  a string reporting the name of the file
 * 
 */
struct FileDataInfo{
	int ftype;
	std::string fname;
	std::string fdir;
	
	FileDataInfo();
	virtual ~FileDataInfo();
	FileDataInfo(const FileDataInfo & other);
	FileDataInfo & operator=(const FileDataInfo & other);
};




/*!
 * 
 *@brief	Namespace utilities to create pin connections
 * 
 * Here are collected enums and methods to create connections between blocks throughout established ports.
 * All the available port types are collected in a special enum, istantiated through BETTER_ENUM macro, and 
 * each type reports a specific TAG in order to be retracked in each moment for compatibility checking purposes.
 * 
 * A type of data is related to each TAG. Same type of data can be related to
 * multiple TAG of ports but with different meaning.
 * A port with a TAG can communicate with other TAG in function of its pre-coded
 * compatibility.
 * The ports defined in the base library MiMMO must have prefix M_ before their name.
 
 * The nomenclature used at the moment says:
 * 
 * 0-9 coordinates of 3D points
 * 10-19 structures for objects handling deformation
 * 20-29 structures for objects building structured lattice
 * 30-39	matrix structures
 * 40-49 FFDlattice setting parameters
 * 80-89 fields & geoemtries which refers to passed as std::pair
 * 99-100	MimmoObject pointer structures;
 * > 100 miscellaneous data
 * 
 * 	All contents are:
 * 
 * - M_COORDS        = 0    Port dedicated to communicate coordinates of points [vector<array<double,3>>].,
 * - M_GLOBAL        = 1    Port dedicated to communicate coordinates of points in a global reference system [vector<array<double,3>>].,
 * - M_LOCAL         = 2    Port dedicated to communicate coordinates of points in a local reference system. [vector<array<double,3>>].,
 * - M_DISPLS        = 10   Port dedicated to communicate displacements of points [vector<array<double,3>>].,
 * - M_GDISPLS       = 11   Port dedicated to communicate displacements of geometry vertex [vector<array<double,3>>].,
 * - M_FILTER        = 12   Port dedicated to communicate a scalar field used as filter function [vector<double>].,
 * - M_VECTORSI      = 17   Port dedicated to communicate a generic list of short integers [vector<short int>].,
 * - M_VECTORLI      = 18   Port dedicated to communicate a generic list of long integers [vector<long int>].,
 * - M_SCALARFIELD   = 19   Port dedicated to communicate a generic scalar field [vector<double>].,
 * - M_POINT         = 20   Port dedicated to communicate the coordinate of a point field [array<double,3>].,
 * - M_AXIS          = 21   Port dedicated to communicate the direction of an axis [array<double,3>].,
 * - M_AXES          = 22   Port dedicated to communicate a reference system [array<array<double,3>,3>].,
 * - M_SPAN          = 23   Port dedicated to communicate the span of elemental object [array<double,3>].,
 * - M_DIMENSION     = 24   Port dedicated to communicate the dimension of lattice mesh [array<int,3>].,
 * - M_INFLIMITS     = 25   Port dedicated to communicate a inferior limits for building elemental shape [array<double,3>].,
 * - M_SHAPE         = 26   Port dedicated to communicate a type of an elemental shape itself [mimmo::ShapeType].,
 * - M_COPYSHAPE     = 27   Port dedicated to communicate a an elemental shape itself, passed by pointer, instantiated elsewhere [BasicShape *].,
 * - M_SHAPEI        = 28   Port dedicated to communicate a type of an elemental shape itself [int].,
 * - M_PLANE         = 29   Port dedicated to communicate planes in a*x+b*y+c*z +d = 0 implicit form [array<double,4>].,	
 * - M_VALUED        = 30   Port dedicated to communicate a scalar value [double].,
 * - M_VALUEI        = 31   Port dedicated to communicate a scalar value [int].,
 * - M_VALUEB        = 32   Port dedicated to communicate a condition value [bool].,
 * - M_BMATRIX       = 33   Port dedicated to communicate a matrix of values [array<array<double,3>,3>].,
 * - M_BCOEFFS       = 34   Port dedicated to communicate a matrix of vector of values [array<array<vector<double>,3>,3>].,
 * - M_VALUESI       = 35   Port dedicated to communicate a scalar short int value [short int].,
 * - M_DEG           = 40   Port dedicated to communicate a set of number of degrees of freedom [array<int,3>] .,
 * - M_RANGE         = 41   Port dedicated to communicate a set of span interval [array<array<double,2>,3>] .,
 * - M_BOOLS3        = 42   Port dedicated to communicate a set of conditions [array<bool,3>].,
 * - M_NURBSCOORDTYPE= 43   Port dedicated to communicate condition to design NURBS on FFDLattice [array<mimmo::CoordType,3>] .,
 * - M_NURBSWEIGHTS  = 44   Port dedicated to communicate condition to exchange NURBS weights on FFDLattice [vector<double>],
 * - M_FILENAME      = 50   Port dedicated to communicate a filename [string].,
 * - M_DIR           = 51   Port dedicated to communicate a directory path [string].,
 * - M_PAIRVECFIELD  = 80   Port dedicated to communicate a vector field on a MimmoObject geometry as std::pair<MimmoObject *, vector<array<double,3>>> ,
 * - M_PAIRSCAFIELD  = 81   Port dedicated to communicate a scalar field on a MimmoObject geometry as std::pair<MimmoObject *, vector<double>> ,
 * - M_VIOLATION     = 82   Port dedicated to communicate a double value, with a BaseManipulation object which refers to as std::pair<BaseManipulation *, double>,
 * - M_GEOM2         = 98   Port dedicated to communicate a pointer to a geometry [MimmoObject *],
 * - M_GEOM          = 99   Port dedicated to communicate a pointer to a geometry [MimmoObject *],
 * - M_VECGEOM       = 100  Port dedicated to communicate a list of pointers to geometries std::vector<MimmoObject *>,
 * - M_MAPGEOM       = 101  Port dedicated to communicate a cells map list of the type std::unordered_map<long, short>,
 * - M_FINFO         = 102  Port dedicated to communicate a list of FileInfoData classes std::vector<FileInfoData>,
 * - M_FINFO2        = 103  Port dedicated to communicate a list of FileInfoData classes std::vector<FileInfoData>.,	 
 * - M_MAPDCELL      = 104  Port dedicated to communicate a map of cell-ids of type <long, std::pair<int, long>>,
 * - M_MAPDVERT      = 105  Port dedicated to communicate a map of vertex-ids of type <long, std::pair<int, long>>,
 * - M_UMGEOSFD      = 106  Port dedicated to communicate a map of MimmoObject* and dvector1D*,
 * - M_UMGEOVFD      = 107  Port dedicated to communicate a map of MimmoObject* and dvecarr3E*,
 * - M_POINT2        = 120  Port dedicated to communicate the coordinate of a point field [array<double,3>].,
 * - M_VALUED2       = 130  Port dedicated to communicate a scalar value [double].,
 * - M_VALUEB2       = 140  Port dedicated to communicate a scalar value [bool].,
 * - M_VALUEB3       = 141  Port dedicated to communicate a scalar value [bool].,
 * - M_VALUEB4       = 142  Port dedicated to communicate a scalar value [bool].,
 * - M_VALUEB5       = 143  Port dedicated to communicate a scalar value [bool].,
 * - M_VALUEI2       = 150  Port dedicated to communicate a scalar value [int].
 * - M_VECPAIRSF     = 200  Port dedicated to communicate a std::vector<std::pair<MimmoObject*, dvector1D*> >.
 * - M_VECPAIRVF     = 201  Port dedicated to communicate a std::vector<std::pair<MimmoObject*, dvecarr3E*> >. 
 */
namespace pin{

/*! Type of allowed conections of the object: bidirectional, only input or only output.*/
enum ConnectionType{
	BOTH 		/**<Bidirectional object. It allows both input and output connections.*/,
	BACKWARD 	/**<Uni-directional backward object. It allows only input connections.*/,
	FORWARD 	/**<Uni-directional forwadr object. It allows only output connections.*/
};

typedef	short int	PortID;


/*! Container TAG specification.
*	Container TAG are used to define wich kind of data a port is able to
*	communicate/receive. A data type is formed by a containerTAG and a dataTAG.
*	The data type is used to check the compatibility between ports.
*
*/
enum class containerTAG{

	SCALAR			/**< TAG related to single value container.*/,
	VECTOR			/**< TAG related to std::vector< . > container.*/,
	ARRAY3			/**< TAG related to std::array< . , 3 > container.*/,
	ARRAY4			/**< TAG related to std::array< . , 4 > container.*/,
	VECVEC			/**< TAG related to std::vector< std::vector< . > > container.*/,
	VECARR3			/**< TAG related to std::vector< std::array< . , 3 > > container.*/,
	VECARR2			/**< TAG related to std::vector< std::array< . , 2 > > container */,	
	VECVECARR2		/**< TAG related to std::vector< std::vector< std::array< . , 2 > > > container */,	
	ARR3ARR3		/**< TAG related to std::array< std::array< . , 3 > , 3 > container.*/,
	ARR3ARR3VEC		/**< TAG related to std::array< std::array< std::array< . , 3 > , 3 > , 3 > container.*/,
	MAP				/**< TAG related to std::map< . , . > container.*/,
	UN_MAP			/**< TAG related to std::unordered_map< . , . > container.*/,
	PAIR			/**< TAG related to std::pair< . , . > container.*/

};

/*! Data TAG specification.
*	Data TAG are used to define wich kind of data a port is able to
*	communicate/receive. A data type is formed by a containerTAG and a dataTAG.
*	The data type is used to check the compatibility between ports.
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
	MIMMO_VECFLOAT_				/**< TAG related to a couple (normally used in pair container) of mimmo::MimmoObject* and std::vector<double>* data.*/,
	MIMMO_VECARR3FLOAT_			/**< TAG related to a couple (normally used in pair container) of mimmo::MimmoObject* and std::vector<std::array<double,3> >* data.*/,
	LONGSHORT					/**< TAG related to a couple (normally used in pair or map containers) of a long and a short element data*/,
	STRINGPAIRINTMIMMO_			/**< TAG related to a couple (normally used in pair or map containers) of a std::string and a pair of int and mimmo::MimmoObject* */,
	LONGPAIRINTLONG				/**< TAG related to a couple (normally used in pair or map containers) of a long and a pair of int and long */,
	PAIRMIMMO_VECFLOAT_			/**< TAG related to a std::pair<mimmo::MimmoObject*, std::vector<double>* > data.*/,
	PAIRMIMMO_VECARR3FLOAT_		/**< TAG related to a std::pair<mimmo::MimmoObject*, std::vector<std::array<double,3> >* > data.*/,
	PAIRMIMMO_OBJFLOAT_     	/**< TAG related to a std::pair<mimmo::BaseManipulation *, double> data.*/,
	SHAPET						/**< TAG related to a mimmo::ShapeType data.*/,
	SHAPE_						/**< TAG related to a mimmo::BasicShape* data.*/,
	COORDT						/**< TAG related to a mimmo::CoordType data.*/,
	POLYDATA_					/**< TAG related to a VTK::vtkPolyData* data.*/,
	TRACKINGPTR_				/**< TAG related to a generic object derived from mimmo::TrackingPointer class */
};


bool addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced = false);

void removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec);

void removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

bool checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

};

};//end namespace mimmo

//#include "MimmoNamespace.tpp"

#endif
