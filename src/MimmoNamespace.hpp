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

namespace mimmo{

class BaseManipulation;

namespace pin{

/*! Type of allowed conections of the object: bidirectional, only input or only output.*/
enum ConnectionType{
	BOTH 		/**<Bidirectional object. It allows both input and output connections.*/,
	BACKWARD 	/**<Uni-directional backward object. It allows only input connections.*/,
	FORWARD 	/**<Uni-directional forwadr object. It allows only output connections.*/
};

/*! Port TAG specification.
*
* A type of data is related to each TAG. Same type of data can be related to
* multiple TAG of ports but with different meaning.
* A port with a TAG can communicate with other TAG in function of its pre-coded
* compatibility.
* The ports defined in the base library MiMMO must have prefix M_ before their name.
* 0-9 coordinates of 3D points
* 10-19 structures for objects handling deformation
* 20-29 structures for objects building structured lattice
* 30-39	matrix structures
* 40-49 FFDlattice setting parameters
* 
* 80-89 fields & geoemtries which refers to passed as std::pair
* 99	MimmoObject pointer structures;
*
*/
enum PortType{

	M_COORDS 			= 0 	/**< Port dedicated to communicate coordinates of points [vector<array<double,3>>].*/,
	M_GLOBAL 			= 1		/**< Port dedicated to communicate coordinates of points in a global reference system [vector<array<double,3>>].*/,
	M_LOCAL 			= 2		/**< Port dedicated to communicate coordinates of points in a local reference system. [vector<array<double,3>>].*/,
	M_DISPLS 			= 10	/**< Port dedicated to communicate displacements of points [vector<array<double,3>>].*/,
	M_GDISPLS 			= 11	/**< Port dedicated to communicate displacements of geometry vertex [vector<array<double,3>>].*/,
	M_FILTER			= 12	/**< Port dedicated to communicate a scalar field used as filter function [vector<double>].*/,
	M_SCALARFIELD		= 19	/**< Port dedicated to communicate a generic scalar field [vector<double>].*/,
	M_POINT				= 20	/**< Port dedicated to communicate the coordinate of a point field [array<double,3>].*/,
	M_AXIS				= 21	/**< Port dedicated to communicate the direction of an axis [array<double,3>].*/,
	M_AXES				= 22	/**< Port dedicated to communicate a reference system [array<array<double,3>,3>].*/,
	M_SPAN				= 23	/**< Port dedicated to communicate the span of elemental object [array<double,3>].*/,
	M_DIMENSION			= 24	/**< Port dedicated to communicate the dimension of lattice mesh [array<int,3>].*/,
	M_INFLIMITS			= 25	/**< Port dedicated to communicate a inferior limits for building elemental shape [array<double,3>].*/,
	M_SHAPE	 			= 26 	/**< Port dedicated to communicate a type of an elemental shape itself [mimmo::ShapeType].*/,
	M_COPYSHAPE			= 27 	/**< Port dedicated to communicate a an elemental shape itself, passed by pointer, instantiated elsewhere [BasicShape *].*/,
	M_SHAPEI	 		= 28 	/**< Port dedicated to communicate a type of an elemental shape itself [int].*/,
	M_VALUED			= 30	/**< Port dedicated to communicate a scalar value [double].*/,
	M_VALUEI			= 31	/**< Port dedicated to communicate a scalar value [int].*/,
	M_VALUEB			= 32	/**< Port dedicated to communicate a condition value [bool].*/,
	M_BMATRIX			= 33	/**< Port dedicated to communicate a matrix of values [array<array<double,3>,3>].*/,
	M_BCOEFFS			= 34	/**< Port dedicated to communicate a matrix of vector of values [array<array<vector<double>,3>,3>].*/,
	M_DEG				= 40	/**< Port dedicated to communicate a set of number of degrees of freedom [array<int,3>] .*/,
	M_RANGE				= 41	/**< Port dedicated to communicate a set of span interval [array<array<double,2>,3>] .*/,
	M_BOOLS3			= 42	/**< Port dedicated to communicate a set of conditions [array<bool,3>].*/,
	M_NURBSCOORDTYPE	= 43	/**< Port dedicated to communicate condition to design NURBS on FFDLattice [array<mimmo::CoordType,3>] .*/,
	M_NURBSWEIGHTS 		= 44 	/**< Port dedicated to communicate condition to exchange NURBS weights on FFDLattice [vector<double>]*/,
	M_FILENAME			= 50	/**< Port dedicated to communicate a filename [string].*/,
	M_DIR				= 51	/**< Port dedicated to communicate a directory path [string].*/,
	M_PAIRVECFIELD 		= 80 	/**< Port dedicated to communicate a vector field on a MimmoObject geometry as std::pair<MimmoObject *, vector<array<double,3>>> */,
	M_PAIRSCAFIELD 		= 81 	/**< Port dedicated to communicate a scalar field on a MimmoObject geometry as std::pair<MimmoObject *, vector<double>> */,
	M_GEOM 				= 99	/**< Port dedicated to communicate a pointer to a geometry [MimmoObject *].*/,

	M_POINT2			= 120	/**< Port dedicated to communicate the coordinate of a point field [array<double,3>].*/,
	M_VALUED2			= 130	/**< Port dedicated to communicate a scalar value [double].*/,

};

typedef	short int	PortID;

bool addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced = false);

void removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec);

void removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

bool checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

};

};//end namespace mimmo

//#include "MimmoNamespace.tpp"

#endif
