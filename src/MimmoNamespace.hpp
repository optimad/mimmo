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

enum PortsType{BOTH, BACKWARD, FORWARD}; 	/**< Type of pins of the object: bidirectional,
												only input or only output.*/

enum PortType{
	COORDS 	= 0		/*! Port dedicated to communicate coordinates of points [vector<array<double,3>>].*/,
	GLOBAL 	= 1		/*! Port dedicated to communicate coordinates of points in a global reference system [vector<array<double,3>>].*/,
	LOCAL 	= 2		/*! Port dedicated to communicate coordinates of points in a local reference system. [vector<array<double,3>>].*/,
	DISPLS 	= 10	/*! Port dedicated to communicate displacements of points [vector<array<double,3>>].*/,
	GDISPLS = 11	/*! Port dedicated to communicate displacements of geometry vertex [vector<array<double,3>>].*/,
	FILTER	= 12	/*! Port dedicated to communicate a scalar field used as filter function [vector<double>].*/,
	POINT	= 20	/*! Port dedicated to communicate the coordinate of a point field [array<double,3>].*/,
	AXIS	= 21	/*! Port dedicated to communicate the direction of an axis [array<double,3>].*/,
	AXES	= 22	/*! Port dedicated to communicate a reference system [array<array<double,3>,3>].*/,
	VALUE	= 30	/*! Port dedicated to communicate a scalar value [double].*/,
	BMATRIX	= 31	/*! Port dedicated to communicate a matrix of values [array<array<double,3>,3>].*/,
	BCOEFFS	= 32	/*! Port dedicated to communicate a matrix of vector of values [array<array<vector<double>,3>,3>].*/,
	DEG		= 40	/*! Port dedicated to communicate a set of number of degrees of freedom [array<double,3>] .*/,
	RANGE	= 41	/*! Port dedicated to communicate a set of span interval [array<double,3>] .*/,
	BOOLS3	= 42	/*! Port dedicated to communicate a set of conditions [array<boolean,3>].*/,
	GEOM 	= 99	/*! Port dedicated to communicate a pointer to a geometry [MimmoObject*].*/,

};
/**< Port TAG specification.
 *
 * A type of data is related to each TAG. Same type of data can be related to
 * multiple TAG of ports but with different meaning.
 * A port with a TAG can communicate with other TAG in function of its pre-coded
 * compatibility.
 *
 */

typedef	short int	PortID;

bool addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced = false);

bool addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortType portS, PortType portR, bool forced = false);

void removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec);

void removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

void removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortType portS, PortType portR);

bool checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

bool checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortType portS, PortType portR);

};

};//end namespace mimmo

//#include "MimmoNamespace.tpp"

#endif
