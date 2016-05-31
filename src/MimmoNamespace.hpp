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
	COORDS	/*! Port dedicated to communicates coordinates of points.*/,
	DISPLS	/*! Port dedicated to communicates displacements of points.*/,
	GDISPLS	/*! Port dedicated to communicates displacements of geometry vertex.*/,
	FILTER	/*! Port dedicated to communicates a scalar field used as filter function.*/,
	BMATRIX	/*! Port dedicated to communicates an array of array 3x3 of double (Bending).*/,
	BCOEFFS	/*! Port dedicated to communicates an array of array 3x3 of vector of double (Bending).*/,
	GLOBAL,
	DEG

};
/**< Port type specification.
 *
 * A type of data is related to each label. Same type of data can be related to
 * multiple type of ports but with different meaning.
 *
 * mimmo::pin::PortType::COORDS -
 * Port dedicated to communicates coordinates of points.
 * A port COORDS communicates a std::vector<std::array<double, 3> >.
 *
 * mimmo::pin::PortType::DISPLS -
 * Port dedicated to communicates displacements of points.
 * A port DISPLS communicates a std::vector<std::array<double, 3> >.
 *
 *  mimmo::pin::PortType::FILTER -
 *  Port dedicated to communicates a scalar field used as filter function.
 *  A port FILTER communicates a std::vector<double>.
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
