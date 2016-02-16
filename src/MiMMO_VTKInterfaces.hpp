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
 \ *---------------------------------------------------------------------------*/

#ifndef __MiMMO_VTKINTERFACES__
#define __MiMMO_VTKINTERFACES__

//std libs
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <map>

// local libs
#include "Operators.hpp"
#include "customOperators.hpp"
#include "bitpit_IO.hpp" //?? what is its new interface in bitpit? please check it out


/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for VTK Output of structured grids  
 *
 *	Class for managing output on file of deformable structured grids as unstructured ones, based on VTK's formats (*.vtu)
 */

class VTK_BASICMESH: public bitpit::VTKUnstructuredGrid
{
protected:
	dvecarr3E * m_points; /*!< pointer to list of vertices */
	ivector2D * m_connectivity; /*!< pointer to cell-vertex connectivity */
public:
	VTK_BASICMESH();
	VTK_BASICMESH(std::string dir_, std::string name_, bitpit::VTKFormat cod_, int nP_, int nC_, int nConn_);
	~VTK_BASICMESH();
	
	void linkData(dvecarr3E &, ivector2D & );
	void unlinkData();
	void flushData(  std::fstream &str, bitpit::VTKFormat codex_, std::string name  ) ; //CRTP
};

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for VTK Output of Cloud points 
 *
 *	Class for managing output on file of point clouds, based on VTK's formats (*.vtu)
 */
class VTK_BASICCLOUD: public bitpit::VTKUnstructuredGrid
{
protected:
	dvecarr3E * m_points; /*!< pointer lo cloud points list */
public:
	VTK_BASICCLOUD();
	VTK_BASICCLOUD(std::string dir_, std::string name_, bitpit::VTKFormat cod_, int nP_);
	~VTK_BASICCLOUD();
	
	void linkData(dvecarr3E &);
	void unlinkData();
	void flushData(  std::fstream &str, bitpit::VTKFormat codex_, std::string name  ) ; //CRTP
};

#endif // __MiMMO_VTKINTERFACES__
