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

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <unordered_map>
#include "Operators.hpp"
#include "customOperators.hpp"
#include "bitpit_IO.hpp" 

namespace mimmo{

/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 *	\authors		Arpa Rocco
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
 *	\authors		Arpa Rocco
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


/*!
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 *	\authors		Arpa Rocco
 *
 *	\brief Class for VTK Output of unstructured grids  
 *
 *	Class for managing output on file of unstructured grids with scalar/vector data attached, based on VTK's formats (*.vtu)
 */

class VTK_DATAMESH: public bitpit::VTKUnstructuredGrid
{
protected:
	dvecarr3E * m_points; /*!< pointer to list of vertices */
	ivector2D * m_connectivity; /*!< pointer to cell-vertex connectivity */
	std::unordered_map< std::string, dvector1D * > m_scalar; /*! map of scalar field added */
	std::unordered_map< std::string, dvecarr3E * > m_vector; /*! map of vector field added */
	bitpit::VTKElementType m_elDM;
public:
	VTK_DATAMESH();
	VTK_DATAMESH(dvecarr3E *, ivector2D * ,std::string dir_, std::string name_, bitpit::VTKFormat cod_);
	VTK_DATAMESH(dvecarr3E *, ivector2D * ,std::string filename);
	~VTK_DATAMESH();
	
	bool linkMeshData(dvecarr3E *, ivector2D * );
	void unlinkMeshData();
	bool addScalarField(std::string &, dvector1D &, bitpit::VTKLocation);
	bool addVectorField(std::string &, dvecarr3E &, bitpit::VTKLocation);
	void removeField(std::string &);
	void removeAllFields();
	void flushData(  std::fstream &str, bitpit::VTKFormat codex_, std::string name  ) ; //CRTP
	
private:
	void flushScalarField(std::string,bitpit::VTKFormat codex_, std::fstream &);
	void flushVectorField(std::string,bitpit::VTKFormat codex_, std::fstream &);
};



}
#endif // __MiMMO_VTKINTERFACES__
