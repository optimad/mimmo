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
#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__

#include "BasicMeshes.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\date			24/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Structured lattice.
 *
 *	Basically, it builds an elemental 3D shape
 *  (box, sphere, cylinder or part of them) around the geometry and set a structured cartesian mesh of control
 *  points on it (lattice). NO displacements for control points and NO NURBS parameters for FFD are present
 *  in this structure, only geometrical information are stored in the object.
 *
 */
class Lattice: public BaseManipulation, public UStructMesh {

protected:
	double		m_np;			/**< Number of control nodes.*/
	ivector1D   m_intMapDOF;    /**< Map of grid nodes -> degrees of freedom of lattice */
	
public:
	Lattice();
	virtual ~Lattice();

	//copy operators/constructors
	Lattice(const Lattice & other);
	Lattice & operator=(const Lattice & other);

	//clean structure;
	void 		clearLattice();

	//internal methods
	int			getNNodes();
	dvecarr3E 	getGlobalCoords();
	dvecarr3E 	getLocalCoords();

	int 		accessDOFFromGrid(int index);
	int 		accessGridFromDOF(int index);

	//plotting wrappers
	void		plotGrid(std::string directory, std::string filename, int counter, bool binary);
	void		plotCloud(std::string directory, std::string filename, int counter, bool binary);

	//building method
	virtual void build();
	void execute();
	
	//nodal displacement utility	
protected:
	void 		resizeMapDof();
private:
	int			reduceDimToDOF(int,int,int, bvector1D &info);

};

}

#endif /* __LATTICE_HPP__ */
