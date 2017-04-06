/*---------------------------------------------------------------------------*\
 * 
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

#include "mimmo_core.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;

/*
 * Test 00005
 * Testing Oriented Bounding Box execution 
 */


/*!
 * Creating surface triangular mesh of a cube given its ordered vertices 
 * and return it in a MimmoObject.
 * 
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \param[in] vertex of ordered vertices z->y->x
 * \return true if successfully created cube mesh 
 */
bool createCube(MimmoObject * mesh, dvecarr3E & vertex){
	
	vertex.push_back(0.25*(vertex[0] + vertex[1] + vertex[2] + vertex[3]));
	vertex.push_back(0.25*(vertex[4] + vertex[5] + vertex[6] + vertex[7]));
	vertex.push_back(0.25*(vertex[0] + vertex[1] + vertex[4] + vertex[5]));
	vertex.push_back(0.25*(vertex[2] + vertex[3] + vertex[7] + vertex[6]));
	vertex.push_back(0.25*(vertex[0] + vertex[2] + vertex[6] + vertex[4]));
	vertex.push_back(0.25*(vertex[5] + vertex[1] + vertex[7] + vertex[3]));
	
	int nC;
	//create vertexlist
	livector2D conn(24, livector1D(3,0));
	
	conn[0]  = {{8,0,1}};
	conn[1]  = {{8,1,3}};
	conn[2]  = {{8,3,2}};
	conn[3]  = {{8,2,0}};
	conn[4]  = {{9,5,4}};
	conn[5]  = {{9,4,6}};
	conn[6]  = {{9,6,7}};
	conn[7]  = {{9,7,5}};
	conn[8]  = {{10,5,1}};
	conn[9]  = {{10,1,0}};
	conn[10] = {{10,0,4}};
	conn[11] = {{10,4,5}};
	conn[12] = {{11,3,7}};
	conn[13] = {{11,7,6}};
	conn[14] = {{11,6,2}};
	conn[15] = {{11,2,3}};
	conn[16] = {{12,2,6}};
	conn[17] = {{12,6,4}};
	conn[18] = {{12,4,0}};
	conn[19] = {{12,0,2}};
	conn[20] = {{13,5,7}};
	conn[21] = {{13,7,3}};
	conn[22] = {{13,3,1}};
	conn[23] = {{13,1,5}};
	
	mesh->getVertices().reserve(14);
	mesh->getCells().reserve(24);
	
	//fill the mimmoObject;
	long cV=0;
	for(auto & val: vertex){
		mesh->addVertex(val, cV);
		cV++;
	}
	
	long cC=0;
	bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
	for(auto & val: conn){
		std::cout<<val<<std::endl;
		mesh->addConnectedCell(val, eltype, cC);
		cC++;
	}
	
	
	bool check = (mesh->getNCells() == 24) && (mesh->getNVertex() == 14);
	
	mesh->buildAdjacencies();
	return check;
}


// =================================================================================== //

int test5() {
	MimmoObject * mesh = new MimmoObject(1);
	Lattice * lattice = new Lattice();
	
	lattice->setShape(ShapeType::CUBE);
	darray3E origin = {{2.0,1.0,0.0}};
	darray3E span= {{1.0,1.0,1.0}};
	iarray3E dim = {{2,2,2}};
	
	lattice->setOrigin(origin);
	lattice->setSpan(span);
	lattice->setDimension(dim);
	
	darray3E axis;
	axis[0] = 0.0;
	axis[1] = 1.0;
	axis[2] = 1.0;
	axis /= norm2(axis);
 	lattice->setRefSystem(1, axis);
	
	
	lattice->execute();
	auto vlist = lattice->getGlobalCoords();
	
	bool check = createCube(mesh,vlist);
		
	OBBox * obb = new OBBox();
	
	obb->setGeometry(mesh);
	obb->execute();
	
	
	check = true;
	
	check = check && (norm2(obb->getOrigin()-origin) <=1.e-15);
	check = check && (norm2(obb->getSpan()-span) <=1.e-15);
	
	obb->plot(".","obbox_t5", 0, 0);
	mesh->getPatch()->write("mesh_t5");
	
	if(!check){
		std::cout<<"Oriented bounding box calculation failed."<<std::endl;
		return 1;
	}else{
		std::cout<<"Oriented bounding box calculation done. See results in obbox_t5.0000.vtu "<<std::endl;
	}	
	
	delete obb;
	delete mesh;
	delete lattice;
	return 0;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
	
#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling mimmo Test routines*/

        int val = test5() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
