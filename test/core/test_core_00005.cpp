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
 * \param[in] vertex of ordered vertices x->y->z
 * \return true if successfully created cube mesh 
 */
bool createCube(MimmoObject * mesh, dvecarr3E & vertex){
	
	int nC;
	//create vertexlist
	livector2D conn;
	conn.reserve(12);
	
	conn[0]  = {{0,1,2}};
	conn[1]  = {{0,2,3}};
	conn[2]  = {{4,5,6}};
	conn[3]  = {{4,6,7}};
	conn[4]  = {{0,1,5}};
	conn[5]  = {{0,5,4}};
	conn[6]  = {{3,2,6}};
	conn[7]  = {{3,6,7}};
	conn[8]  = {{1,2,6}};
	conn[9]  = {{1,6,5}};
	conn[10] = {{0,3,7}};
	conn[11] = {{0,7,4}};
		
	mesh->getVertices().reserve(8);
	mesh->getCells().reserve(12);
	
	//fill the mimmoObject;
	long cV=0;
	for(auto & val: vertex){
		mesh->addVertex(val, cV);
		cV++;
	}
	
	long cC=0;
	bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
	for(auto & val: conn){
		mesh->addConnectedCell(val, eltype, cC);
		cC++;
	}
	
	
	bool check = (mesh->getNCells() == 12) && (mesh->getNVertex() == 8);
	
	mesh->buildAdjacencies();
	return check;
}


// =================================================================================== //

int test5() {
	
	
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
