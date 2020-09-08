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
 \ *---------------------------------------------------------------------------*/

#include "mimmo_core.hpp"
#include <bitpit_common.hpp>

/*
 * Test 00003
 * Testing BasicShapes and skdTree searching.
 * extract sub-patch contained in a Cylinder
 */


/*!
 * Creating surface triangular mesh and return it in a MimmoObject.
 * Pidding the M of mimmo
 *
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \return true if successfully created mesh
 */
bool createMimmoMesh(mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh){

	double dx = 0.25, dy = 0.25;
	int nV, nC;
	//create vertexlist
	dvecarr3E vertex(35,{{0.0,0.0,0.0}});
	livector2D conn(48, livector1D(3));

	for(int i=0; i<7; ++i){
		for(int j=0; j<5; j++){
			nV = 5*i + j;
			vertex[nV][0] = i*dx;
			vertex[nV][1] = j*dy;
		}
	}

	for(int j=0; j<4; ++j){
		for(int i=0; i<3; ++i){
			nC = 8*i + 2*j;

			conn[nC][0] = 5*i + j;
			conn[nC][1] = 5*(i+1) + j;
			conn[nC][2] = 5*i + j+1;

			conn[nC+1][0] = 5*(i+1) + j;
			conn[nC+1][1] = 5*(i+1) + j+1;
			conn[nC+1][2] = 5*i + j+1;
		}
	}

	for(int j=0; j<4; ++j){
		for(int i=3; i<6; ++i){
			nC = 8*i + 2*j;

			conn[nC][0] = 5*i + j;
			conn[nC][1] = 5*(i+1) + j;
			conn[nC][2] = 5*(i+1) + j+1;

			conn[nC+1][0] = 5*i + j;
			conn[nC+1][1] = 5*(i+1) + j+1;
			conn[nC+1][2] = 5*i + j+1;
		}
	}

	mesh->getVertices().reserve(35);
	mesh->getCells().reserve(48);

	//fill the mimmoObject;
	long cV=0;
	for(auto & val: vertex){
		mesh->addVertex(val, cV);
		cV++;
	}

	long cC=0;
	bitpit::ElementType eltype = bitpit::ElementType::TRIANGLE;
	for(auto & val: conn){
		mesh->addConnectedCell(val, eltype, cC);
		cC++;
	}

	bool check = (mesh->getNCells() == 48) && (mesh->getNVertices() == 35);

	mesh->updateAdjacencies();
	return check;
}

/*!
 * Extract a submesh from a target MimmoObject mesh.
 *
 * \param[in] mesh pointer to a MimmoObject target mesh .
 * \param[in] list of cells to extract for target mesh.
 * \return unique_ptr to new mimmoobject mesh
 */
mimmo::MimmoSharedPointer<mimmo::MimmoObject> createSubMesh(mimmo::MimmoSharedPointer<mimmo::MimmoObject> original, livector1D & list){

	mimmo::MimmoSharedPointer<mimmo::MimmoObject> result(new mimmo::MimmoObject(original->getType()));

	//fill the mimmoObject;
	long pid;
	bitpit::ElementType eltype = bitpit::ElementType::TRIANGLE;

	livector1D vertlist = original->getVertexFromCellList(list);

	for(auto & val: vertlist){
		result->addVertex(original->getVertexCoords(val), val);
	}

	for(auto & val: list){
		pid = original->getCells()[val].getPID();
		result->addConnectedCell(original->getCellConnectivity(val), eltype, pid, val);
	}

	result->updateAdjacencies();

	return result;
}


// =================================================================================== //

int test3() {

    mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh(new mimmo::MimmoObject());
	livector1D list;

	//create the target mesh;
	bool check = createMimmoMesh(mesh);

	if(!check){
		std::cout<<"ERROR.Not able to create MimmoObject mesh"<<std::endl;
        return 1;
	}else{
		mesh->getPatch()->write("original_t3");
		std::cout<<"Target Mesh written to file original_t3.vtu"<<std::endl;
	}

	mimmo::Cylinder * shape = new mimmo::Cylinder();
	shape->setOrigin({{-0.00001, 0.5,0.0}});
	shape->setSpan({{0.51, BITPIT_PI, 1.2}});
	shape->setInfLimits({{0.0,-0.5*BITPIT_PI, 0.0}});
	shape->setRefSystem(2, {{0.0,1.0,0.0}});

	list.reserve(16);
	for(int i=0; i<16; ++i)	list.push_back(i);

	//extract a sub mesh from target.
	auto listex = shape->includeGeometry(mesh);
	check = (listex.size() == list.size());

	if(check){
		std::sort(listex.begin(), listex.end());
		int counter = 0;
		for(auto &val : list){
			check = check && (val == listex[counter]);
			++counter;
		}
	}

	mimmo::MimmoSharedPointer<mimmo::MimmoObject> subPatch;
	if(!check){
		std::cout<<"ERROR.Not able to extract sub-patch included in the cylindrical shape w/ bvTree"<<std::endl;
        delete shape;
        return 1;
	}else{
		subPatch = std::move(createSubMesh(mesh, listex));
		subPatch->getPatch()->write("subpatch_t3");
		std::cout<<"Sub-Patch included in cylinder extracted w/ bvTree and written to file subpatch_t3.vtu"<<std::endl;
	}

	delete shape;

	return 0;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
        MPI_Init(&argc, &argv);
#endif

    int val = 1;

    /**<Calling mimmo Test routines*/
    try{
        val = test3() ;
    }
    catch(std::exception & e){
        std::cout<<"test_core_00003 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
