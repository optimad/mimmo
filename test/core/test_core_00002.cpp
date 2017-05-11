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
using namespace std;
using namespace bitpit;
using namespace mimmo;

/*
 * Test 00002
 * Testing base Data Structure: MimmoObject 
 */


/*!
 * Creating surface triangular mesh and return it in a MimmoObject.
 * Pidding the M of mimmo
 * 
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \param[in,out] list of cells pidded w 1 to create M.
 * \return true if successfully created mesh
 */
bool createMimmoMesh(MimmoObject * mesh, livector1D & list){
	
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
	
	//fill the list of marked simplicies;
	list.reserve(20);
	for(int k=0; k<8; ++k)	{
		list.push_back(8+k);
		list.push_back(32+k);
	}
	
	list.push_back(21);
	list.push_back(22);
	list.push_back(29);
	list.push_back(30);
	
	mesh->getVertices().reserve(35);
	mesh->getCells().reserve(48);
	
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
	
	for(auto & val: list){
		mesh->setPIDCell(val, 1);
	}
	
	
	bool check = (mesh->getNCells() == 48) && (mesh->getNVertex() == 35);
	
	mesh->buildAdjacencies();
	return check;
}

/*!
 * Extract a submesh from a target MimmoObject mesh.
 * 
 * \param[in] mesh pointer to a MimmoObject target mesh .
 * \param[in] list of cells to extract for target mesh.
 * \return unique_ptr to new mimmoobject mesh
 */
std::unique_ptr<MimmoObject> createSubMesh(MimmoObject * original, livector1D & list){
	
	std::unique_ptr<MimmoObject> result(new MimmoObject(original->getType()));
	
	//fill the mimmoObject;
	short pid;
	bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
	
	livector1D vertlist = original->getVertexFromCellList(list);
	
	for(auto & val: vertlist){
		result->addVertex(original->getVertexCoords(val), val);
	}
	
	for(auto & val: list){
		pid = original->getCells()[val].getPID();
		result->addConnectedCell(original->getCellConnectivity(val), eltype, pid, val);
	}
	
	result->buildAdjacencies();
	
	return std::move(result);
}




// =================================================================================== //

int test2() {
	
	MimmoObject * mesh = new MimmoObject();
	livector1D list;
	
	//create the target mesh;
	bool check = createMimmoMesh(mesh, list);
	
	if(!check){
		std::cout<<"ERROR.Not able to create MimmoObject mesh"<<std::endl;
		return 1;
	}else{
		mesh->getPatch()->write("original_t2");
		std::cout<<"Target Mesh written to file original_t2.vtu"<<std::endl;
	}
	
	//extract a sub mesh from target.
	auto listpid = mesh->extractPIDCells(1);
	
	check = (listpid.size() == list.size());
	
	if(check){
		std::sort(list.begin(), list.end());
		std::sort(listpid.begin(), listpid.end());
		
		int counter = 0;
		for(auto &val : list){
			check = check && (val == listpid[counter]);
			++counter;
		}
	}
	
	std::unique_ptr<MimmoObject> subPatch;
	if(!check){
		std::cout<<"ERROR.Not able to extract sub-patch of MimmoObject mesh"<<std::endl;
		return 1;
	}else{
		subPatch = std::move(createSubMesh(mesh, listpid));
		subPatch->getPatch()->write("subpatch_t2");
		std::cout<<"Sub-Patch extracted from target Mesh written to file subpatch_t2.vtu"<<std::endl;
	}
	
	livector1D listBVsp;
	listBVsp.reserve(22);
	for(int i=0; i<10; ++i){
		listBVsp.push_back(5+i);
		listBVsp.push_back(20+i);
	}
	listBVsp.push_back(17);
	listBVsp.push_back(18);
	
	livector1D listBVextra = subPatch->extractBoundaryVertexID();
	
	check = (listBVextra.size() == listBVsp.size());
	
	if(check){
		std::sort(listBVsp.begin(), listBVsp.end());
		std::sort(listBVextra.begin(), listBVextra.end());
		
		int counter = 0;
		for(auto &val : listBVsp){
			check = check && (val == listBVextra[counter]);
			++counter;
		}
	}

	if(!check){
		std::cout<<"ERROR.Not able to extract boundary of the sub-patch of MimmoObject mesh"<<std::endl;
		return 1;
	}else{
		std::cout<<"Sub-Patch boundary extraction successfully"<<std::endl;
	}

	
	MimmoObject * mesh_2 = new MimmoObject();
    mesh_2->setHARDCopy(mesh);
    mesh_2->getPatch()->write("hardcopy");
	
	
	delete mesh;
    delete mesh_2;
	
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

        int val = test2() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
