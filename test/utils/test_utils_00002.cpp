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

#include "mimmo_utils.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;


/*!
 * Creating surface triangular mesh and return it in a MimmoObject.
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \return true if successfully created mesh
 */
bool createMimmoMesh(MimmoObject * mesh){
    
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
    
    bool check = (mesh->getNCells() == 48) && (mesh->getNVertex() == 35);
    
    mesh->buildAdjacencies();
    mesh->getPatch()->write("support");
    return check;
}
// =================================================================================== //
/*!
 * Test: testing CreateSeedsOnSurface utility 
 */
int test2() {
	
    MimmoObject * m1 = new MimmoObject();
    if(!createMimmoMesh(m1)) {
        delete m1;
        return 1;
    }    
    darray3E normal = {{0,0,1}};
    
    CreateSeedsOnSurface * cseed = new CreateSeedsOnSurface();
    cseed->setGeometry(m1);
    cseed->setSeed({{0.0,0.0,0.0}});
    cseed->setNPoints(2);
    cseed->setEngine(0);
    cseed->exec();
    cseed->setPlotInExecution(true);
    
    bool check = std::abs(dotProduct(cseed->getPoints()[1], normal)) < 1.e-18;
    
//     std::cout<<cseed->getPoints()[1]<<std::endl;
    
    cseed->setEngine(1);
    cseed->exec();
    
    darray3E target = {{1.5,1.0,0.0}};
    
//     std::cout<<cseed->getPoints()[1]<<std::endl;
    
    check = check && (norm2(cseed->getPoints()[1]- target) < 1.e-18);
    
    delete m1;
    delete cseed;
    std::cout<<"test passed :" <<check<<std::endl;
    return int(!check);
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
