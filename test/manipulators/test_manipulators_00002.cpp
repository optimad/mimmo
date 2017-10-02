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

#include "mimmo_manipulators.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;



// =================================================================================== //
/*!
 * Testing RBF manipulator 
 */

int test2() {
    
    //create a mimmoobject containing a single triangle.
    MimmoObject * mesh = new MimmoObject(1);
    dvecarr3E p(3, {{0.0,0.0,0.0}});
    livector1D conn(3,0);
    
    p[1][0] = 1.0;
    p[2][1] = 1.0;
    conn[1] = 1;
    conn[2] = 2;
    
    int counter = 0;
    for(auto &val: p){
        mesh->addVertex(val, counter);
        ++counter;
    }
    mesh->addConnectedCell(conn, bitpit::ElementType::TRIANGLE, 0, 0);
//     mesh->getPatch()->write("undeformed");
    
    double area     = (static_cast<SurfaceKernel * >(mesh->getPatch()))->evalCellArea(0);
    
    
    dvecarr3E rbfpoints, rbfdispls;
    rbfpoints.push_back(p[1]);
    rbfpoints.push_back(p[2]);
    rbfdispls.push_back({{1.0,0.0,0.0}});
    rbfdispls.push_back({{0.0,1.0,0.0}});
    
    MRBF * mrbf = new MRBF();
    mrbf->setGeometry(mesh);
    mrbf->setNode(rbfpoints);
    mrbf->setDisplacements(rbfdispls);
    mrbf->setSupportRadiusValue(0.3);
    mrbf->exec();
    
    
    Apply * applier = new Apply();
    applier->setGeometry(mesh);
    applier->setInput(mrbf->getDisplacements());
    
    applier->exec();

    //recover normal of the triangle rotated, and area of the scaled;
    double area2     = (static_cast<SurfaceKernel * >(mesh->getPatch()))->evalCellArea(0);
    
//     mesh->getPatch()->write("deformed");
    //check phase
    bool check = ( (std::abs(area/area2) - 0.25) <= 1.e-18);
    
    delete mesh;
    delete mrbf;
    delete applier;
    
    std::cout<<"test passed: "<<check<<std::endl;
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
