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
 * Testing Global Manipulators-> Rotation and Scaling 
 */

int test1() {
    
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
    
    //recover normal of the triangle, and area;
    darray3E normal = (static_cast<SurfaceKernel * >(mesh->getPatch()))->evalFacetNormal(0);
    double area     = (static_cast<SurfaceKernel * >(mesh->getPatch()))->evalCellArea(0);
    
    MimmoObject * mesh2 = new MimmoObject(1);
    mesh2->setHARDCopy(mesh);
    
    RotationGeometry * rot = new RotationGeometry();
    rot->setGeometry(mesh);
    rot->setAxis({{0.0,0.0,0.0}},{{1.0,0.0,0.0}});
    rot->setRotation(M_PI/3.);
    rot->exec();
    
    ScaleGeometry * scale = new ScaleGeometry();
    scale->setGeometry(mesh);
    scale->setScaling({{0.5,0.5,0.5}});
    scale->exec();
    
    
    Apply * applier = new Apply();
    applier->setGeometry(mesh);
    applier->setInput(rot->getDisplacements());
    
    Apply * applier2 = new Apply();
    applier2->setGeometry(mesh2);
    applier2->setInput(scale->getDisplacements());
    
    applier->exec();
    applier2->exec();


    //recover normal of the triangle rotated, and area of the scaled;
    darray3E normal2 = (static_cast<SurfaceKernel * >(mesh->getPatch()))->evalFacetNormal(0);
    double area2     = (static_cast<SurfaceKernel * >(mesh2->getPatch()))->evalCellArea(0);
    
    
    //check phase
    bool check = true;
    
    check = check && ( (std::abs(std::acos(dotProduct(normal2,normal))) - M_PI/3.0) <= 1.e-18);
    check = check && ( (std::abs(area/area2) - 4.0) <= 1.e-18);
    
//     mesh->getPatch()->write("rotated");
//     mesh2->getPatch()->write("scaled");
    
    delete mesh;
    delete mesh2;
    delete rot;
    delete scale;
    delete applier;
    delete applier2;
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

        int val = test1() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
