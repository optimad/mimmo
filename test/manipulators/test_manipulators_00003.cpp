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
#include <bitpit_common.hpp>

// =================================================================================== //
/*!
 * Testing Lattice manipulator
 */

int test3() {

    //create a mimmoobject containing a single triangle.
	mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh(new mimmo::MimmoObject(1));
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
    darray3E normal = (static_cast<bitpit::SurfaceKernel * >(mesh->getPatch()))->evalFacetNormal(0);
    darray3E bbmin, bbmax;
    mesh->getBoundingBox(bbmin, bbmax);

    darray3E origin = 0.5*(bbmin + bbmax);
    darray3E span   = bbmax - bbmin;
    for(auto & val : span){
        if (val<1.e-18) val = 2.E-03;
    }
    double zdispl = 0.5*span[0]*std::tan(BITPIT_PI/3.0);

    dvecarr3E displ(8,{{0.0,0.0,0.0}});

    displ[0][2] = -1.0*zdispl;
    displ[1][2] = -1.0*zdispl;
    displ[4][2] = -1.0*zdispl;
    displ[5][2] = -1.0*zdispl;

    displ[2][2] = zdispl;
    displ[3][2] = zdispl;
    displ[6][2] = zdispl;
    displ[7][2] = zdispl;

    mimmo::FFDLattice * latt = new mimmo::FFDLattice();
    latt->setGeometry(mesh);
    latt->setShape(mimmo::ShapeType::CUBE);
    latt->setOrigin(origin);
    latt->setSpan(span);
    latt->setDimension(ivector1D(3,2));
    latt->setDisplacements(displ);
    latt->exec();

    mimmo::Apply * applier = new mimmo::Apply();
    applier->setGeometry(mesh);
    applier->setInput(latt->getDeformation());
    applier->exec();

    //recover normal of the triangle rotated, and area of the scaled;
    darray3E normal2 = (static_cast<bitpit::SurfaceKernel * >(mesh->getPatch()))->evalFacetNormal(0);

    //check phase
    bool check =( (std::abs(std::acos(dotProduct(normal2,normal))) - BITPIT_PI/3.0) <= 1.E-18);

    delete latt;
    delete applier;

    std::cout<<"test passed: "<<check<<std::endl;
    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
		/**<Calling mimmo Test routines*/
        int val =1;
        try{
            val = test3() ;
        }

        catch(std::exception & e){
            std::cout<<"test_manipulators_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
