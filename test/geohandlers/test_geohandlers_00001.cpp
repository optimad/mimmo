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

#include "mimmo_geohandlers.hpp"

// =================================================================================== //
/*!
 * Testing geohandlers module. Stitching geometries, Splitting and Overlapping fields defined on them
 */
int test1() {

    //define 3 single triangle mesh
	mimmo::MimmoObject * m1 = new mimmo::MimmoObject(1);
	mimmo::MimmoObject * m2 = new mimmo::MimmoObject(1);
	mimmo::MimmoObject * m3 = new mimmo::MimmoObject(1);

    dvecarr3E points(5, {{0.0,0.0,0.0}});
    points[1] = {{1.0,0.0,0.0}};
    points[2] = {{2.0,0.0,0.0}};
    points[3] = {{0.5,1.0,0.0}};
    points[4] = {{1.5,1.0,0.0}};

    livector1D conn(3, 0);

    m1->addVertex(points[0],0);
    m1->addVertex(points[1],1);
    m1->addVertex(points[3],3);
    conn[0] = 0; conn[1] = 1; conn[2] = 3;
    m1->addConnectedCell(conn, bitpit::ElementType::TRIANGLE, long(12), long(4));

    m2->addVertex(points[1],1);
    m2->addVertex(points[3],3);
    m2->addVertex(points[4],4);
    conn[0] = 1; conn[1] = 4; conn[2] = 3;
    m2->addConnectedCell(conn, bitpit::ElementType::TRIANGLE, long(21), long(9));


    m3->addVertex(points[1],1);
    m3->addVertex(points[2],2);
    m3->addVertex(points[4],4);
    conn[0] = 1; conn[1] = 2; conn[2] = 4;
    m3->addConnectedCell(conn, bitpit::ElementType::TRIANGLE, long(44), long(12));


    //stitch geometries
    mimmo::StitchGeometry * stitch1 = new mimmo::StitchGeometry(1);
    stitch1->forceRePID(true);
    mimmo::StitchGeometry * stitch2 = new mimmo::StitchGeometry(1);

    stitch1->addGeometry(m1);
    stitch1->addGeometry(m2);

    stitch2->addGeometry(m2);
    stitch2->addGeometry(m3);

    stitch1->exec();
    stitch2->exec();

    if(stitch1->getGeometry()->getNCells() !=2 || stitch2->getGeometry()->getNCells() !=2){
        delete m1;
        delete m2;
        delete m3;
        delete stitch1;
        delete stitch2;
        return 1;
   }

    bool check = true;

    std::cout<<"test passed :"<<check<<std::endl;
    m1->getPatch()->write("t1");
    m2->getPatch()->write("t2");
    m3->getPatch()->write("t3");

    stitch1->getGeometry()->getPatch()->write("mesh1");
    stitch2->getGeometry()->getPatch()->write("mesh2");

    delete m1;
    delete m2;
    delete m3;
    delete stitch1;
    delete stitch2;

    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
        int val = 1;
		try{
            /**<Calling mimmo Test routines*/
            val = test1() ;
        }
        catch(std::exception & e){
            std::cout<<"test_geohandlers_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
