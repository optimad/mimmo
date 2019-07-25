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
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;


// =================================================================================== //
/*!
 * Test: testing OBbox utility
 */
int test3() {

    MimmoGeometry * reader1 = new MimmoGeometry();
    reader1->setIOMode(IOMode::READ);
    reader1->setReadDir("geodata");
    reader1->setReadFilename("stanfordBunny2");
    reader1->setReadFileType(FileType::STL);
    reader1->execute();

    MimmoGeometry * reader2 = new MimmoGeometry();
    reader2->setIOMode(IOMode::READ);
    reader2->setReadDir("geodata");
    reader2->setReadFilename("sphere2");
    reader2->setReadFileType(FileType::STL);
    reader2->execute();

    OBBox * box1 = new OBBox();
    box1->setGeometry(reader1->getGeometry());
    box1->setGeometry(reader2->getGeometry());

    box1->exec();
    box1->plot(".","obbox", 0, false);

    auto span1 = box1->getSpan();

    box1->setForceAABB(true);
    box1->exec();
    box1->plot(".","obbox", 1, false);

    auto span2 = box1->getSpan();

    delete reader1;
    delete reader2;
    delete box1;

    std::cout<<"AABB volume: "<<span2[0]*span2[1]*span2[2]<<std::endl;
    std::cout<<"OBB volume:  "<<span1[0]*span1[1]*span1[2]<<std::endl;

    std::cout<<"test passed "<<std::endl;
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
            std::cout<<"test_utils_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
