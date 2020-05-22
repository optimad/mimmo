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

#include "mimmo_iogeneric.hpp"

// =================================================================================== //
/*!
 * Reading a file with MimmoGeometry
 */
int test1() {

	mimmo::MimmoGeometry * reader = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    reader->setReadDir("geodata");
    reader->setReadFilename("prism");
    reader->setReadFileType(FileType::STL);
    reader->exec();

    bool check = reader->getGeometry()->getNCells() == 12288;
    check = check && reader->getGeometry()->getNVertices() == 6146;

    mimmo::MimmoGeometry * readerCopy = new mimmo::MimmoGeometry();
    *readerCopy = *reader;
    check = check && (readerCopy->getGeometry() == reader->getGeometry());

    mimmo::MimmoSharedPointer<mimmo::MimmoObject> objHC = reader->getGeometry()->clone();
    check = check && objHC->getNCells() == 12288;
    check = check && objHC->getNVertices() == 6146;

    std::cout<<"test1 passed :"<<check<<std::endl;

    delete reader;
    delete readerCopy;
    return int(!check);
}

/*!
 * Reading mixed type vtu with MimmoGeometry
 */
int test2() {

	mimmo::MimmoGeometry * reader1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    reader1->setReadDir("geodata");
    reader1->setReadFilename("mixedP2D");
    reader1->setReadFileType(FileType::SURFVTU);
    reader1->exec();

    bool check = reader1->getGeometry()->getNCells() == 20;
    check = check && reader1->getGeometry()->getNVertices() == 19;

    reader1->getGeometry()->getPatch()->write("surface");

    mimmo::MimmoGeometry * reader2 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    reader2->setReadDir("geodata");
    reader2->setReadFilename("mixedP3D");
    reader2->setReadFileType(FileType::VOLVTU);
    reader2->exec();

    check = check && (reader2->getGeometry()->getNCells() == 75);
    check = check && (reader2->getGeometry()->getNVertices() == 52);

    reader2->getGeometry()->getPatch()->write("volume");

    delete reader1;
    delete reader2;
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
        val = std::max(val,test2());
    }
    catch(std::exception & e){
        std::cout<<"test_iogeneric_00001 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
