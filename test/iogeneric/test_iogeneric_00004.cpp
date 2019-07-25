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
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;



// =================================================================================== //
/*!
 * Reading a file ascii vtu with MimmoGeometry
 */
int test1() {

    MimmoGeometry * reader = new MimmoGeometry();
    reader->setIOMode(IOMode::READ);
    reader->setReadDir("geodata");
    reader->setReadFilename("curve");
    reader->setReadFileType(FileType::CURVEVTU);
    reader->setCodex(0);
    reader->exec();

    bool check = reader->getGeometry()->getPatch()->getCellCount() == 412;
    check = check && reader->getGeometry()->getPatch()->getVertexCount() == 414;
    std::cout<<"test1 passed :"<<check<<std::endl;

    MimmoGeometry * writer = new MimmoGeometry();
    writer->setGeometry(reader->getGeometry());
    writer->setIOMode(IOMode::WRITE);
    writer->setWriteDir(".");
    writer->setWriteFilename("curve_copy");
    writer->setWriteFileType(FileType::CURVEVTU);
    writer->setCodex(0);
    writer->exec();

    delete reader;
    delete writer;
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
            std::cout<<"test_iogeneric_00004 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
