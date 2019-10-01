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

#include "mimmo_iocgns.hpp"
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;



// =================================================================================== //

int test1() {

    IOCGNS * cgnsI = new IOCGNS();
    cgnsI->setMode(IOCGNS::IOCGNS_Mode::READ);
    cgnsI->setDir("geodata");
    cgnsI->setFilename("grid");

    cgnsI->execute();

    // std::cout<<cgnsI->getGeometry()->getPatch()->getVertexCount()<<std::endl;
    // std::cout<<cgnsI->getSurfaceBoundary()->getPatch()->getVertexCount()<<std::endl;
    //
    // std::cout<<cgnsI->getGeometry()->getPatch()->getCellCount()<<std::endl;
    // std::cout<<cgnsI->getSurfaceBoundary()->getPatch()->getCellCount()<<std::endl;

    bool check = true;
    check = check && ( cgnsI->getGeometry()->getPatch()->getVertexCount()== 201306);
    check = check && ( cgnsI->getSurfaceBoundary()->getPatch()->getVertexCount()== 11320);
    check = check && ( cgnsI->getGeometry()->getPatch()->getCellCount()== 643873);
    check = check && ( cgnsI->getSurfaceBoundary()->getPatch()->getCellCount()== 18856);

    delete cgnsI;

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

        int val = test1() ;

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
