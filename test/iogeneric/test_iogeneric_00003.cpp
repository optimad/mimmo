/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
 * Reading a generic input from file with GenericDispls
 */
int test3() {

	mimmo::GenericDispls * ginput = new mimmo::GenericDispls();
    ginput->setReadDir("input");
    ginput->setReadFilename("generic_displ_00001.txt");
    ginput->exec();

    auto result = ginput->getLabels();

    bool check = (result.size() == 2);
    check = check && (result[0] == 12 && result[1]==21);

    std::cout<<"test passed :"<<check<<std::endl;

    delete ginput;
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
            std::cout<<"test_iogeneric_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
