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
#include "bitpit.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;


// =================================================================================== //


void test00001() {

    /* Creation of Generic output block to read a set of
     * coordinates of cloud points as generic input in csv format.
     */
    GenericInput * read = new GenericInput(true, true);
    read->setReadFromFile(true);
    read->setReadDir("input");
    read->setFilename("generic_input_00001.csv");
    read->setCSV(true);

    /* Creation of Generic output block to write a set of
     * coordinates of cloud points as generic output in csv format.
     */
    GenericOutput * write = new GenericOutput();
    write->setFilename("generic_output_00001.csv");
    write->setCSV(true);

    /* Setup pin connections.
     */
    addPin(read, write, M_COORDS, M_COORDS);

    /* Setup execution chain.
     */
    Chain ch0;
    ch0.addObject(read);
    ch0.addObject(write);

    /* Execution of chain.
     * Use debug flag true to print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
    delete read;
    delete write;

    read = NULL;
    write = NULL;

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

        test00001();

#if ENABLE_MPI==1
    }

    MPI::Finalize();
#endif

    return 0;
}

