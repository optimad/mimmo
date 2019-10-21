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

/*!
 * \example genericinput_example_00002.cpp
 *
 * \brief Example of reading/writing a generic csv file of MimmoPiercedVector data.

 * Using: GenericInputMPVData, GenericOutputMPVData
 *
 * <b>To run</b>: ./genericinput_example_00002 \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //


void test00002() {

    /* Creation of Generic output block to read a set of
     * coordinates of cloud points as generic input in csv format.
     */
	mimmo::GenericInputMPVData * read = new mimmo::GenericInputMPVData(true);
    read->setReadDir("input");
    read->setFilename("generic_inputMPVData_00002.csv");
    read->setBinary(false);

    /* Creation of Generic output block to write a set of
     * coordinates of cloud points as generic output in csv format.
     */
    mimmo::GenericOutputMPVData * write = new mimmo::GenericOutputMPVData();
    write->setWriteDir("./");
    write->setFilename("generic_outputMPVData_00002.csv");
    write->setCSV(true);
    write->setBinary(false);

    /* Setup pin connections.
     */
    mimmo::pin::addPin(read, write, M_VECTORFIELD, M_VECTORFIELD);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
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

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        /**<Calling mimmo Test routines*/
        try{
            test00002();
        }
        catch(std::exception & e){
            std::cout<<"genericinput_example_00002 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;
}
