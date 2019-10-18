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
#include "MRBF.hpp"
#include <exception>
using namespace mimmo;


/*!
 * \example genericinput_example_00004.cpp
 *
 * \brief Example of reading point positions and displacements related to a RBF manipulator set

 * Using: IOCloudPoints, MRBF
 *
 * <b>To run</b>: ./genericinput_example_00003 \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //


void test00004() {

    /* Reading plane
     */
     /* Reading plane
      */
     MimmoGeometry * read = new MimmoGeometry();
     read->setIOMode(IOMode::CONVERT);
     read->setReadDir("geodata");
     read->setReadFilename("plane4");
     read->setReadFileType(FileType::STL);
     read->setWriteDir(".");
     read->setWriteFilename("./genericinput_output_00004.0001");
     read->setWriteFileType(FileType::SURFVTU);


    /* Reading rbf set of points and their diplacements from file
     */
    IOCloudPoints * iocp = new IOCloudPoints(true);
    iocp->setName("genericinput_example_00004_RBFSet");
    iocp->setReadDir("input");
    iocp->setReadFilename("generic_iocloud.txt");
    iocp->setPlotInExecution(true);
    /*!Create RBF manipulator, specyfing only the RBF shape and support Radius
    */
    MRBF * rbf = new MRBF();
    rbf->setFunction(bitpit::RBFBasisFunction::C1C1);
    rbf->setSupportRadiusValue(0.6);
    rbf->setApply(true);
    /* Setup pin connections.
     */
    pin::addPin(read, rbf, M_GEOM, M_GEOM);
    pin::addPin(iocp, rbf, M_COORDS, M_COORDS);
    pin::addPin(iocp, rbf, M_DISPLS, M_DISPLS);

    /* Setup execution chain.
     */
    Chain ch0;
    ch0.addObject(read);
    ch0.addObject(iocp);
    ch0.addObject(rbf);


    /* Execution of chain.
     * Use debug flag true to print out the execution steps.
     */
    ch0.exec(true);

    /*
        print the deformed geometry
    */
    read->getGeometry()->getPatch()->write("./genericinput_output_00004.0002");

    /* Clean up & exit;
     */
    delete read;
    delete iocp;
    delete rbf;

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
            test00004();
        }
        catch(std::exception & e){
            std::cout<<"genericinput_example_00004 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;
}
