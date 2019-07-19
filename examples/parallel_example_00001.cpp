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
#include "mimmo_parallel.hpp"
#include "mimmo_iogeneric.hpp"
#include "mimmo_iocgns.hpp"
#include "bitpit.hpp"
#include <exception>
#include <mpi.h>
using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmo::pin;

// =================================================================================== //
/*!
	\example parallel_example_00001.cpp

	\brief Example of usage of parallel partition block to partitioning an input geometry.

	Parallel manipulation block used: Partition.

	<b>To run</b>: ./parallel_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


//TODO - Find  a suitable volume mesh for Partition demo!!!!!!!
void test00001() {

    // /* Creation of mimmo containers.
    //  * Input and output MimmoGeometry are instantiated
    //  * as two different objects (no loop in chain are permitted).
    //  */
    // /* Create IO_CGNS object to import input file. */
    // IOCGNS * mimmo0 = new IOCGNS();
    // mimmo0->setMode(IOCGNS::IOCGNS_Mode::READ);
    // mimmo0->setDir("geodata");
    // mimmo0->setFilename("StaticMixer");
    //
    // MimmoGeometry * mimmo1 = new MimmoGeometry();
    // mimmo1->setIOMode(IOMode::WRITE);
    // mimmo1->setWriteDir(".");
    // mimmo1->setWriteFileType(FileType::VOLVTU);
    // mimmo1->setWriteFilename("parallel_Voutput_00001.0001");
    //
    // MimmoGeometry * mimmo2 = new MimmoGeometry();
    // mimmo2->setIOMode(IOMode::WRITE);
    // mimmo2->setWriteDir(".");
    // mimmo2->setWriteFileType(FileType::SURFVTU);
    // mimmo2->setWriteFilename("parallel_Soutput_00001.0001");
    //
    // /* Instantiation of a Partition object with default patition method space filling curve.
    //  * Plot Optional results during execution active for Partition block.
    //  */
    // Partition* partition = new Partition();
    // partition->setPlotInExecution(true);
    //
    // /* Setup pin connections.
    //  */
    //
    // /* Add pin with port TAG ONLY
    //  */
    //
    // addPin(mimmo0, partition, M_GEOM, M_GEOM);
    // addPin(mimmo0, partition, M_GEOM2, M_GEOM2);
    // addPin(partition, mimmo1, M_GEOM, M_GEOM);
    // addPin(partition, mimmo2, M_GEOM2, M_GEOM);
    //
    // /* Setup execution chain.
    //  */
    // Chain ch0;
    // ch0.addObject(mimmo0);
    // ch0.addObject(partition);
    // ch0.addObject(mimmo1);
    // ch0.addObject(mimmo2);
    //
    // //force the chain to plot all the optional results of its children...
    // ch0.setPlotDebugResults(false);
    // //...in the path specified by the User.
    // ch0.setOutputDebugResults(".");
    //
    // /* Execution of chain.
    //  * Use debug flag true to full print out the execution steps.
    //  */
    // ch0.exec(true);
    //
    // /* Clean up & exit;
    //  */
    // delete partition;
    // delete mimmo0;
    // delete mimmo1;
    // delete mimmo2;
    //
    // partition	= NULL;
    // mimmo0  	= NULL;
    // mimmo1  	= NULL;
    // mimmo2  	= NULL;
    //

}

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);

    {
#endif
        /**<Calling mimmo Test routine*/
        try{
            test00001() ;
        }
        catch(std::exception & e){
            std::cout<<"parallel_example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    }

    MPI_Finalize();
#endif

    return 0;
}
