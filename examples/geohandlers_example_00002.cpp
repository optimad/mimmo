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

#include "mimmo_geohandlers.hpp"
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

// =================================================================================== //
/*!
	\example geohandlers_example_00002.cpp

	\brief Example of usage of stitching block of two input geometries.

	Using: MimmoGeometry, StitchGeometry, Chain, Partition(MPI version)

	<b>To run</b>              : ./geohandlers_example_00002 \n
    <b>To run (MPI version)</b>: mpirun -np X geohandlers_example_00002 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //

void test00002() {

    /*
        Read a sphere from STL file. Convert mode is to save the just read geometry in
        another file with name geohandlers_output_00002.0000.stl
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
	mimmo0->setReadDir("geodata");
	mimmo0->setReadFilename("sphere2");
	mimmo0->setReadFileType(FileType::STL);
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00002.0000");

    /*
        Read the Stanford bunny from STL file. Convert mode is to save the just read geometry in
        another file with name geohandlers_output_00002.0001.stl
     */
    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo1->setReadDir("geodata");
    mimmo1->setReadFilename("stanfordBunny2");
    mimmo1->setReadFileType(FileType::STL);
    mimmo1->setWriteDir("./");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("geohandlers_output_00002.0001");

    /*
        Write the stiched geometry to STL file.
     */
    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo2->setWriteDir("./");
    mimmo2->setWriteFileType(FileType::STL);
    mimmo2->setWriteFilename("geohandlers_output_00002.0002");

#if MIMMO_ENABLE_MPI

    /* Block to distribute among processors the sphere geometry */
    mimmo::Partition* partition0 = new mimmo::Partition();
    partition0->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition0->setPlotInExecution(true);

    /* Block to distribute among processors the bunny geometry */
    mimmo::Partition* partition1 = new mimmo::Partition();
    partition1->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition1->setPlotInExecution(true);
#endif

    /*
       Stitcher of multiple geometries.
     * Plot Optional results during execution active for Stitcher block.
     */
    mimmo::StitchGeometry * stitcher = new mimmo::StitchGeometry(1);
	stitcher->setPlotInExecution(true);
	stitcher->setOutputPlot(".");

    /*
        Setup block pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, partition0, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition0, stitcher, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, partition1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition1, stitcher, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(mimmo0, stitcher, M_GEOM, M_GEOM);
	mimmo::pin::addPin(mimmo1, stitcher, M_GEOM, M_GEOM);
#endif
	mimmo::pin::addPin(stitcher, mimmo2, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
	mimmo::Chain ch0;
	ch0.addObject(mimmo0);
	ch0.addObject(mimmo1);
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition0);
    ch0.addObject(partition1);
#endif
	ch0.addObject(stitcher);
	ch0.addObject(mimmo2);

    /*
        Execution the chain.
        Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
    delete mimmo0;
    delete mimmo1;
#if MIMMO_ENABLE_MPI
    delete partition0;
    delete partition1;
#endif
    delete mimmo2;
    delete stitcher;

	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
		/**<Calling core function*/
        try{
            test00002() ;
        }
        catch(std::exception & e){
            std::cout<<"geohandlers_example_00002 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
