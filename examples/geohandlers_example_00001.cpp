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
#include <bitpit_common.hpp>

// =================================================================================== //
/*!
	\example geohandlers_example_00001.cpp

	\brief Example of usage of selection and refinement of a target geometry sub-patch.

	Using: MimmoGeometry, SelectionByBox, SelectionBySphere, RefineGeometry, Chain,
           Partition (only MPI ).

	<b>To run</b>              : ./geohandlers_example_00001 \n
    <b>To run (MPI Version)</b>: mpirun -np X geohandlers_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //

void test00001() {

    /*
        Read a target geometry from file. CONVERT option will let the block to write
        the just read file in another file, immediately after the reading.
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setName("mimmo0");
    mimmo0->setReadDir("./geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00001.0000");

    /*
        Write a geometry to file
    */
    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::SURFVTU);
    mimmo1->setWriteFilename("geohandlers_output_00001.0001");

    /*
        Write a geometry to file
    */
    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo2->setWriteDir(".");
    mimmo2->setWriteFileType(FileType::SURFVTU);
    mimmo2->setWriteFilename("geohandlers_output_00001.0002");

#if MIMMO_ENABLE_MPI
    /*
        Distribute a target mesh geometry among processors
    */
    mimmo::Partition* partition = new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition->setPlotInExecution(true);
#endif

    /*
        Select a portion of geometry incapsulated inside a Box/Cube.
        Span and origin of the box are needed.
     */
    mimmo::SelectionByBox * boxSel = new mimmo::SelectionByBox();
	boxSel->setOrigin({{-0.5,-0.5,0.2}});
	boxSel->setSpan(0.6,0.6,0.6);

    /*
        Select a portion of geometry incapsulated inside a Sphere.
        Radius, angular spans and origin of the sphere are needed.
     */
	mimmo::SelectionBySphere * sphSel = new mimmo::SelectionBySphere();
	sphSel->setOrigin({{-0.5, 0.5,0.2}});
	sphSel->setSpan(0.34, 2*BITPIT_PI, BITPIT_PI);

    /*
       Block to refine triangles of a target mesh
       Refine Engine type and number of refinement steps are needed.
       Smoothing steps >0 activate a n-steps procedure to attempt smoothing/regularizing
       the final refined mesh.
     */
	mimmo::RefineGeometry * refine = new mimmo::RefineGeometry();
    refine->setRefineType(mimmo::RefineType::REDGREEN);
    refine->setRefineSteps(2);
    refine->setSmoothingSteps(1);

    /*
        Define block pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, partition, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, boxSel, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, sphSel, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(mimmo0, boxSel, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, sphSel, M_GEOM, M_GEOM);
#endif
    mimmo::pin::addPin(boxSel, refine, M_GEOM, M_GEOM);
    mimmo::pin::addPin(refine, mimmo1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(sphSel, mimmo2, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition);
#endif
    ch0.addObject(boxSel);
    ch0.addObject(refine);
    ch0.addObject(sphSel);
    ch0.addObject(mimmo1);
    ch0.addObject(mimmo2);

    /*
        Execute the chain.
        Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

    /*
        Clean up & exit;
     */
	delete boxSel;
	delete sphSel;
    delete mimmo0;
#if MIMMO_ENABLE_MPI
    delete partition;
#endif
    delete refine;
    delete mimmo1;
    delete mimmo2;

	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
		try{
            /**<Calling core function*/
            test00001();
        }
        catch(std::exception & e){
            std::cout<<"geohandlers_example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
