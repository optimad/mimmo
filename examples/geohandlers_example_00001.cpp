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

#include "mimmo_geohandlers.hpp"
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif
#include <bitpit_common.hpp>

// =================================================================================== //
/*!
	\example geohandlers_example_00001.cpp

	\brief Example of usage of selection block to select a sub-patch of an input geometry.

	Geometry handler block used: SelectionByBox, SelectionBySphere.

	<b>To run</b>: ./geohandlers_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //

void test00001() {

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setName("mimmo0");
    mimmo0->setReadDir("./geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00001.0000");

    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::SURFVTU);
    mimmo1->setWriteFilename("geohandlers_output_00001.0001");

    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo2->setWriteDir(".");
    mimmo2->setWriteFileType(FileType::SURFVTU);
    mimmo2->setWriteFilename("geohandlers_output_00001.0002");

#if MIMMO_ENABLE_MPI
    /* Instantiation of a Partition object with default partition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition* partition = new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition->setPlotInExecution(true);
#endif

    /* Instantiation of a Selection By Box block.
     * Setup of span and origin of cube.
     */
    mimmo::SelectionByBox * boxSel = new mimmo::SelectionByBox();
	boxSel->setOrigin({{-0.5,-0.5,0.2}});
	boxSel->setSpan(0.6,0.6,0.6);

    /* Instantiation of a Selection By Sphere block.
     * Setup of span and origin of sphere.
     */
	mimmo::SelectionBySphere * sphSel = new mimmo::SelectionBySphere();
	sphSel->setOrigin({{-0.5, 0.5,0.2}});
	sphSel->setSpan(0.34, 2*BITPIT_PI, BITPIT_PI);

    /* Instantiation of a Refine Geometry block.
     * Setup refining.
     */
	mimmo::RefineGeometry * refine = new mimmo::RefineGeometry();
    refine->setRefineType(mimmo::RefineType::REDGREEN);
    refine->setRefineSteps(2);
    refine->setSmoothingSteps(1);

    /* Setup pin connections.
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

    /* Setup execution chain.
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

    /* Execution of chain.
     * Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
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

    boxSel = NULL;
    sphSel = NULL;
    mimmo0 = NULL;
#if MIMMO_ENABLE_MPI
    partition = NULL;
#endif
    refine = NULL;
    mimmo1 = NULL;
    mimmo2 = NULL;

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
            /**<Calling mimmo Test routine*/
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
