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


#include "bitpit.hpp"
#include "mimmo_geohandlers.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmo::pin;

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
    MimmoGeometry * mimmo0 = new MimmoGeometry();
    mimmo0->setName("mimmo0");
    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("./geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00001.0000");

    MimmoGeometry * mimmo1 = new MimmoGeometry();
    mimmo1->setIOMode(IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("geohandlers_output_00001.0001");

    MimmoGeometry * mimmo2 = new MimmoGeometry();
    mimmo2->setIOMode(IOMode::WRITE);
    mimmo2->setWriteDir(".");
    mimmo2->setWriteFileType(FileType::STL);
    mimmo2->setWriteFilename("geohandlers_output_00001.0002");

    /* Instantiation of a Selection By Box block.
     * Setup of span and origin of cube.
     */
    SelectionByBox      * boxSel = new SelectionByBox();
	boxSel->setOrigin({{-0.5,-0.5,0.2}});
	boxSel->setSpan(0.6,0.6,0.6);
	
    /* Instantiation of a Selection By Sphere block.
     * Setup of span and origin of sphere.
     */
    SelectionBySphere   * sphSel = new SelectionBySphere();
	sphSel->setOrigin({{-0.5, 0.5,0.2}});
	sphSel->setSpan(0.34, 2*M_PI, M_PI);

    /* Setup pin connections.
     */
    pin::addPin(mimmo0, boxSel, M_GEOM, M_GEOM);
    pin::addPin(mimmo0, sphSel, M_GEOM, M_GEOM);
    pin::addPin(boxSel, mimmo1, M_GEOM, M_GEOM);
    pin::addPin(sphSel, mimmo2, M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(boxSel);
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
    delete mimmo1;
    delete mimmo2;

    boxSel = NULL;
    sphSel = NULL;
    mimmo0 = NULL;
    mimmo1 = NULL;
    mimmo2 = NULL;

	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
	
#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling mimmo Test routine*/
		test00001();
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

