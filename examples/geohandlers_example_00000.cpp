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

// =================================================================================== //
/*!
	\example geohandlers_example_00000.cpp

	\brief Example of usage of selection block to select a sub-patch of an input geometry.

	Geometry handler block used: SelectionByCylinder.

	<b>To run</b>: ./geohandlers_example_00000 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //

void test00001() {

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry();
    mimmo0->setName("mimmo0");
    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("./geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("plane3");
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00000.0000");

    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry();
    mimmo1->setIOMode(IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("geohandlers_output_00000.0001");

    mimmo::SelectionByCylinder * sel3 = new mimmo::SelectionByCylinder();
    sel3->setOrigin({{0.5, 0.0, 0.0}});
    sel3->setSpan({{0.75, 2.0*M_PI, 0.4}});
    sel3->setInfLimits({{0.25, 0.0, 0.0}});
    sel3->setRefSystem({{0,1,0}}, {{0,0,1}}, {{1,0,0}});
    sel3->setPlotInExecution(true);

    /* Setup pin connections.
     */
    mimmo::pin::addPin(mimmo0, sel3, M_GEOM, M_GEOM);
    mimmo::pin::addPin(sel3, mimmo1, M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(sel3);
    ch0.addObject(mimmo1);

    /* Execution of chain.
     * Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
	delete sel3;
    delete mimmo0;
    delete mimmo1;

    sel3 = NULL;
    mimmo0 = NULL;
    mimmo1 = NULL;

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
