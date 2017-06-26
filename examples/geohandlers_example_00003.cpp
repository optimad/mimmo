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
#include "mimmo_manipulators.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmo::pin;

// =================================================================================== //
/*!
	\example geohandlers_example_00003.cpp

	\brief Example of usage of selection block to select a sub-patch of an input geometry.

    Geometry handler block used: SelectionByMap, ReconstructVector

	<b>To run</b>: ./geohandlers_example_00003 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //

void test00003() {

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
    MimmoGeometry * mimmo0 = new MimmoGeometry();
    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00003.0000");

    MimmoGeometry * mimmo1 = new MimmoGeometry();
    mimmo1->setIOMode(IOMode::CONVERT);
    mimmo1->setReadDir("geodata");
    mimmo1->setReadFileType(FileType::STL);
    mimmo1->setReadFilename("plane1");
    mimmo1->setWriteDir("./");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("geohandlers_output_00003p1.0000");

    MimmoGeometry * mimmo2 = new MimmoGeometry();
    mimmo2->setIOMode(IOMode::CONVERT);
    mimmo2->setReadDir("geodata");
    mimmo2->setReadFileType(FileType::STL);
    mimmo2->setReadFilename("plane2");
    mimmo2->setWriteDir("./");
    mimmo2->setWriteFileType(FileType::STL);
    mimmo2->setWriteFilename("geohandlers_output_00003p2.0000");

    MimmoGeometry * mimmo3 = new MimmoGeometry();
    mimmo3->setIOMode(IOMode::WRITE);
    mimmo3->setWriteDir(".");
    mimmo3->setWriteFileType(FileType::STL);
    mimmo3->setWriteFilename("geohandlers_output_00003.0001");


    /* Instantiation of two Selection By Map block.
     * The planes are used as selection objects with an offset
     * defined by the user.
     */
    SelectionByMapping  * mapSel1 = new SelectionByMapping();
    SelectionByMapping  * mapSel2 = new SelectionByMapping();
    mapSel1->setTolerance(1.0e-01);
    mapSel1->setPlotInExecution(true);
    mapSel2->setTolerance(1.0e-01);
    mapSel2->setPlotInExecution(true);

    /* Creation of a two points with x coordinates -0.5 and 0.5 (y=z=0.0);
     */
    dvecarr3E rbfNodes1(1, {{-0.5, 0.0, 0.0}});
    dvecarr3E rbfNodes2(1, {{0.5, 0.0, 0.0}});

    /* Creation of a two displacements for the rbf control points.
     */
    dvecarr3E displ1(1, {{-0.25, 0.0, 0.0}});
    dvecarr3E displ2(1, {{0.25, 0.0, 0.0}});

    /* Instantiation of two MRBF objects.
     * Set rbf points and displacements defined above.
     * Plot Optional results during execution active for MRBF block.
     */
    MRBF* mrbf1 = new MRBF();
    mrbf1->setMode(MRBFSol::NONE);
    mrbf1->setSupportRadius(0.4);
    mrbf1->setPlotInExecution(true);
    mrbf1->setNode(rbfNodes1);
    mrbf1->setDisplacements(displ1);

    MRBF* mrbf2 = new MRBF();
    mrbf2->setMode(MRBFSol::NONE);
    mrbf2->setSupportRadius(0.4);
    mrbf2->setPlotInExecution(true);
    mrbf2->setNode(rbfNodes2);
    mrbf2->setDisplacements(displ2);

    /* Create reconstruct vector block and set to recontruct over the whole
     * input geometry the displacements fields
     * given by the two rbf blocks on two separate patches.
     */
    ReconstructVector* recon = new ReconstructVector();


    /* Create applier block.
     * It applies the deformation displacements to the original input geometry.
     */
    Apply* applier = new Apply();

    /* Setup pin connections.
     */
    addPin(mimmo0, mapSel1, PortType::M_GEOM, PortType::M_GEOM);
    addPin(mimmo0, mapSel2, PortType::M_GEOM, PortType::M_GEOM);
    addPin(mimmo0, applier, PortType::M_GEOM, PortType::M_GEOM);
    addPin(mimmo1, mapSel1, PortType::M_GEOM, PortType::M_GEOM2);
    addPin(mimmo2, mapSel2, PortType::M_GEOM, PortType::M_GEOM2);
    addPin(mapSel1, mrbf1, PortType::M_GEOM, PortType::M_GEOM);
    addPin(mapSel2, mrbf2, PortType::M_GEOM, PortType::M_GEOM);
    addPin(mimmo0, recon, PortType::M_GEOM, PortType::M_GEOM);
    addPin(mrbf1, recon, PortType::M_GDISPLS, PortType::M_GDISPLS);
    addPin(mrbf2, recon, PortType::M_GDISPLS, PortType::M_GDISPLS);
    addPin(recon, applier, PortType::M_GDISPLS, PortType::M_GDISPLS);
    addPin(applier, mimmo3, PortType::M_GEOM, PortType::M_GEOM);


    /* Setup execution chain.
     */
    Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(mimmo1);
    ch0.addObject(mimmo2);
    ch0.addObject(mapSel1);
    ch0.addObject(mapSel2);
    ch0.addObject(mrbf1);
    ch0.addObject(mrbf2);
    ch0.addObject(recon);
    ch0.addObject(applier);
    ch0.addObject(mimmo3);

    /* Execution of chain.
     * Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
    delete mimmo0;
    delete mimmo1;
    delete mimmo2;
    delete mimmo3;
    delete mapSel1;
    delete mapSel2;
    delete mrbf1;
    delete mrbf2;
    delete recon;
    delete applier;

    mimmo0 = NULL;
    mimmo1 = NULL;
    mimmo2 = NULL;
    mimmo3 = NULL;
    mapSel1 = NULL;
    mapSel2 = NULL;
    mrbf1 = NULL;
    mrbf2 = NULL;
    recon = NULL;
    applier = NULL;

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
		test00003();
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

