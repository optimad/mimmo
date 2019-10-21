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
#include "mimmo_manipulators.hpp"

// =================================================================================== //
/*!
	\example geohandlers_example_00003.cpp

	\brief Example of usage of selection block to select sub-patches from an input geometry,
           define fields onto them, and reconstruct the final field on the original geometry.

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
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry();
    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00003.0000");

    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry();
    mimmo1->setIOMode(IOMode::CONVERT);
    mimmo1->setReadDir("geodata");
    mimmo1->setReadFileType(FileType::STL);
    mimmo1->setReadFilename("plane1");
    mimmo1->setWriteDir("./");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("geohandlers_output_00003p1.0000");

    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry();
    mimmo2->setIOMode(IOMode::CONVERT);
    mimmo2->setReadDir("geodata");
    mimmo2->setReadFileType(FileType::STL);
    mimmo2->setReadFilename("plane2");
    mimmo2->setWriteDir("./");
    mimmo2->setWriteFileType(FileType::STL);
    mimmo2->setWriteFilename("geohandlers_output_00003p2.0000");

    mimmo::MimmoGeometry * mimmo3 = new mimmo::MimmoGeometry();
    mimmo3->setIOMode(IOMode::WRITE);
    mimmo3->setWriteDir(".");
    mimmo3->setWriteFileType(FileType::STL);
    mimmo3->setWriteFilename("geohandlers_output_00003.0001");

    mimmo::MimmoGeometry * mimmo4 = new mimmo::MimmoGeometry();
    mimmo4->setIOMode(IOMode::WRITE);
    mimmo4->setWriteDir(".");
    mimmo4->setWriteFileType(FileType::STL);
    mimmo4->setWriteFilename("geohandlers_output_00003.0002");

    mimmo::MimmoGeometry * mimmo5 = new mimmo::MimmoGeometry();
    mimmo5->setIOMode(IOMode::WRITE);
    mimmo5->setWriteDir(".");
    mimmo5->setWriteFileType(FileType::STL);
    mimmo5->setWriteFilename("geohandlers_output_00003.0003");


    /* Instantiation of two Selection By Map block.
     * The planes are used as selection objects with an offset
     * defined by the user.
     */
    mimmo::SelectionByMapping  * mapSel1 = new mimmo::SelectionByMapping();
    mimmo::SelectionByMapping  * mapSel2 = new mimmo::SelectionByMapping();
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
    mimmo::MRBF* mrbf1 = new mimmo::MRBF();
    mrbf1->setMode(mimmo::MRBFSol::NONE);
    mrbf1->setSupportRadius(0.4);
    mrbf1->setPlotInExecution(true);
    mrbf1->setNode(rbfNodes1);
    mrbf1->setDisplacements(displ1);

    mimmo::MRBF* mrbf2 = new mimmo::MRBF();
    mrbf2->setMode(mimmo::MRBFSol::NONE);
    mrbf2->setSupportRadius(0.4);
    mrbf2->setPlotInExecution(true);
    mrbf2->setNode(rbfNodes2);
    mrbf2->setDisplacements(displ2);

    /* Create reconstruct vector block and set to recontruct over the whole
     * input geometry the displacements fields
     * given by the two rbf blocks on two separate patches.
     */
    mimmo::ReconstructVector* recon = new mimmo::ReconstructVector();

    /* Create applier block.
     * It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /* Create extract vector field block and set to extract by ID over an
     * input geometry the displacements fields
     * given by the two reconstruction block on two unified patches.
     */
    mimmo::ExtractVectorField* extr = new mimmo::ExtractVectorField();
    extr->setMode(mimmo::ExtractMode::ID);

    /* Create applier extraction block.
     * It applies the extracted deformation displacements
     * to the selected input geometry.
     */
    mimmo::Apply* applierextr = new mimmo::Apply();

    /* Create extract vector field block and set to extract by MAPPING over an
     * input geometry the displacements fields
     * given by the two reconstruction block on two unified patches.
     */
    mimmo::ExtractVectorField* extr2 = new mimmo::ExtractVectorField();
    extr2->setMode(mimmo::ExtractMode::MAPPING);
    extr2->setTolerance(1.0e-01);

    /* Create applier extraction2 block.
     * It applies the extracted deformation displacements
     * to the selected input geometry.
     */
    mimmo::Apply* applierextr2 = new mimmo::Apply();

    /* Setup pin connections.
     */
    mimmo::pin::addPin(mimmo0, mapSel1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, mapSel2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, mapSel1, M_GEOM, M_GEOM2);
    mimmo::pin::addPin(mimmo2, mapSel2, M_GEOM, M_GEOM2);
    mimmo::pin::addPin(mapSel1, mrbf1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mapSel2, mrbf2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, recon, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mrbf1, recon, M_GDISPLS, M_VECTORFIELD);
    mimmo::pin::addPin(mrbf2, recon, M_GDISPLS, M_VECTORFIELD);
    mimmo::pin::addPin(recon, applier, M_VECTORFIELD, M_GDISPLS);
    mimmo::pin::addPin(applier, mimmo3, M_GEOM, M_GEOM);

    mimmo::pin::addPin(mapSel1, extr, M_GEOM, M_GEOM);
    mimmo::pin::addPin(recon, extr, M_VECTORFIELD, M_VECTORFIELD);
    mimmo::pin::addPin(mapSel1, applierextr, M_GEOM, M_GEOM);
    mimmo::pin::addPin(extr, applierextr, M_VECTORFIELD, M_GDISPLS);
    mimmo::pin::addPin(applierextr, mimmo4, M_GEOM, M_GEOM);

    mimmo::pin::addPin(mapSel2, extr2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(recon, extr2, M_VECTORFIELD, M_VECTORFIELD);
    mimmo::pin::addPin(mapSel2, applierextr2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(extr2, applierextr2, M_VECTORFIELD, M_GDISPLS);
    mimmo::pin::addPin(applierextr2, mimmo5, M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
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
    ch0.addObject(extr);
    ch0.addObject(applierextr);
    ch0.addObject(mimmo4);
    ch0.addObject(extr2);
    ch0.addObject(applierextr2);
    ch0.addObject(mimmo5);

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
    delete mimmo4;
    delete mimmo5;
    delete mapSel1;
    delete mapSel2;
    delete mrbf1;
    delete mrbf2;
    delete recon;
    delete extr;
    delete extr2;
    delete applier;
    delete applierextr;
    delete applierextr2;

    mimmo0 = NULL;
    mimmo1 = NULL;
    mimmo2 = NULL;
    mimmo3 = NULL;
    mimmo4 = NULL;
    mimmo5 = NULL;
    mapSel1 = NULL;
    mapSel2 = NULL;
    mrbf1 = NULL;
    mrbf2 = NULL;
    recon = NULL;
    extr = NULL;
    extr2 = NULL;
    applier = NULL;
    applierextr = NULL;
    applierextr2 = NULL;

	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
		/**<Calling mimmo Test routine*/
		try{
            test00003();
        }
        catch(std::exception & e){
            std::cout<<"geohandlers_example_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
