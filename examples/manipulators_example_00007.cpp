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

#include "mimmo_manipulators.hpp"
#include "mimmo_iogeneric.hpp"
#include "mimmo_utils.hpp"
#include <random>

// =================================================================================== //
/*!
	\example manipulators_example_00007.cpp

	\brief Example of usage of radial basis function block/ with variable support radius
     to manipulate an input geometry.

    Geometry deformation block used: MRBF.
    Utils block used: ProjectCloud.

	<b>To run</b>: ./manipulators_example_00007 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00007() {

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
    mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("manipulators_output_00007.0000");

    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("manipulators_output_00007.0001");

    /* Creation of a random distribution of 10 points with coordinates between [-0.5, 0.5]
       and a random distribution of 10 support radii between [0.05, 0.35]
     */
    int np = 10;
    dvecarr3E rbfNodes(10);
    dvector1D supportRadii(10);
    std::minstd_rand rgen, sgen;
    rgen.seed(160);
    sgen.seed(333);
    double dist = (rgen.max()-rgen.min());
    double dist2 = (sgen.max()-sgen.min());

    for (int i=0; i<np; i++){
        for (int j=0; j<3; j++){
            rbfNodes[i][j] = 1.0*( double(rgen() - rgen.min() ) / dist ) - 0.5;
        }
        supportRadii[i] = 0.3 * ( double(sgen() - sgen.min() ) / dist2 ) + 0.05;
    }

    /* Set Generic input block with the
     * nodes defined above.
     */
    mimmo::GenericInput* inputn = new mimmo::GenericInput();
    inputn->setInput(rbfNodes);

    /* Creation of a projection block aimed to project
     * the point cloud previously defined over the input geometry.
     * The rbf control nodes will be defined on the surface of the
     * geometry as output of this block.
     *
     */
    mimmo::ProjectCloud* proj = new mimmo::ProjectCloud();

    /* Instantiation of a MRBF object with a distribution of 10 random control nodes projected
     * ont he input surface.
     * Plot Optional results during execution active for MRBF block.
     */
    mimmo::MRBF* mrbf = new mimmo::MRBF(mimmo::MRBFSol::NONE);
    mrbf->setFunction(bitpit::RBFBasisFunction::WENDLANDC2);
    mrbf->setPlotInExecution(true);


    /* Creation of a set of displacements of the control nodes of the radial basis functions.
     * Use a radial displacements from a center point placed in axes origin.
     */
    dvecarr3E displ(np, darray3E{0.0, 0.0, 0.0});
    darray3E center({0.0, 0.0, 0.0});
    for (int i=0; i<np; i++){
        displ[i] = rbfNodes[i] - center;
        displ[i] /= 2.0*norm2(displ[i]);
    }


    /* Set Generic input block with the
     * displacements defined above.
     */
    mimmo::GenericInput* input = new mimmo::GenericInput();
    input->setInput(displ);

    /* Set Generic input block with the
     * supportRadii defined above.
     */
    mimmo::GenericInput* inputSR = new mimmo::GenericInput();
    inputSR->setInput(supportRadii);


    /* Create applier block.
     * It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /* Setup pin connections.
     */
    mimmo::pin::addPin(mimmo0, mrbf, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, proj, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(inputn, proj, M_COORDS, M_COORDS);
    mimmo::pin::addPin(proj, mrbf, M_COORDS, M_COORDS);
    mimmo::pin::addPin(input, mrbf, M_DISPLS, M_DISPLS);
    mimmo::pin::addPin(inputSR, mrbf, M_DATAFIELD, M_DATAFIELD);
    mimmo::pin::addPin(mrbf, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applier, mimmo1, M_GEOM, M_GEOM);

    /* Setup execution chain.
     * The object can be insert in the chain in random order.
     * The chain object recover the correct order of execution from
     * the pin connections.
     */
    mimmo::Chain ch0;
    ch0.addObject(input);
    ch0.addObject(inputn);
    ch0.addObject(inputSR);
    ch0.addObject(mimmo0);
    ch0.addObject(proj);
    ch0.addObject(applier);
    ch0.addObject(mrbf);
    ch0.addObject(mimmo1);

    /* Execution of chain.
     * Use debug flag true to print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
    delete mrbf;
    delete proj;
    delete applier;
    delete input;
    delete inputSR;
    delete mimmo0;
    delete mimmo1;

    proj    = NULL;
    mrbf    = NULL;
    applier = NULL;
    input   = NULL;
    inputSR = NULL;
    mimmo0  = NULL;
    mimmo1  = NULL;

    return;
}

int	main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        try{
            /**<Calling mimmo Test routine*/
            test00007() ;
        }
        catch(std::exception & e){
            std::cout<<"manipulators_example_00007 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;
}
