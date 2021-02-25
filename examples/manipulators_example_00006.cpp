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

#include "mimmo_manipulators.hpp"
#include "mimmo_iogeneric.hpp"
#include "mimmo_utils.hpp"
#include <random>
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

// =================================================================================== //
/*!
	\example manipulators_example_00006.cpp

	\brief Example of usage of RBF block to manipulate an input geometry.

    Using: MimmoGeometry, GenericInput, CreatePointCloud, ProjPatchOnSurface,
           MRBF, Apply, Chain, Partition(MPI version).

	<b>To run</b>              : ./manipulators_example_00006 \n
    <b>To run (MPI version)</b>: mpirun -np X manipulators_example_00006 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00006() {

    /*
        Read a sphere from STL. Convert mode is to save the just read geometry in
        another file with name manipulators_output_00006.0000.stl
     */
    mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("manipulators_output_00006.0000");

    /* write the deformed mesh to file */
    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("manipulators_output_00006.0001");

#if MIMMO_ENABLE_MPI
    /*
        Distribute mesh among processes
     */
    mimmo::Partition* partition= new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
#endif

    /*
        Creation of a random distribution of 10 points with coordinates between [-0.5, 0.5]
     */
    int np = 10;
    dvecarr3E rbfNodes(10);
    std::minstd_rand rgen;
    rgen.seed(160);
    double dist = (rgen.max()-rgen.min());
    for (int i=0; i<np; i++){
        for (int j=0; j<3; j++)
            rbfNodes[i][j] = 1.0*( double(rgen() - rgen.min() ) / dist ) - 0.5;
    }

    /*
        Set a Generic input block with the
        nodes defined above.
     */
    mimmo::GenericInput* inputn = new mimmo::GenericInput();
    inputn->setInput(rbfNodes);

    /*
        Creation of a set of displacements of the control nodes of the radial basis functions.
        Use a radial displacements from a center point placed in axes origin.
     */
    dvecarr3E displ(np, darray3E{0.0, 0.0, 0.0});
    darray3E center({0.0, 0.0, 0.0});
    for (int i=0; i<np; i++){
        displ[i] = rbfNodes[i] - center;
        displ[i] /= 2.0*norm2(displ[i]);
    }

    /*
        Set Generic input block with the
        displacements defined above.
     */
    mimmo::GenericInput* input = new mimmo::GenericInput();
    input->setInput(displ);

    /*
       Create a Point Cloud that will use the outputs of previous generic inputs.
      Nodes and displacements will be stored as a point cloud mesh with a vector
      field attached
     */
    mimmo::CreatePointCloud * rbfPC = new mimmo::CreatePointCloud();

    /*
        This block will project the rbf point cloud onto the target surface
     */
    mimmo::ProjPatchOnSurface* proj = new mimmo::ProjPatchOnSurface();
    proj->setWorkingOnTarget(true);
    proj->setPlotInExecution(true);

    /*
        RBF manipulator, will use RBF point cloud and its vectorfield of displacements
        to build up itself.
        It requires the definition of RBF function and support radius
     */
    mimmo::MRBF* mrbf = new mimmo::MRBF(mimmo::MRBFSol::NONE);
    mrbf->setFunction(bitpit::RBFBasisFunction::WENDLANDC2);
    mrbf->setSupportRadiusReal(0.1);
    mrbf->setPlotInExecution(true);

    /*
        Create applier block.
        It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /* Setup pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, partition, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, mrbf, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(mimmo0, mrbf, M_GEOM, M_GEOM);
#endif
    mimmo::pin::addPin(mimmo0, proj, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(inputn, rbfPC, M_COORDS, M_COORDS);
    mimmo::pin::addPin(input, rbfPC, M_DISPLS, M_DISPLS);
    mimmo::pin::addPin(rbfPC, proj, M_GEOM, M_GEOM2);
    mimmo::pin::addPin(rbfPC, mrbf, M_VECTORFIELD, M_VECTORFIELD);
    mimmo::pin::addPin(proj, mrbf, M_GEOM, M_GEOM2);
    mimmo::pin::addPin(mrbf, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applier, mimmo1, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
    mimmo::Chain ch0;
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition);
#endif
    ch0.addObject(input);
    ch0.addObject(inputn);
    ch0.addObject(rbfPC);
    ch0.addObject(mimmo0);
    ch0.addObject(proj);
    ch0.addObject(applier);
    ch0.addObject(mrbf);
    ch0.addObject(mimmo1);

    /*
        Execute the chain.
     * Use debug flag true to print out the execution steps.
     */
    ch0.exec(true);

    /*
        Clean up & exit;
     */
#if MIMMO_ENABLE_MPI
    delete partition;
#endif
    delete mrbf;
    delete proj;
    delete applier;
    delete input;
    delete rbfPC;
    delete mimmo0;
    delete mimmo1;

    return;
}

int	main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        try{
            /**<Calling core function*/
            test00006() ;
        }
        catch(std::exception & e){
            std::cout<<"manipulators_example_00006 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;
}
