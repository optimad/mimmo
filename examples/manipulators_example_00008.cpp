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
	\example manipulators_example_00008.cpp

	\brief Example of usage of radial basis function MRBF block employing scalar input dofs to retrieve
    an output scalar field on the geometry. The field will be interpreted as the entity of displacement
    on the local normal and will be applied to geometry to obtain the final deformation.
    Nodal RBF are distributed on surface with CreateSeedsOnSurface block.

    Using: MimmoGeometry, CreateSeedsOnSurface, GenericInput, MRBF, Apply, Chain, Partition(MPI version).

	<b>To run</b>              : ./manipulators_example_00008 \n
    <b>To run(MPI version)</b> : mpirun -np X manipulators_example_00008 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00008() {

    /*
        Read a pipe sample from vtu mesh. Convert mode is to save the just read geometry in
        another file with name manipulators_output_00008.0000.stl
     */
    mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::SURFVTU);
    mimmo0->setReadFilename("pipeSample");
    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::SURFVTU);
    mimmo0->setWriteFilename("manipulators_output_00008.0000");

    /* write the deformed mesh to file */
    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::SURFVTU);
    mimmo1->setWriteFilename("manipulators_output_00008.0001");


#if MIMMO_ENABLE_MPI
    /*
        Distribute mesh among processes
    */
    mimmo::Partition* partition= new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition->setPlotInExecution(true);
#endif

    /*
        Seed 10 3D points on the surface. Thi will be used as RBF nodes to
        create the manipulator
    */
    mimmo::CreateSeedsOnSurface * seeder = new mimmo::CreateSeedsOnSurface();
    seeder->setNPoints(10);
    seeder->setEngineENUM(mimmo::CSeedSurf::RANDOM);
    seeder->setMassCenterAsSeed(true);
    seeder->setRandomFixed(true);
    seeder->setRandomSignature(12457834);
    seeder->setPlotInExecution(true);

    /*
        Read d.o.f. for RBF Cloud from plain txt file.
       For each node is associated a scalar DOF, in order to obtain a "scalar"
       field of displacements from MRBF
     */
    mimmo::GenericInput* inputDof = new mimmo::GenericInput(true, false);
    inputDof->setReadDir("input");
    inputDof->setFilename("manipulators_input_00008.txt");

    /*
        MRBF manipulator. Setting function and support radius.
        Plot Optional results during execution active for MRBF block.
     */
    mimmo::MRBF* mrbf = new mimmo::MRBF(mimmo::MRBFSol::NONE);
    mrbf->setFunction(bitpit::RBFBasisFunction::C1C2);
    mrbf->setSupportRadiusReal(0.6);
    mrbf->setPlotInExecution(true);

    /*
        It applies the MRBF output scalar field onto the original input geometry, using
       the local vertex normals of the original geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /*
        Setup pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, partition, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, seeder, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, mrbf, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, applier, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(mimmo0, seeder, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, mrbf, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);
#endif

    mimmo::pin::addPin(inputDof, mrbf, M_DATAFIELD, M_DATAFIELD);
    mimmo::pin::addPin(seeder, mrbf, M_COORDS, M_COORDS);
    mimmo::pin::addPin(mrbf, applier, M_SCALARFIELD, M_SCALARFIELD);
    mimmo::pin::addPin(applier, mimmo1, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(inputDof);
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition);
#endif
    ch0.addObject(mimmo0);
    ch0.addObject(seeder);
    ch0.addObject(applier);
    ch0.addObject(mrbf);
    ch0.addObject(mimmo1);

    /*
        Execute the chain.
        Use debug flag true to print out the execution steps.
     */
    ch0.exec(true);

    /*
        Clean up & exit;
     */
    delete mrbf;
    delete seeder;
    delete applier;
    delete inputDof;
    delete mimmo0;
    delete mimmo1;
#if MIMMO_ENABLE_MPI
    delete partition;
#endif

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
            test00008() ;
        }
        catch(std::exception & e){
            std::cout<<"manipulators_example_00008 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;
}
