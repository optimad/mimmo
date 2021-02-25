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
#include "mimmo_parallel.hpp"
#include "mimmo_manipulators.hpp"
#include "mimmo_iogeneric.hpp"
#include <bitpit_common.hpp>

// =================================================================================== //
/*!
	example parallel_example_00002.cpp

	brief Example of Free form manipulation applied to a target partitioned geometry.

	Using: MimmoGeometry, Partition, FFDLattice, GenericInput, Apply, Chain.

	<b>To run</b>: mpirun -np X parallel_example_00002 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00002() {

    /*
        Read a STL sphere
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");

    /*
        Write the partitioned geometry in parallel vtu format
     */
    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::SURFVTU);
    mimmo1->setWriteFilename("parallel_output_00002.0001");

    /*
        Write the partitioned deformed geometry in parallel vtu format
     */
    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo2->setWriteDir(".");
    mimmo2->setWriteFileType(FileType::SURFVTU);
    mimmo2->setWriteFilename("parallel_output_00002.0002");

    /*
        Distribute the target mesh among processes.
     */
    mimmo::Partition* partition = new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition->setPlotInExecution(false);

    /*
        FFDLattice manipulator - Sphere shaped - creation .
        It will require the sphere dimensions and lattice number of
        nodes and nurbs degrees
     */
    mimmo::FFDLattice* lattice = new mimmo::FFDLattice();
    darray3E origin = {0.0, 0.0,0.0};
    darray3E span;
    span[0]= 3.01;
    span[1]= 2*BITPIT_PI;
    span[2]= BITPIT_PI;

    /*
        Set number of nodes of the mesh (dim) and degree of nurbs functions (deg).
     */
    iarray3E dim, deg;
    dim[0] = 30;
    dim[1] = 30;
    dim[2] = 30;
    deg[0] = 2;
    deg[1] = 2;
    deg[2] = 2;

    lattice->setLattice(origin,span,mimmo::ShapeType::SPHERE,dim, deg);

    /*
        Change reference system to work in local spherical coordinates.
        Set coordinates as CLAMPED (continuity in origins of angles).
     */
    lattice->setRefSystem(2, darray3E{0,1,0});
    lattice->setCoordType(mimmo::CoordType::CLAMPED, 2);
    lattice->setPlotInExecution(false);

    /*
        Build mesh of lattice outside the execution chain
        to use it during setup the displacements.
     */
    lattice->build();

    /*
        Use random values to set the displacements of the control nodes at a longitude
     * angle smaller than PI and expansion on radius direction for nodes with longitude
     * angle greater than PI.
     */
    int ndeg = lattice->getNNodes();
    dvecarr3E displ(ndeg, darray3E{0,0,0});
    time_t Time = time(NULL);
    srand(Time);
    for (int i=0; i<ndeg; i++){
        int l1,l2,l3;
        int index = lattice->accessGridFromDOF(i);
        lattice->accessPointIndex(index,l1,l2,l3);
        if(l1 > 0 && lattice->getLocalPoint(l1,l2,l3)[1] < BITPIT_PI){
            displ[i][0] = 1.0*( (double) (rand()) / RAND_MAX );
        }
        if( (l1 > 0 && lattice->getLocalPoint(l1,l2,l3)[1] >= BITPIT_PI)
                || lattice->getLocalPoint(l1,l2,l3)[1] == 0){
            displ[i][0] = 1.25;
        }

    }

    /*
        Set Generic input block with the
        displacements defined above.
     */
    mimmo::GenericInput* input = new mimmo::GenericInput();
    input->setInput(displ);

    /*
        Create applier block.
        It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /*
        Setup pin connections.
     */
    mimmo::pin::addPin(mimmo0, partition, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, mimmo1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, lattice, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(input, lattice, M_DISPLS, M_DISPLS);
    mimmo::pin::addPin(lattice, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applier, mimmo2, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(input);
    ch0.addObject(partition);
    ch0.addObject(applier);
    ch0.addObject(lattice);
    ch0.addObject(mimmo0);
    ch0.addObject(mimmo1);
    ch0.addObject(mimmo2);

    /*
        Execute the chain.
        Use debug flag false to avoid printing intermediate results of the execution steps.
     */
    ch0.exec(true);

    /*
        Clean up & exit;
     */
    delete lattice;
    delete applier;
    delete input;
    delete mimmo0;
    delete partition;
    delete mimmo1;
    delete mimmo2;

    return;
}

int	main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);

    {
#endif
        try{
            /**<Calling core function*/
            test00002() ;
        }
        catch(std::exception & e){
            std::cout<<"parallel_example_00002 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    }
    MPI_Finalize();
#endif

    return 0;
}
