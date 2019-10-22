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
#include "mimmo_parallel.hpp"
#include "mimmo_manipulators.hpp"
#include "mimmo_iogeneric.hpp"

// =================================================================================== //
/*!
	\example parallel_example_00002.cpp

	\brief Example of usage of free form deformation block to manipulate an input geometry.

	Geometry deformation block used: FFD (shape sphere).

	<b>To run</b>: ./parallel_example_00002 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00002() {

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry();
    mimmo0->setIOMode(IOMode::READ);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");

    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry();
    mimmo1->setIOMode(IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::SURFVTU);
    mimmo1->setWriteFilename("parallel_output_00002.0001");

    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry();
    mimmo2->setIOMode(IOMode::WRITE);
    mimmo2->setWriteDir(".");
    mimmo2->setWriteFileType(FileType::SURFVTU);
    mimmo2->setWriteFilename("parallel_output_00002.0002");

    /* Instantiation of a Partition object with default patition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition* partition = new mimmo::Partition();
    partition->setPlotInExecution(true);

    /* Instantiation of a FFDobject with spherical shape.
     * Setup of origin and span (radius and span angles) of sphere.
     * Plot Optional results during execution active for FFD block.
     */
    mimmo::FFDLattice* lattice = new mimmo::FFDLattice();
    darray3E origin = {0.0, 0.0,0.0};
    darray3E span;
    span[0]= 3.01;
    span[1]= 2*M_PI;
    span[2]= M_PI;

    /* Set number of nodes of the mesh (dim) and degree of nurbs functions (deg).
     */
    iarray3E dim, deg;
    dim[0] = 30;
    dim[1] = 30;
    dim[2] = 30;
    deg[0] = 2;
    deg[1] = 2;
    deg[2] = 2;

    lattice->setLattice(origin,span,mimmo::ShapeType::SPHERE,dim, deg);

    /* Change reference system to work in local spherical coordinates.
     * Set coordinates as CLAMPED (continuity in origins of angles).
     */
    lattice->setRefSystem(2, darray3E{0,1,0});
    lattice->setCoordType(mimmo::CoordType::CLAMPED, 2);
    lattice->setPlotInExecution(false);

    /* Build mesh of lattice outside the execution chain
     * to use it during setup the displacements.
     */
    lattice->build();

    /* Creation of Generic input block and set it with the
     * displacements of the control nodes of the lattice.
     * Use random values to set the displacements of the control nodes at a longitude
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
        if(l1 > 0 && lattice->getLocalPoint(l1,l2,l3)[1] < M_PI){
            displ[i][0] = 1.0*( (double) (rand()) / RAND_MAX );
        }
        if( (l1 > 0 && lattice->getLocalPoint(l1,l2,l3)[1] >= M_PI)
                || lattice->getLocalPoint(l1,l2,l3)[1] == 0){
            displ[i][0] = 1.25;
        }

    }

    /* Set Generic input block with the
     * displacements defined above.
     */
    mimmo::GenericInput* input = new mimmo::GenericInput();
    input->setInput(displ);

    /* Create applier block.
     * It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /* Setup pin connections.
     */
    mimmo::pin::addPin(mimmo0, partition, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, mimmo1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, lattice, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(input, lattice, M_DISPLS, M_DISPLS);
    mimmo::pin::addPin(lattice, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applier, mimmo2, M_GEOM, M_GEOM);

    /* Setup execution chain.
     * The object can be insert in the chain in random order.
     * The chain object recover the correct order of execution from
     * the pin connections.
     */
    mimmo::Chain ch0;
    ch0.addObject(input);
    ch0.addObject(partition);
    ch0.addObject(applier);
    ch0.addObject(lattice);
    ch0.addObject(mimmo0);
    ch0.addObject(mimmo1);
    ch0.addObject(mimmo2);

    /* Execution of chain.
     * Use debug flag false to avoid printing intermediate results of the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
    delete lattice;
    delete applier;
    delete input;
    delete mimmo0;
    delete partition;
    delete mimmo1;
    delete mimmo2;

    lattice = NULL;
    applier = NULL;
    input   = NULL;
    mimmo0 = NULL;
    partition = NULL;
    mimmo1 = NULL;
    mimmo2 = NULL;

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
            /**<Calling mimmo Test routine*/
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
