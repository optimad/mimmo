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
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

// =================================================================================== //
/*!
	\example manipulators_example_00002.cpp

	\brief Example of usage of free form deformation Lattice to manipulate an input geometry.

	Using: MimmoGeometry, FFDLattice , GenericInput, Apply, Chain, Partition(MPI version).

	<b>To run</b>              : ./manipulators_example_00002 \n
    <b>To run (MPI version)</b>: mpirun -np X manipulators_example_00002 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00002() {

    /*
        Read a sphere from STL file. Convert mode is to save the just read geometry in
        another file with name manipulators_output_00002.0000.stl
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setBuildSkdTree(true);

    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("manipulators_output_00002.0000");

    /*
        Write final deformed mesh on file
     */
    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("manipulators_output_00002.0001");

#if MIMMO_ENABLE_MPI
    /*
        Distribute target mesh among processes.
     */
    mimmo::Partition* partition= new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
#endif

    /*
       Create a FFDLattice manipulator, shaped as a box.
       Span and origin of the box are required.
       Dimensions and nurbs degrees for each spatial directions must be provided.
       Plot Optional results during execution active for FFD block.
     */
    mimmo::FFDLattice* lattice = new mimmo::FFDLattice();
    darray3E origin = {0.0, 0.0, 0.0};
    darray3E span;
    span[0]= 1.2;
    span[1]= 1.2;
    span[2]= 1.2;

    /*
       Set number of nodes of the mesh (dim) and degree of nurbs functions (deg).
     */
    iarray3E dim, deg;
    dim[0] = 20;
    dim[1] = 20;
    dim[2] = 20;
    deg[0] = 2;
    deg[1] = 2;
    deg[2] = 2;

    lattice->setLattice(origin, span, mimmo::ShapeType::CUBE, dim, deg);

    /*
        Reading displacements associated ot the lattice's nodes from external
        plain file.
     */
    mimmo::GenericInput* input = new mimmo::GenericInput();
    input->setReadFromFile(true);
    input->setReadDir("input");
    input->setFilename("manipulators_input_00002.txt");

    /*
        It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /*
        Setup pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, partition, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, lattice, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(mimmo0, lattice, M_GEOM, M_GEOM);
#endif
    mimmo::pin::addPin(input, lattice, M_DISPLS, M_DISPLS);
    mimmo::pin::addPin(lattice, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(applier, mimmo1, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition);
#endif
    ch0.addObject(input);
    ch0.addObject(lattice);
    ch0.addObject(applier);
    ch0.addObject(mimmo1);

    //force the chain to plot all the optional results of its children...
    ch0.setPlotDebugResults(true);
    //...in the path specified by the User.
    ch0.setOutputDebugResults(".");

    /*
        Execute the chain.
     * Use debug flag true to full print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
#if MIMMO_ENABLE_MPI
    delete partition;
#endif
    delete lattice;
    delete applier;
    delete input;
    delete mimmo0;
    delete mimmo1;

}

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        /**<Calling core function*/
        try{
            test00002() ;
        }
        catch(std::exception & e){
            std::cout<<"manipulators_example_00002 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return  0;
}
