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

// =================================================================================== //
/*!
	\example manipulators_example_00002.cpp

	\brief Example of usage of free form deformation block to manipulate an input geometry.

	Geometry deformation block used: FFDLattice (Shape->cube).

	<b>To run</b>: ./manipulators_example_00002 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00002() {

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry();

    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setBuildSkdTree(true);

    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("manipulators_output_00002.0000");

    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry();
    mimmo1->setIOMode(IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("manipulators_output_00002.0001");

    /* Instantiation of a FFDobject with default shape cube.
     * Setup of span and origin of cube.
     * Plot Optional results during execution active for FFD block.
     */
    mimmo::FFDLattice* lattice = new mimmo::FFDLattice();
    darray3E origin = {0.0, 0.0, 0.0};
    darray3E span;
    span[0]= 1.2;
    span[1]= 1.2;
    span[2]= 1.2;

    /* Set number of nodes of the mesh (dim) and degree of nurbs functions (deg).
     */
    iarray3E dim, deg;
    dim[0] = 20;
    dim[1] = 20;
    dim[2] = 20;
    deg[0] = 2;
    deg[1] = 2;
    deg[2] = 2;

    lattice->setLattice(origin, span, mimmo::ShapeType::CUBE, dim, deg);

    /* Creation of Generic input block to read the
     * displacements of the control nodes of the lattice.
     */
    mimmo::GenericInput* input = new mimmo::GenericInput();
    input->setReadFromFile(true);
    input->setReadDir("input");
    input->setFilename("manipulators_input_00002.txt");

    /* Create applier block.
     * It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /* Setup pin connections.
     */
    std::cout << " --- create pin ---" << std::endl;
    std::cout << " " << std::endl;
    /* Add pin with port TAG ONLY
     */

    std::cout << " add pin info : " << std::boolalpha << mimmo::pin::addPin(mimmo0, lattice, M_GEOM, M_GEOM) << std::endl;
    std::cout << " add pin info : " << std::boolalpha << mimmo::pin::addPin(input, lattice, M_DISPLS, M_DISPLS) << std::endl;
    std::cout << " add pin info : " << std::boolalpha << mimmo::pin::addPin(lattice, applier, M_GDISPLS, M_GDISPLS) << std::endl;
    std::cout << " add pin info : " << std::boolalpha << mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM) << std::endl;
    std::cout << " add pin info : " << std::boolalpha << mimmo::pin::addPin(applier, mimmo1, M_GEOM, M_GEOM) << std::endl;
    std::cout << " " << std::endl;

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(input);
    ch0.addObject(lattice);
    ch0.addObject(applier);
    ch0.addObject(mimmo1);

    //force the chain to plot all the optional results of its children...
    ch0.setPlotDebugResults(true);
    //...in the path specified by the User.
    ch0.setOutputDebugResults(".");

    /* Execution of chain.
     * Use debug flag true to full print out the execution steps.
     */
    std::cout << " " << std::endl;
    std::cout << " --- execution start ---" << std::endl;
    ch0.exec(true);
    std::cout << " --- execution done --- " << std::endl;
    std::cout << " " << std::endl;

    /* Clean up & exit;
     */
    delete lattice;
    delete applier;
    delete input;
    delete mimmo0;
    delete mimmo1;

    lattice = NULL;
    applier = NULL;
    input 	= NULL;
    mimmo0  = NULL;
    mimmo1  = NULL;

}

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        /**<Calling mimmo Test routine*/
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
