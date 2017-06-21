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
#include "bitpit.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmo::pin;

// =================================================================================== //
/*!
	\example manipulators_example_00003.cpp

	\brief Example of usage of free form deformation block to manipulate an input geometry.

	Geometry deformation block used: FFD (shape cylinder).

	<b>To run</b>: ./manipulators_example_00003 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void test00003() {

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
    MimmoGeometry * mimmo0 = new MimmoGeometry();
    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("catpipe");
    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("manipulators_output_00003.0000");

    MimmoGeometry * mimmo1 = new MimmoGeometry();
    mimmo1->setIOMode(IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("manipulators_output_00003.0001");

    /* Instantiation of a FFDobject with cylindrical shape.
     * Setup of span and origin of cylinder.
     * Plot Optional results during execution active for FFD block.
     */
    FFDLattice* lattice = new FFDLattice();
    darray3E origin = {-1537.5, -500.0, 3352.5};
    darray3E span;
    span[0]= 100.0;
    span[1]= 2*M_PI;
    span[2]= 1000.0;

    /* Set number of nodes of the mesh (dim) and degree of nurbs functions (deg).
     */
    iarray3E dim, deg;
    dim[0] = 2;
    dim[1] = 15;
    dim[2] = 20;

    deg[0] = 1;
    deg[1] = 2;
    deg[2] = 10;

    lattice->setLattice(origin,span,ShapeType::CYLINDER,dim, deg);

    /* Change reference system to work in local cylindrical coordinates.
     */
    lattice->setRefSystem(2, darray3E{0,-1,0});
    lattice->setDisplGlobal(false);
    lattice->setPlotInExecution(true);

    /* Build mesh of lattice outside the execution chain
     * to use it during setup the displacements.
     */
    lattice->build();

    /* Creation of Generic input block and set it with the
     * displacements of the control nodes of the lattice.
     * Use a polynomial law in local coordinates to define the displacements.
     * The used expression guarantees the continuity of the surface at the interface
     * between deformed and the undeformed parts.
     */
    int ndof = lattice->getNNodes();
    dvecarr3E displ(ndof, darray3E{0,0,0});
    int lt = 2*dim[2]/3;
    double lx;
    double a, b, c, max;
    c = 6;
    b = -15;
    a = 10;
    max = 100;
    for (int i=0; i<ndof; i++){
        int l1,l2,l3;
        int index = lattice->accessGridFromDOF(i);
        lattice->accessPointIndex(index,l1,l2,l3);
        if (l3 < lt && l1 > 0){
            lx = double(l3)/double(lt);
            displ[i][0] = max*(c*pow(lx,5) + b*pow(lx,4) + a*pow(lx,3));
        }
        else if (l1 > 0){
            displ[i][0] = max;
        }
    }

    /* Set Generic input block with the
     * displacements defined above.
     */
    GenericInput* input = new GenericInput();
    input->setReadFromFile(false);
    input->setInput(displ);

    /* Set Generic output block to write the
     * displacements defined above.
     */
    GenericOutput * output = new GenericOutput();
    output->setFilename("manipulators_output_00003.csv");
    output->setCSV(true);

    /* Create applier block.
     * It applies the deformation displacements to the original input geometry.
     */
    Apply* applier = new Apply();

    /* Setup pin connections.
     */
    pin::addPin(mimmo0, lattice, PortType::M_GEOM, PortType::M_GEOM);
    pin::addPin(input, lattice, PortType::M_DISPLS, PortType::M_DISPLS);
    pin::addPin(input, output, PortType::M_DISPLS, PortType::M_DISPLS);
    pin::addPin(mimmo0, applier, PortType::M_GEOM, PortType::M_GEOM);
    pin::addPin(lattice, applier, PortType::M_GDISPLS, PortType::M_GDISPLS);
    pin::addPin(applier, mimmo1, PortType::M_GEOM, PortType::M_GEOM);

    /* Setup execution chain.
     */
    Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(input);
    ch0.addObject(output);
    ch0.addObject(lattice);
    ch0.addObject(applier);
    ch0.addObject(mimmo1);

    /* Execution of chain.
     * Use debug flag false (default) to avoid to to print out the execution steps.
     */
    cout << " " << endl;
    cout << " --- execution start --- " << endl;
    ch0.exec();
    cout << " --- execution done --- " << endl;
    cout << " " << endl;

    /* Clean up & exit;
     */
    delete applier;
    delete lattice;
    delete input;
    delete output;
    delete mimmo0;
    delete mimmo1;

    lattice = NULL;
    applier = NULL;
    input   = NULL;
    output   = NULL;
    mimmo0  = NULL;
    mimmo1  = NULL;

    return;
}

int	main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
    MPI::Init(argc, argv);

    {
#endif
        /**<Calling mimmo Test routine*/
        test00003() ;

#if ENABLE_MPI==1
    }

    MPI::Finalize();
#endif

    return(1);
}
