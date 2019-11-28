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

#include "mimmo_iogeneric.hpp"
#include "FFDLattice.hpp"

/*!
 * \example genericinput_example_00003.cpp
 *
 * \brief Example of reading displacements related to a 3x2x2 lattice.

 * Using: GenericDispls, FFDLattice
 *
 * <b>To run</b>: ./genericinput_example_00003 \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //


void test00003() {

    /* Reading plane
     */
	mimmo::MimmoGeometry * read = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    read->setReadDir("geodata");
    read->setReadFilename("plane4");
    read->setReadFileType(FileType::STL);
    read->setWriteDir(".");
    read->setWriteFilename("./genericinput_output_00003.0001");
    read->setWriteFileType(FileType::SURFVTU);

    /*!Creating Bezier FFD Lattice 3x2x2 around a box centered in 0,0,-0.5
     with span 2,2,1
    */
    mimmo::FFDLattice * latt = new mimmo::FFDLattice();
    latt->setName("genericinput_example_00003_Lattice");
    latt->setShape(mimmo::ShapeType::CUBE);
    latt->setOrigin({{0.0,0.0,-0.5}});
    latt->setSpan({{2.0,2.0,1.0}});
    latt->setDimension(iarray3E({{3,2,2}}));
    latt->setDegrees(iarray3E{2,1,1});
    latt->build();
    latt->setPlotInExecution(true);
    latt->setApply(true);

    /* Creation of GenericDispls block to read
       displacement from file
     */
    mimmo::GenericDispls * iodispls = new mimmo::GenericDispls(true);
    iodispls->setReadDir("input");
    iodispls->setReadFilename("generic_displ_00002.txt");

    /* Setup pin connections.
     */
    mimmo::pin::addPin(read, latt, M_GEOM, M_GEOM);
    mimmo::pin::addPin(iodispls, latt, M_DISPLS, M_DISPLS);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(read);
    ch0.addObject(latt);
    ch0.addObject(iodispls);


    /* Execution of chain.
     * Use debug flag true to print out the execution steps.
     */
    ch0.exec(true);

    /*
        print the deformed geometry
    */
    read->getGeometry()->getPatch()->write("./genericinput_output_00003.0002");

    /* Clean up & exit;
     */
    delete read;
    delete iodispls;
    delete latt;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        /**<Calling mimmo Test routines*/
        try{
            test00003();
        }
        catch(std::exception & e){
            std::cout<<"genericinput_example_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;
}
