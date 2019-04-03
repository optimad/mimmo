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
#include "mimmo_iocgns.hpp"
#include "mimmo_geohandlers.hpp"
#include "mimmo_manipulators.hpp"
#include "mimmo_propagators.hpp"
#include <exception>

using namespace mimmo;

/*!
 * \example iocgns_example_00001.cpp
 *
 * \brief Reading of a CGNS volume mesh and boundary extraction.
 *
 * Using: IOCGNS, CGNSPidExtractor, SelectionByBox
 *
 * Depends on mimmo optional module geohandlers
 *
 * <b>To run</b>: ./iocgns_example_00001 \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //

void example00001() {

    /* Create IO_CGNS object to import input file. */
    IOCGNS * cgnsI = new IOCGNS();
    cgnsI->setRead(true);
    cgnsI->setReadDir("geodata");
    cgnsI->setReadFilename("grid");

    /* Create IO_CGNS object to export output file. */
    IOCGNS * cgnsO = new IOCGNS();
    cgnsO->setRead(false);
    cgnsO->setWriteDir(".");
    cgnsO->setWriteFilename("iocgns_output_00001");

    /* Create CGNS PID extractor object to test input file.
     * Extraction of PID = 1,2,3 (All boundaries in input file).
     */
    CGNSPidExtractor * cgnsExtr = new CGNSPidExtractor();
    cgnsExtr->setPID({1, 2, 3});
    cgnsExtr->setForcedToTriangulate(false);
    cgnsExtr->setPlotInExecution(true);

    /* Create CGNS PID extractor object to test input file.
     * Extraction of PID = 1,2 (Wing wall and outer part boundaries
     * where imposing Dirichlet conditions).
     */
    CGNSPidExtractor * cgnsDirichlet = new CGNSPidExtractor();
    cgnsDirichlet->setPID({1, 2});
    cgnsDirichlet->setForcedToTriangulate(false);
    cgnsDirichlet->setPlotInExecution(true);

    /* Create CGNS PID extractor object to test input file.
     * Extraction of PID = 3 (Simmetry plane
     * where imposing Slip/impermeability conditions).
     */
    CGNSPidExtractor * cgnsSlip = new CGNSPidExtractor();
    cgnsSlip->setPID({3});
    cgnsSlip->setForcedToTriangulate(false);
    cgnsSlip->setPlotInExecution(true);

    /* Instantiation of a Selection By Box block.
     * Setup of span and origin of cube.
     */
    SelectionByBox      * boxSel = new SelectionByBox();
    boxSel->setOrigin({{700., 800., 0.}});
    boxSel->setSpan(1500.,1600.,100.);
    boxSel->setPlotInExecution(true);


    /* Create PINs. */
    addPin(cgnsI, cgnsExtr, M_GEOM2, M_GEOM);
    addPin(cgnsI, cgnsDirichlet, M_GEOM2, M_GEOM);
    addPin(cgnsI, cgnsSlip, M_GEOM2, M_GEOM);

    addPin(cgnsDirichlet, boxSel, M_GEOM, M_GEOM);


    /* Create and execute chain. */
    Chain ch0;
    ch0.addObject(cgnsI);
    ch0.addObject(cgnsExtr);
    ch0.addObject(cgnsDirichlet);
    ch0.addObject(cgnsSlip);
    ch0.addObject(boxSel);


    ch0.exec(true);

    cgnsI->getGeometry()->getPatch()->write("bulk");
    cgnsI->getSurfaceBoundary()->getPatch()->write("boundary");

    /* Destroy objects. */
    delete cgnsI;
    delete cgnsExtr;
    delete cgnsDirichlet;
    delete cgnsSlip;
    delete boxSel;

    return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
    MPI::Init(argc, argv);

    {
#endif
        try{
            /**< Call mimmo example routine. */
            example00001();
        }
        catch(std::exception & e){
            std::cout<<"iocgns_example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if ENABLE_MPI==1
    }

    MPI::Finalize();
#endif

    return 0;
}
