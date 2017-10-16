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

using namespace mimmo;

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
     * Extraction of PID = 1,2 (Wing wall and symmetry plane boundaries
     * in input file).
     */
    CGNSPidExtractor * cgnsExtr2 = new CGNSPidExtractor();
    cgnsExtr2->setPID({1, 2});
    cgnsExtr2->setForcedToTriangulate(false);
    cgnsExtr2->setPlotInExecution(true);

    /* Instantiation of a Selection By Box block.
     * Setup of span and origin of cube.
     */
    SelectionByBox      * boxSel = new SelectionByBox();
    boxSel->setOrigin({{700., 800., 0.}});
    boxSel->setSpan(1500.,1600.,100.);
    boxSel->setPlotInExecution(true);

    /* Creation of rotation block.
     */
    RotationGeometry* rotation = new RotationGeometry();
    rotation->setDirection(darray3E{0.,1.,0.});
    rotation->setRotation((M_PI/6));

    /* Create reconstruct vector block and set to reconstruct over the whole
     * input surface geometry the displacements field.
     */
    ReconstructVector* recon = new ReconstructVector();

    /* Create propagate vector block and set to propagate over the whole
     * input volume geometry the displacements field.
     */
    PropagateVectorField* prop = new PropagateVectorField();
    prop->setWeightConstant(1.0);
//    prop->setSmoothingSteps(50);
    prop->setPlotInExecution(true);
    prop->setConvergence(true);
    prop->setTolerance(1.0e-03);
    prop->setDumpingFactor(1.0);
    prop->setDumpingRadius(3000);
    prop->setSolver(true);

    /* Create propagate vector block and set to propagate over the whole
     * input volume geometry the displacements field.
     */
    ExtractVectorField* extrF = new ExtractVectorField();
    extrF->setMode(1);
    extrF->setPlotInExecution(true);

    /* Create applier block.
     * It applies the deformation displacements
     * to the selected input volume geometry.
     */
    Apply* applier = new Apply();

    /* Create applier block.
     * It applies the deformation displacements
     * to the selected input surface geometry.
     */
    Apply* applierS = new Apply();

    /* Create PINs. */
    addPin(cgnsI, cgnsExtr, M_GEOM2, M_GEOM);
    addPin(cgnsExtr2, boxSel, M_GEOM, M_GEOM);
    addPin(boxSel, rotation, M_GEOM, M_GEOM);
    addPin(rotation, recon, M_GDISPLS, M_VECTORFIELD);
    addPin(cgnsI, cgnsExtr2, M_GEOM2, M_GEOM);
    addPin(cgnsExtr2, recon, M_GEOM, M_GEOM);
    addPin(cgnsI, prop, M_GEOM, M_GEOM);
    addPin(cgnsExtr2, prop, M_GEOM, M_GEOM2);
    addPin(recon, prop, M_VECTORFIELD, M_GDISPLS);
    addPin(prop, applier, M_GDISPLS, M_GDISPLS);
    addPin(cgnsI, applier, M_GEOM, M_GEOM);

    addPin(cgnsExtr, extrF, M_GEOM, M_GEOM);
    addPin(prop, extrF, M_GDISPLS, M_VECTORFIELD);
    addPin(extrF, applierS, M_VECTORFIELD, M_GDISPLS);
    addPin(cgnsExtr, applierS, M_GEOM, M_GEOM);

    addPin(applier, cgnsO, M_GEOM, M_GEOM);
    addPin(applierS, cgnsO, M_GEOM, M_GEOM2);
    addPin(cgnsI, cgnsO, M_BCCGNS, M_BCCGNS);

    /* Create and execute chain. */
    Chain ch0;
    ch0.addObject(cgnsI);
    ch0.addObject(cgnsExtr);
    ch0.addObject(cgnsExtr2);
    ch0.addObject(boxSel);
    ch0.addObject(rotation);
    ch0.addObject(recon);
    ch0.addObject(prop);
    ch0.addObject(extrF);
    ch0.addObject(applier);
    ch0.addObject(applierS);
    ch0.addObject(cgnsO);
    ch0.exec(true);

    /* Destroy objects. */
    delete cgnsI;
    delete cgnsExtr;
    delete cgnsExtr2;
    delete cgnsO;
    delete boxSel;
    delete rotation;
    delete recon;
    delete prop;
    delete extrF;
    delete applier;
    delete applierS;

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
        /**< Call mimmo example routine. */
        example00001();

#if ENABLE_MPI==1
    }

    MPI::Finalize();
#endif

    return 0;
}

