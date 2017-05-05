/*---------------------------------------------------------------------------*\
 * 
 *  MIMIC
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License Commercial (//TODO Temporary header of license)
 *  This file is part of MIMIC.
 *
 *  MIMIC is a commercial software: you do not own rights to redistribute it 
 * 	and/or modify it both in source or pre-build formats
 *  Please contact Optimad offices for any further informations				
 *
 *  You should have received a copy of the mimic Commercial License
 *  along with MIMIC, as well as the key to unlock the software.
 *
 \ *----------------*-----------------------------------------------------------*/

#include "mimmo_iovtk.hpp"

using namespace mimmo;

// =================================================================================== //

void test00001() {


    /* Import vtk. */
    IOVTKScalar* iovtk = new IOVTKScalar();
    iovtk->setRead(true);
    iovtk->setReadDir("geodata");
    iovtk->setReadFilename("iovtk");
    iovtk->setScaling(1.0);
    iovtk->setNormalize(false);

    /* Export vtk with polydata */
    IOVTKScalar* iovtk2 = new IOVTKScalar();
    iovtk2->setWrite(true);
    iovtk2->setWriteDir(".");
    iovtk2->setWriteFilename("iovtk_output_00001");

    /* Export vtk without polydata */
    IOVTKScalar* iovtk3 = new IOVTKScalar();
    iovtk3->setWrite(true);
    iovtk3->setWriteDir(".");
    iovtk3->setWriteFilename("iovtk_output_00002");

    /* Create PINs. */
    addPin(iovtk,iovtk2, PortType::M_GEOM, PortType::M_GEOM);
    addPin(iovtk,iovtk2, PortType::M_SCALARFIELD, PortType::M_SCALARFIELD);
    addPin(iovtk,iovtk2, PortType::M_POLYDATA_, PortType::M_POLYDATA_);
    addPin(iovtk,iovtk3, PortType::M_GEOM, PortType::M_GEOM);
    addPin(iovtk,iovtk3, PortType::M_SCALARFIELD, PortType::M_SCALARFIELD);

    /* Create and execute chain. */
    Chain ch0;
    ch0.addObject(iovtk);
    ch0.addObject(iovtk2);
    ch0.addObject(iovtk3);
    ch0.exec(true);

    /* Destroy objects. */
    delete iovtk;
    delete iovtk2;
    delete iovtk3;
    iovtk   = NULL;
    iovtk2  = NULL;
    iovtk3  = NULL;

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
        /**<Calling mimmo Test routines*/

        test00001() ;

#if ENABLE_MPI==1
    }

    MPI::Finalize();
#endif

    return 0;

}
