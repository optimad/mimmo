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
    addPin(iovtk,iovtk2, M_GEOM, M_GEOM);
    addPin(iovtk,iovtk2, M_SCALARFIELD, M_SCALARFIELD);
    addPin(iovtk,iovtk2, M_POLYDATA_, M_POLYDATA_);
    addPin(iovtk,iovtk3, M_GEOM, M_GEOM);
    addPin(iovtk,iovtk3, M_SCALARFIELD, M_SCALARFIELD);

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
