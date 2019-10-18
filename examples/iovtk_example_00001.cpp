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
#include <exception>

using namespace mimmo;

/*!
 * \example iovtk_example_00001.cpp
 *
 * \brief Example of reading/writing a vtk polydata file w/ scalar field attached.
 *
 * Using: IOVTKScalar
 *
 * <b>To run</b>: ./iovtk_example_00001 \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

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
    pin::addPin(iovtk,iovtk2, M_GEOM, M_GEOM);
    pin::addPin(iovtk,iovtk2, M_SCALARFIELD, M_SCALARFIELD);
    pin::addPin(iovtk,iovtk2, M_POLYDATA_, M_POLYDATA_);
    pin::addPin(iovtk,iovtk3, M_GEOM, M_GEOM);
    pin::addPin(iovtk,iovtk3, M_SCALARFIELD, M_SCALARFIELD);

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

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        /**<Calling mimmo Test routines*/
        try{
            test00001() ;
        }
        catch(std::exception & e){
            std::cout<<"iovtk_example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;

}
