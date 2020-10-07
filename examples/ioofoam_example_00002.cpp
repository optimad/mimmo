/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2020 OPTIMAD engineering Srl
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
#include "IOOFOAM.hpp"
#include "mimmo_manipulators.hpp"

// =================================================================================== //

/*!
 * \example ioofoam_example_00002.cpp
 *
 * \brief Example of reading,morphing and writing of a OpenFOAM case mesh.
 *
 * PART1 - reading and writing of a parallel OpenFOAM case mesh.

 * Mesh is read from an OpenFOAM case.
 * In writing, boundary mesh is written on output.
 *
 * Using: IOOFOAM

 * <b>To run</b>: ./ioofoam_example_00002 \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void OFOAM_readwrite() {

    mimmo::IOOFOAM * reader = new mimmo::IOOFOAM(false);
    reader->setDir("geodata/ofoam.00002");

    mimmo::MimmoGeometry * surfacewriter = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    surfacewriter->setWriteDir(".");
    surfacewriter->setWriteFileType(FileType::SURFVTU);
    surfacewriter->setWriteFilename("ioofoam.00002");

    mimmo::MimmoGeometry * volumewriter = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    volumewriter->setWriteDir(".");
    volumewriter->setWriteFileType(FileType::VOLVTU);
    volumewriter->setWriteFilename("ioofoam.00003");

    mimmo::IOOFOAM * writer = new mimmo::IOOFOAM(true);
    writer->setDir("geodata/ofoam.00002");

    mimmo::pin::addPin(reader, writer, M_GEOMOFOAM2, M_GEOMOFOAM2);
    mimmo::pin::addPin(reader, surfacewriter, M_GEOMOFOAM2, M_GEOM);
    mimmo::pin::addPin(reader, volumewriter, M_GEOMOFOAM, M_GEOM);
    mimmo::pin::addPin(reader, writer, M_GEOMOFOAM, M_GEOMOFOAM);
    mimmo::pin::addPin(reader, writer, M_UMAPIDS, M_UMAPIDS);


    mimmo::Chain c0;
    c0.addObject(reader);
    c0.addObject(surfacewriter);
    c0.addObject(volumewriter);
    c0.addObject(writer);
    c0.exec(true);

    delete reader;
    delete surfacewriter;
    delete volumewriter;
    delete writer;
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
             OFOAM_readwrite() ;
        }
        catch(std::exception & e){
            std::cout<<"ioofoam_example_00002 exit with the following errors :"<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
