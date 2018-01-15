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

#include "IOOFOAM.hpp"
#include "MeshSelection.hpp"
#include "ReconstructFields.hpp"
#include "mimmo_manipulators.hpp"
#include "mimmo_utils.hpp"
#include "openFoamFiles_native.hpp"
#include <exception>

// =================================================================================== //

/*!
 * \example ioofoam_example_00001.cpp
 * 
 * \brief Example of reading/writing of a Openfoam case mesh.
 * 
 * Mesh is read from a path given by the User. A RBF point is created on a
 * randomly picked-patch boundary of the mesh and moved normally to the patch itself.
 * The bulk volume mesh is deformed accordingly.
 * In writing, moved bulk points update those on the target mesh.
 *
 * Using: IOOFOAM, FFDLattice, PropagateField, ReconstructVector, SelectionByPID, Apply, OBBox
 * Depends on geohandlers, utils
 * 
 * <b>To run</b>: ./ioofoam_example_00001  < path to OpenFoam case > \n
 * 
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void OFOAM_manip(std::string & path) {

    //First Chain;
    
    mimmo::IOOFOAM * reader = new mimmo::IOOFOAM(IOOFMode::READ);
    reader->setDir(path);

    mimmo::IOOFOAM * writer = new mimmo::IOOFOAM(IOOFMode::WRITEPOINTSONLY);
    writer->setDir(path);
    writer->setOverwrite(false);

    mimmo::pin::addPin(reader, writer, M_GEOM, M_GEOM);
    mimmo::pin::addPin(reader, writer, M_GEOM, M_GEOM2);
    
    mimmo::Chain c0;
    c0.addObject(reader);
    c0.addObject(writer);
    c0.exec(true);

    delete reader;
    delete writer;
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
        if (argc < 2) {
            std::cout<<"please insert a valid path to OpenFoam case to manipulate the mesh"<<std::endl;
            return 1;
        }
        if(argv[1] == "--help"){
            std::cout<<"Insert after executable command a valid path to OpenFoam case to manipulate the mesh"<<std::endl;
            return 1;
        }
        try{
            std::string path(argv[1]);
             OFOAM_manip(path) ;
        }
        catch(std::exception & e){
            std::cout<<"test_ioofoam_00001 exit with the following errors :"<<e.what()<<std::endl;
            return 1;
        }

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif

	return 0;
}
