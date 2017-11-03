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
#include <exception>

// =================================================================================== //

int test1() {
//     //char * rootPath ("./ascii");
//     Foam::Time *foamRunTime = 0;
//     Foam::fvMesh *foamMesh = 0;
//     
//     mimmo::foamUtilsNative::initializeCase("/home/rocco/Desktop/exampleOpenFoamMeshes/ascii", &foamRunTime, &foamMesh);
//     
//     std::cout<<"Number of boundary patches : "<<foamMesh->boundaryMesh().size()<<std::endl;
//     
//     Foam::label nCells = foamMesh->nCells();
//     std::cout<<"Number of cells :  " <<nCells<<std::endl;
//     auto cellShapes =foamMesh->cellShapes(); 
//     auto cells = foamMesh->cells();
//     auto faces = foamMesh->faces();
//     for(Foam::label i=0; i<nCells; ++i){
//         std::cout<<"Cell "<<i<<"  ================================================="<<std::endl;
//         cellShapes[i].collapse();
//         std::cout<<"model     :"<<'\t'<<cellShapes[i].model().name()<<std::endl;
//         std::cout<<"nPoints   :"<<'\t'<<cellShapes[i].nPoints()<<std::endl;
//         std::cout<<"nFaces    :"<<'\t'<<cellShapes[i].nFaces()<<std::endl;
//         std::cout<<"nEdges    :"<<'\t'<<cellShapes[i].nEdges()<<std::endl;
//         std::cout<<"faces:"<<std::endl;
//         auto listlabelFaces = cellShapes[i].meshFaces(faces, cells[i]);
//         std::cout<<"ORDERED MODEL faces:      "<<'\t';
//         for(Foam::List<Foam::label>::iterator it=listlabelFaces.begin(); it!=listlabelFaces.end(); ++it){
//             std::cout<<*it<<'\t';
//         }
//         std::cout<<std::endl;
//         std::cout<<"RAW FROM CELL faces:      "<<'\t';
//         Foam::label nFaces = cells[i].size();
//         auto listCellFaces = cells[i];
//         for(Foam::label j=0; j<nFaces; ++j ){
//             std::cout<<listCellFaces[j]<<'\t';
//         }
//         std::cout<<std::endl;
//         
//         std::cout<<"RAW FROM CELLSHAPE PP   : "<<'\t';
//         Foam::label sizewhatCellShapes = cellShapes[i].size();
//         auto whatCellShapes = cellShapes[i];
//         for(Foam::label j=0; j<sizewhatCellShapes; ++j ){
//             std::cout<<whatCellShapes[j]<<'\t';
//         }
//         std::cout<<std::endl;
//         
//         std::cout<<"RAW FROM CELL      PP   : "<<'\t';
//         auto pplabels = cells[i].labels(faces);
//         Foam::label pp = pplabels.size();
//         for(Foam::label j=0; j<pp; ++j ){
//             std::cout<<pplabels[j]<<'\t';
//         }
//         std::cout<<std::endl;
//         
//     /*    
//         for (const auto & idface: cells[i]){
//             std::cout<<"vertex-face connectivity:"<<'\t';
//             
//             for(Foam::label j =0; j<faces[id].size(); ++j){
//                 std::cout<<face[j]<<'\t';
//             }
//             std::endl;
//         }*/
//         exit(1);
//         std::cout<<"=============================================================="<<std::endl;
//     }
    

    mimmo::IOOFOAM * reader = new mimmo::IOOFOAM();
    reader->setRead(true);
    reader->setReadDir("/home/rocco/Desktop/exampleOpenFoamMeshes/ascii");
    
    reader->execute();
    
    delete reader;
    
    return  0;
    
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

        int val = 1;
        try{
            val = test1() ;
        }
        catch(std::exception & e){
            std::cout<<"test_ioofoam_00001 exit with the following errors :"<<e.what()<<std::endl;
            return 1;
        }

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif

	return val;
}
