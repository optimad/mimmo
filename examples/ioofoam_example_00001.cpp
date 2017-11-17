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
 * \brief Example of reading/writing of a Openfoam case mesh manipulation.
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

    mimmo::SelectionByPID * pidsel = new mimmo::SelectionByPID();
    pidsel->setGeometry(reader->getBoundaryGeometry());
    pidsel->setPID(3);

    mimmo::OBBox * box = new mimmo::OBBox();
    box->setForceAABB(true);


    dvecarr3E displ(8, {{0.0,0.0,0.0}});
    displ[0][1] = 0.03;
    displ[1][1] = 0.03;
    displ[4][1] = 0.03;
    displ[5][1] = 0.03;
    std::array<int,3> dim = {{2,2,2}};
    mimmo::FFDLattice * manip = new mimmo::FFDLattice();
    manip->setShape(mimmo::ShapeType::CUBE);
    manip->setDimension(dim);
    manip->setDegrees(dim);
    manip->setDisplacements(displ);


    mimmo::PropagateVectorField * prop = new mimmo::PropagateVectorField();
    prop->setSolver(true);
//    prop->setSmoothingSteps(100);
    prop->setDumping(true);
    prop->setDumpingInnerDistance(0.08);
    prop->setDumpingOuterDistance(0.3);
    prop->setPlotInExecution(true);
    
    
    mimmo::Apply * applyBulk = new mimmo::Apply();
    mimmo::Apply * applyBoundary = new mimmo::Apply();

    mimmo::IOOFOAM * writer = new mimmo::IOOFOAM(IOOFMode::WRITEPOINTSONLY);
    writer->setDir(path);
    writer->setOverwrite(false);

    
    mimmo::pin::addPin(reader, pidsel, M_GEOM2, M_GEOM);
    mimmo::pin::addPin(pidsel, box, M_GEOM, M_GEOM );
    mimmo::pin::addPin(box, manip, M_POINT, M_POINT );
    mimmo::pin::addPin(box, manip, M_AXES, M_AXES );
    mimmo::pin::addPin(reader, manip, M_GEOM2, M_GEOM );
    
    mimmo::pin::addPin(reader, prop, M_GEOM, M_GEOM);
    mimmo::pin::addPin(reader, prop, M_GEOM2, M_GEOM2);
    mimmo::pin::addPin(manip, prop, M_GDISPLS, M_GDISPLS);
    
    mimmo::pin::addPin(reader, applyBulk, M_GEOM, M_GEOM);
    mimmo::pin::addPin(prop,   applyBulk, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(reader, applyBoundary, M_GEOM2, M_GEOM);
    mimmo::pin::addPin(manip,  applyBoundary, M_GDISPLS, M_GDISPLS);

    mimmo::pin::addPin(applyBulk, writer, M_GEOM, M_GEOM);
    mimmo::pin::addPin(applyBoundary, writer, M_GEOM, M_GEOM2);
    
    
    mimmo::Chain c0;
    c0.addObject(reader);
    c0.addObject(pidsel);
    c0.addObject(box);
    c0.setPlotDebugResults(true);
    c0.exec(true);
    
    reader->getBoundaryGeometry()->getPatch()->write("boundary");
    
    darray3E spanBox = box->getSpan();
    spanBox *= 1.1;
    spanBox[1] = 0.1;
    manip->setSpan(spanBox);
    
    //post process value of the box and create the manipulation chain
    mimmo::Chain c1;
    c1.addObject(manip);
    c1.addObject(prop);
    c1.addObject(applyBulk);
    c1.addObject(applyBoundary);
    c1.addObject(writer);

    c1.setPlotDebugResults(true);
    c1.exec(true);

    delete reader;
    delete pidsel;
    delete box;
    delete manip;
    delete prop;
    delete applyBulk;
    delete applyBoundary;
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
