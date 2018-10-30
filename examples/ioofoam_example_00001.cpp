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
#include "mimmo_manipulators.hpp"
#include "openFoamFiles_native.hpp"
#include <exception>

// =================================================================================== //

/*!
 * \example ioofoam_example_00001.cpp
 * 
 * \brief Example of reading/writing of a OpenFOAM case mesh.
 * 
 * Mesh is read from an OpenFOAM case. A FFD deformation is applied.
 * The bulk volume mesh is deformed accordingly.
 * In writing, moved bulk points update those on the target mesh.
 *
 * Using: IOOFOAM
 * 
 * <b>To run</b>: ./ioofoam_example_00001 \n
 * 
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

void OFOAM_manip() {

    mimmo::IOOFOAM * reader = new mimmo::IOOFOAM(IOOFMode::READ);
    reader->setDir("geodata/OFOAM");

    mimmo::FFDLattice * ffd = new mimmo::FFDLattice();
    ffd->setShape(mimmo::ShapeType::CUBE);

    darray3E origin = {{0.367, 0., 0.}};
    darray3E span   = {{0.2, 0.1, 0.1}};
    dvecarr3E displ(12,{{0.0,0.0,0.0}});
    displ[4][1] = 0.045;
    displ[5][1] = 0.045;
    displ[6][1] = -0.045;
    displ[7][1] = -0.045;

    ffd->setOrigin(origin);
    ffd->setSpan(span);
    ffd->setDimension(ivector1D({{3,2,2}}));
    ffd->setDegrees(iarray3E({{2,2,2}}));
    ffd->setDisplacements(displ);
    ffd->setPlotInExecution(true);

    mimmo::Apply * applier = new mimmo::Apply();

    mimmo::IOOFOAM * writer = new mimmo::IOOFOAM(IOOFMode::WRITEPOINTSONLY);
    writer->setOverwrite(false);

    mimmo::pin::addPin(reader, ffd, M_GEOM, M_GEOM);
    mimmo::pin::addPin(reader, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(ffd, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applier, writer, M_GEOM, M_GEOM);
    
    mimmo::Chain c0;
    c0.addObject(reader);
    c0.addObject(ffd);
    c0.addObject(applier);
    c0.addObject(writer);
    c0.exec(true);

    delete reader;
    delete ffd;
    delete applier;
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
        try{
             OFOAM_manip() ;
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
