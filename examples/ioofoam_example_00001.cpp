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

#include "mimmo_iogeneric.hpp"
#include "IOOFOAM.hpp"
#include "mimmo_manipulators.hpp"

// =================================================================================== //

/*!
 * \example ioofoam_example_00001.cpp
 *
 * \brief Example of reading,morphing and writing of a OpenFOAM case mesh.
 *
 * PART1 - reading,morphing and writing of a OpenFOAM case mesh.

 * Mesh is read from an OpenFOAM case. A FFD deformation is applied.
 * The bulk volume mesh is deformed accordingly.
 * In writing, moved bulk points update those on the target mesh.

   PART 2 -  Example of reading field from a OpenFOAM case mesh

   Mesh and its boundary pressure scalar field are read from an OpenFOAM case.
   The applier get the scalar field and convert it in a vectorfield of
   geometrical displacements using the local boundary mesh normals.
   Then apply the displacements to the boundary mesh and save the deformed version in a vtu file.
 *
 * Using: IOOFOAM, FFDLattice, IOOFOAMScalarField, Apply

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

    mimmo::pin::addPin(reader, ffd, M_GEOMOFOAM, M_GEOM);
    mimmo::pin::addPin(reader, applier, M_GEOMOFOAM, M_GEOM);
    mimmo::pin::addPin(ffd, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applier, writer, M_GEOM, M_GEOMOFOAM);

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

void OFOAM_sensi() {

    mimmo::IOOFOAM * reader = new mimmo::IOOFOAM(IOOFMode::READ);
    reader->setDir("geodata/OFOAM");

    mimmo::IOOFOAMScalarField * fieldreader = new mimmo::IOOFOAMScalarField();
    fieldreader->setDir("geodata/OFOAM");
    fieldreader->setFieldName("p");

    mimmo::Apply * applier = new mimmo::Apply();
    applier->setScaling(0.1);

    mimmo::MimmoGeometry * writer = new mimmo::MimmoGeometry();
    writer->setIOMode(IOMode::WRITE);
    writer->setWriteDir(".");
    writer->setWriteFileType(FileType::SURFVTU);
    writer->setWriteFilename("ofoam_sensi_output");


    mimmo::pin::addPin(reader, fieldreader, M_GEOMOFOAM2, M_GEOMOFOAM2);
    mimmo::pin::addPin(reader, fieldreader, M_UMAPIDS, M_UMAPIDS);
    mimmo::pin::addPin(reader, applier, M_GEOMOFOAM2, M_GEOM);
    mimmo::pin::addPin(fieldreader, applier, M_SCALARFIELD2, M_SCALARFIELD);
    mimmo::pin::addPin(applier, writer, M_GEOM, M_GEOM);


    mimmo::Chain c0;
    c0.addObject(reader);
    c0.addObject(fieldreader);
    c0.addObject(applier);
    c0.addObject(writer);
    c0.exec(true);

    delete reader;
    delete fieldreader;
    delete applier;
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
             OFOAM_sensi() ;
        }
        catch(std::exception & e){
            std::cout<<"test_ioofoam_00001 PART1 exit with the following errors :"<<e.what()<<std::endl;
            return 1;
        }
        try{
             OFOAM_manip() ;
        }
        catch(std::exception & e){
            std::cout<<"test_ioofoam_00001 PART2 exit with the following errors :"<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
