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


#include "mimmo_core.hpp"
#include "mimmo_iogeneric.hpp"
#include "bitpit.hpp"
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmo::pin;

// =================================================================================== //
/*!
	\example core_example_00001.cpp

	\brief Example of interpolating data over an unstructured non-homogeneous mesh.

	Using: MimmoGeometry, MimmoPiercedVector

	<b>To run</b>: ./core_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00001() {

	/* Creation of mimmo containers.
	 */
	MimmoGeometry * mimmo0 = new MimmoGeometry();

	mimmo0->setIOMode(IOMode::READ);
	mimmo0->setReadDir("geodata");
	mimmo0->setReadFileType(FileType::STL);
	mimmo0->setReadFilename("plane3");

	mimmo0->execute();

	/*
	 * Creation of a synthetic point field.
	 */
	MimmoPiercedVector<double> pointField;
	pointField.initialize(mimmo0->getGeometry(), MPVLocation::POINT, 1.);
	std::array<double,3> center({{0.5,0.,0.}});
	for (bitpit::Vertex & vertex : mimmo0->getGeometry()->getPatch()->getVertices()){
		std::array<double,3> coords = vertex.getCoords();
		double value = norm2(coords-center);
		pointField[vertex.getId()] = value;
	}

	{
		std::vector<double> field(mimmo0->getGeometry()->getNVertices());
		int count = 0;
		for (bitpit::Vertex & vertex : mimmo0->getGeometry()->getPatch()->getVertices()){
			field[count] = pointField[vertex.getId()];
			count++;
		}
		mimmo0->getGeometry()->getPatch()->getVTK().addData("pointField", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, field);
		mimmo0->getGeometry()->getPatch()->write("core_example_00001.0000");
		mimmo0->getGeometry()->getPatch()->getVTK().removeData("pointField");
	}

	/*
	 * Interpolation of synthetic point field on cells.
	 */
	double p = 5.;
	MimmoPiercedVector<double> cellField = pointField.pointDataToCellData(p);
	{
		std::vector<double> field(mimmo0->getGeometry()->getNCells());
		int count = 0;
		for (bitpit::Cell & cell : mimmo0->getGeometry()->getPatch()->getCells()){
			field[count] = cellField[cell.getId()];
			count++;
		}
		mimmo0->getGeometry()->getPatch()->getVTK().addData("cellField", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, field);
		mimmo0->getGeometry()->getPatch()->write("core_example_00001.0001");
		mimmo0->getGeometry()->getPatch()->getVTK().removeData("cellField");
	}

	/*
	 * Back Interpolation of synthetic field on points.
	 */
	p = 5.;
	MimmoPiercedVector<double> pointField2 = cellField.cellDataToPointData(p);
	{
		std::vector<double> field(mimmo0->getGeometry()->getNVertices());
		int count = 0;
		for (bitpit::Vertex & vertex : mimmo0->getGeometry()->getPatch()->getVertices()){
			field[count] = pointField2[vertex.getId()];
			count++;
		}
		mimmo0->getGeometry()->getPatch()->getVTK().addData("pointField", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, field);
		mimmo0->getGeometry()->getPatch()->write("core_example_00001.0002");
		mimmo0->getGeometry()->getPatch()->getVTK().removeData("pointField");
	}

	/*
	 * Interpolation of synthetic point field on boundary interfaces.
	 */
	mimmo0->getGeometry()->buildInterfaces();
	p = 5.;
	MimmoPiercedVector<double> interfaceField = pointField.pointDataToBoundaryInterfaceData(p);
	{
		std::vector<double> field(mimmo0->getGeometry()->getPatch()->getInterfaceCount());
		int count = 0;
		for (bitpit::Interface & interface : mimmo0->getGeometry()->getPatch()->getInterfaces()){
			if (interface.isBorder()){
				field[count] = interfaceField[interface.getId()];
				log::cout() << "field on boundary interface " << count << " : " << field[count] << std::endl;
				count++;
			}
		}
	}

	/* Clean up & exit;
	 */
	delete mimmo0;
	return;

}


int main(int argc, char *argv[]) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc, &argv);

	{
#endif

		/**<Change the name of mimmo logger file (default mimmo.log)
		 * before initialization of BaseManipulation objects*/
		mimmo::setLogger("mimmo");

		/**<Calling mimmo Test routines*/
		try{
			test00001() ;
		}
		catch(std::exception & e){
			std::cout<<"utils_example_00001 exited with an error of type : "<<e.what()<<std::endl;
			return 1;
		}

#if BITPIT_ENABLE_MPI==1
	}

	MPI_Finalize();
#endif

	return (1);
}


