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
#include <exception>


// =================================================================================== //
/*!
 * //testing creation of PointCloud and 3DCurve from a raw list of inputs
 */
int test1() {

    mimmo::GenericInput * readpoints = new mimmo::GenericInput("./input", "rawpoints.csv", true);
    readpoints->exec();

    int np = (int)readpoints->getResult<dvecarr3E>().size();
    MPI_Comm communicator;
    MPI_Comm_dup(MPI_COMM_WORLD, &communicator);
    MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_INT, MPI_MAX, communicator);

    mimmo::GenericInput * readfield = new mimmo::GenericInput("./input", "rawfield.csv", true);
    readfield->exec();

    mimmo::CreatePointCloud * pc = new mimmo::CreatePointCloud();
    pc->setName("mimmo_iogeneric_parallel_00000_PointCloud");
    pc->setRawPoints(readpoints->getResult<dvecarr3E>());
    pc->setRawVectorField(readfield->getResult<dvecarr3E>());
    pc->setPlotInExecution(true);
    pc->exec();

    int refCells = 40;
    mimmo::Create3DCurve * curve = new mimmo::Create3DCurve();
    curve->setName("mimmo_iogeneric_parallel_00000_3DCurve");
    curve->setRawPoints(readpoints->getResult<dvecarr3E>());
    curve->setRawVectorField(readfield->getResult<dvecarr3E>());
    curve->setClosedLoop(true);
    curve->setNCells(refCells);
    curve->setPlotInExecution(true);
    curve->exec();

    bool check = (pc->getGeometry()->getNGlobalVertices() == np);
    check = check && (curve->getGeometry()->getNGlobalCells() == refCells);
  

    delete readpoints;
    delete readfield;
    delete pc;
    delete curve;

    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
    int val = 1;
    try{
        /**<Calling mimmo Test routines*/
        val = test1() ;
    }
    catch(std::exception & e){
        std::cout<<"test_iogeneric_parallel_00000 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
