/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#include "Lattice.hpp"
#include <bitpit_common.hpp>

// =================================================================================== //
/*!
	\example core_example_00002.cpp

	\brief Example of Lattice creation.

	Using: Lattice.

	<b>To run</b>              : ./core_example_00002 \n
    <b>To run (MPI version)</b>: mpirun -np X core_example_00002 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //

void test00002() {

    /*
        Create a Lattice of cylindrical shape
    */
    mimmo::Lattice * latt3 = new mimmo::Lattice();
    latt3->setShape(mimmo::ShapeType::CYLINDER);
    latt3->setOrigin({{0.5, 0.0, 0.0}});
    latt3->setSpan({{0.75, 2.0*BITPIT_PI, 0.4}});
    latt3->setInfLimits({{0.25, 0.0, 0.0}});
    latt3->setRefSystem({{0,1,0}}, {{0,0,1}}, {{1,0,0}});
    latt3->setDimension(iarray3E({{2, 35, 2}}));
    latt3->setPlotInExecution(true);
    latt3->exec();

    /*
        Create a Lattice of wedge shape (prism with triangular basis)
    */
    mimmo::Lattice * lattw = new mimmo::Lattice();
    lattw->setShape(mimmo::ShapeType::WEDGE);
    lattw->setOrigin({{1.0, 0.0, 0.0}});
    lattw->setSpan({{1.0, 0.75, 0.4}});
    lattw->setRefSystem(2, {{0,1,0}});
    lattw->setDimension(iarray3E({{5, 5, 4}}));
    lattw->setPlotInExecution(true);
    lattw->exec();


    /* Clean up & exit;
     */
	delete lattw;
    delete latt3;

	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
		try{
            /**<Calling core function*/
            test00002();
        }
        catch(std::exception & e){
            std::cout<<"core_example_00002 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
