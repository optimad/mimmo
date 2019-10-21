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

#include "mimmo_utils.hpp"
#include "mimmo_iogeneric.hpp"

// =================================================================================== //
/*!
	\example utils_example_00001.cpp

	\brief Example of seeding a point cloud over an unstructured non-homogeneous mesh.

	Using: MimmoGeometry, CreateSeedsOnSurface

	<b>To run</b>: ./utils_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00001() {

    /* Creation of mimmo containers.
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry();

    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::SURFVTU);
    mimmo0->setReadFilename("mixedP2D");

    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::SURFVTU);
    mimmo0->setWriteFilename("utils_mesh_00001.0000");

    /*Creation of the seeder. Use Level set engine to place 9 points on
     * surface, starting from a seed point in the absolute origin.
     */
    mimmo::CreateSeedsOnSurface * cseed = new mimmo::CreateSeedsOnSurface();
    cseed->setSeed({{0.0,0.0,0.0}});
    cseed->setNPoints(9);
    cseed->setEngine(1);
    cseed->setPlotInExecution(true);

    /* Setup pin connections.
     */
    mimmo::pin::addPin(mimmo0, cseed, M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(cseed);

    /* Execution of chain.
     * Use debug flag false to avoid to print out the execution steps on console.
     */
    ch0.exec();

    /* Clean up & exit;
     */
    delete mimmo0;
    delete cseed;
    return;

}


int main(int argc, char *argv[]) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI==1
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

#if MIMMO_ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;
}
