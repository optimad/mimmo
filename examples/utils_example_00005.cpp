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
#if MIMMO_ENABLE_MPI
    #include "Partition.hpp"
#endif

// =================================================================================== //
/*!
	\example utils_example_00005.cpp

	\brief Project external curve on the Stanford Bunny surface.

	Using: MimmoGeometry, ProjPatchOnSurface

	<b>To run</b>: ./utils_example_00005 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00005() {

    /* Reading target surface (bunny)
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    bitpit::Logger & log = mimmo0->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);

    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("stanfordBunny3");

    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::SURFVTU);
    mimmo0->setWriteFilename("utils_mesh_00005_surface");
    mimmo0->execute();

    /* Reading 3D Curve to be projected
     */
	mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);

    mimmo1->setReadDir("geodata");
    mimmo1->setReadFileType(FileType::CURVEVTU);
    mimmo1->setReadFilename("curve");

    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::CURVEVTU);
    mimmo1->setWriteFilename("utils_mesh_00005_curve");
    mimmo1->execute();


    log<<"read geometries from file"<<std::endl;

    mimmo::MimmoSharedPointer<mimmo::MimmoObject> curve3D;
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> targetSurf;

#if MIMMO_ENABLE_MPI
    //partition of bunny and 3D curve
    mimmo::Partition * part0 = new mimmo::Partition();
    part0->setGeometry(mimmo0->getGeometry());
    part0->setPlotInExecution(true);
    part0->exec();

    //partition of bunny and 3D curve
    mimmo::Partition * part1 = new mimmo::Partition();
    part1->setGeometry(mimmo1->getGeometry());
    part1->setPlotInExecution(true);
    part1->exec();

    curve3D = part1->getGeometry();
    targetSurf = part0->getGeometry();

    log<<"geometries distributed on rank"<<std::endl;

#else
    curve3D = mimmo1->getGeometry();
    targetSurf = mimmo0->getGeometry();
#endif
    /*
        Project the curve onto target surface
     */
    mimmo::ProjPatchOnSurface * proj = new mimmo::ProjPatchOnSurface();
    proj->setGeometry(targetSurf);
    proj->setPatch(curve3D);
    proj->setWorkingOnTarget(false);
    proj->execute();

    //write projected curve
    proj->getProjectedElement()->getPatch()->write("utils_example_00005_projectedCurve");

    log<<"curve projected on bunny"<<std::endl;

    /* Clean up & exit;
     */
    delete mimmo0;
    delete mimmo1;
#if MIMMO_ENABLE_MPI
    delete part0;
    delete part1;
#endif
    delete proj;


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
            test00005() ;
        }
        catch(std::exception & e){
            std::cout<<"utils_example_00005 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;
}
