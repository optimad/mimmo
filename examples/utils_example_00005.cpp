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

#include "mimmo_utils.hpp"
#if MIMMO_ENABLE_MPI
    #include "Partition.hpp"
#endif

// =================================================================================== //
/*!
	\example utils_example_00005.cpp

	\brief Project an external curve on the Stanford Bunny surface.

	Using: MimmoGeometry, ProjPatchOnSurface, Partition(MPI version)

	<b>To run</b>              : ./utils_example_00005 \n
    <b>To run(MPI version)</b> : mpirun -np X utils_example_00005 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00005() {

    /*
        Reading target surface (bunny). Convert is used to rewrite it as-it-is
        from reading in a new vtu format.
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

    /*
        Reading 3D Curve to be projected. Convert mode is used to rewrite it as-it-is once read.
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
    //distribute bunny and curve geometries among processes
    mimmo::Partition * part0 = new mimmo::Partition();
    part0->setName("PartitionedBunny");
    part0->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    part0->setGeometry(mimmo0->getGeometry());
    part0->setPlotInExecution(true);
    part0->exec();

    mimmo::Partition * part1 = new mimmo::Partition();
    part1->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    part1->setName("PartitionedCurve");
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
        Project the curve onto the bunny target surface
     */
    mimmo::ProjPatchOnSurface * proj = new mimmo::ProjPatchOnSurface();
    proj->setGeometry(targetSurf);
    proj->setPatch(curve3D);
    proj->setWorkingOnTarget(false);
    proj->execute();

    //write projected curve
    proj->getProjectedElement()->getPatch()->write("utils_example_00005_projectedCurve");

    log<<"curve projected on bunny"<<std::endl;

    //re-read and write again the resulting projected curve to test further the IO
    mimmo::MimmoGeometry * mimmo3 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);

    mimmo3->setReadDir(".");
    mimmo3->setReadFileType(FileType::CURVEVTU);
    mimmo3->setReadFilename("utils_example_00005_projectedCurve");

    mimmo3->setWriteDir(".");
    mimmo3->setWriteFileType(FileType::CURVEVTU);
    mimmo3->setWriteFilename("utils_mesh_00005_projectedCurveREREAD");
    mimmo3->execute();

    log<<"test IO read and write on the newly projected curve done"<<std::endl;

    /* Clean up & exit;
     */
    delete mimmo0;
    delete mimmo1;
    delete mimmo3;
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

        mimmo::setLogger("mimmo");

        /**<Calling core function*/
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
