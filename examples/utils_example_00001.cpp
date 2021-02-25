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
	\example utils_example_00001.cpp

	\brief Example of seeding a point cloud over an unstructured non-homogeneous mesh.

	Using: MimmoGeometry, CreateSeedsOnSurface

	<b>To run</b>: ./utils_example_00001 \n
    <b>To run</b>: mpirun -np X utils_example_00001 \n

    Please note MPI version of the example partition a very coarse support mesh.
    Do not exceed with number of processors (X=2 is fine).

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00001() {

    /*
        Read a target polygonal mesh from file. CONVERT option will let the block to write
        the just read file in another file, immediately after the reading.
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    bitpit::Logger & log = mimmo0->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);

    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::SURFVTU);
    mimmo0->setReadFilename("mixedP2D");

    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::SURFVTU);
    mimmo0->setWriteFilename("utils_mesh_00001.0000");

    mimmo0->exec();
    log<<"Reading and converting the initial mesh"<<std::endl;


    mimmo::MimmoSharedPointer<mimmo::MimmoObject> target;

#if MIMMO_ENABLE_MPI
    /*
        Distribute the mesh among the processes
    */
    mimmo::Partition * part = new mimmo::Partition();
    part->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    part->setGeometry(mimmo0->getGeometry());
    part->setPlotInExecution(true);
    part->exec();

    log<<"Partitioning the initial mesh"<<std::endl;

    target = part->getGeometry();
#else
    target = mimmo0->getGeometry();
#endif

    /*
     Create a custom sensitivity field on target points
     use a linear decay in spherical coordinates from center 0.5,0.5,0.
    */
    mimmo::MimmoPiercedVector<double> sensitivity;
    sensitivity.initialize(target, mimmo::MPVLocation::POINT, 1.);

    darray3E center({{0.5,0.5,0}});
    {
        darray3E pMin,pMax;
        target->getBoundingBox(pMin, pMax, true);
        double maxD = norm2(pMax-pMin)/3.;
        double localvalue;
        for(bitpit::Vertex & v : target->getVertices()){
            localvalue =  norm2(v.getCoords() - center);
            if(localvalue < maxD){
                sensitivity[v.getId()] += -1.0*localvalue/maxD;
            }else{
                sensitivity[v.getId()] = 0.0;
            }
        }
    }
    log<<"Creating artificial sensitivity"<<std::endl;


    /*
        Creation of the seeder. Use Level set engine to place 9 points on
        surface, starting from a seed point in the absolute origin.
     */
    mimmo::CreateSeedsOnSurface * cseed = new mimmo::CreateSeedsOnSurface();
    cseed->setSeed(center);
    cseed->setNPoints(20);
    cseed->setSensitivityMap(&sensitivity);
    cseed->setEngine(0);
    cseed->setRandomFixed(true);
    cseed->setRandomSignature(1599826209);
    cseed->setPlotInExecution(true);
    cseed->setGeometry(target);

    cseed->exec();

    log<<"Seeding points on surface mesh"<<std::endl;

    /*
        plot on screen the id to reproduce random distribution if a Random engine
        is selected into cseed block.
    */
    if(cseed->getEngineENUM() == mimmo::CSeedSurf::RANDOM){
        log<<"utils_example_0001: Random seeder used signature : "<<cseed->getRandomSignature()<<std::endl;
    }

    /*
        Clean up & exit;
     */
    delete mimmo0;
    delete cseed;
#if MIMMO_ENABLE_MPI
    delete part;
#endif
    return;

}


int main(int argc, char *argv[]) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI==1
    MPI_Init(&argc, &argv);

    {
#endif

        /**<Calling core function*/
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
