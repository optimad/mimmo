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
#include <mimmo_iocgns.hpp>
#include <mimmo_utils.hpp>
#if MIMMO_ENABLE_MPI
    #include <mimmo_parallel.hpp>
#endif

/*!
 * \example iocgns_example_00002.cpp
 *
 * \brief Reading of a CGNS volume mesh and check it with MeshChecker class.
 *
 * Using: IOCGNS, MeshChecker, Partition(MPI version).
 *
 * Depends on mimmo optional module geohandlers
 *
 * <b>To run</b>             : ./iocgns_example_00002  \n
 * <b>To run(MPI version)</b>: mpirun -np X iocgns_example_00002  \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //

void example00002() {

    /*
       Create IOCGNS object to import the cgns mesh.
       Bulk volume and relative boundary mesh will be exposed
    */
	mimmo::IOCGNS * cgnsI = new mimmo::IOCGNS(mimmo::IOCGNS::IOCGNS_Mode::READ);
    cgnsI->setDir("geodata");
    cgnsI->setFilename("grid");

#if MIMMO_ENABLE_MPI
    /*
        Distribute bulk/boundary meshes among processes
     */
    mimmo::Partition *partition = new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition->setPlotInExecution(true);
#endif

    /*
        Check the status of the mesh elements: min/max volume, skewness
        skewness of the boundary cell elements etc...
    */
    mimmo::MeshChecker* checkmesh = new mimmo::MeshChecker();
    checkmesh->setPlotInExecution(true);
    checkmesh->setMinimumVolumeTolerance(1.0E-1);
    checkmesh->setMaximumVolumeTolerance(1.0E9);
    checkmesh->setMaximumSkewnessTolerance(80.0);
    checkmesh->setMaximumBoundarySkewnessTolerance(90.0);
    checkmesh->setMinimumFaceValidityTolerance(0.2);
    checkmesh->setMinimumVolumeChangeTolerance(5.0E-2);

    /* Create block connections. */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(cgnsI, partition, M_GEOM, M_GEOM)  ;
    mimmo::pin::addPin(cgnsI, partition, M_GEOM2, M_GEOM2)  ;
    mimmo::pin::addPin(partition, checkmesh, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(cgnsI, checkmesh, M_GEOM, M_GEOM)  ;
#endif

    /* Create and execute chain. */
    mimmo::Chain ch0;
    ch0.addObject(cgnsI);
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition);
#endif
    ch0.addObject(checkmesh);

    ch0.exec(true);

    /* Destroy objects. */
    delete cgnsI;
#if MIMMO_ENABLE_MPI
    delete partition;
#endif
    delete checkmesh;

    return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);

    {
#endif
        try{
            /**< Call core function. */
            example00002();
        }
        catch(std::exception & e){
            std::cout<<"iocgns_example_00002 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
    }

    MPI_Finalize();
#endif

    return 0;
}
