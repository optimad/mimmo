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
 * Using: IOCGNS, MeshChecker.
 *
 * Depends on mimmo optional module geohandlers
 *
 * <b>To run</b>: ./iocgns_example_00002  \n
 *
 * <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //

void example00002() {

    /* Create IO_CGNS object to import input file. */
	mimmo::IOCGNS * cgnsI = new mimmo::IOCGNS();
    cgnsI->setMode(mimmo::IOCGNS::IOCGNS_Mode::READ);
    cgnsI->setDir("geodata");
    cgnsI->setFilename("grid");

#if MIMMO_ENABLE_MPI
    /* Instantiation of a Partition object with default patition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition *partition = new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition->setPlotInExecution(true);
#endif

    mimmo::MeshChecker* checkmesh = new mimmo::MeshChecker();
    checkmesh->setPlotInExecution(true);
    checkmesh->setMinimumVolumeTolerance(1.0E-1);
    checkmesh->setMaximumVolumeTolerance(1.0E9);
    checkmesh->setMaximumSkewnessTolerance(80.0);
    checkmesh->setMaximumBoundarySkewnessTolerance(90.0);
    checkmesh->setMinimumFaceValidityTolerance(0.2);
    checkmesh->setMinimumVolumeChangeTolerance(5.0E-2);

    /* Create PINs. */
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
            /**< Call mimmo example routine. */
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
