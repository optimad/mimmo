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

#include "mimmo_geohandlers.hpp"
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

// =================================================================================== //
/*!
	\example geohandlers_example_00004.cpp

	\brief Using Pid Selection, Clipping and  Surface Triangulator to extract a mesh

	Geometry handler block used: SelectionByPID, CLipGeometry, SurfaceTriangulator.

	<b>To run</b>: ./geohandlers_example_00004 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //

/*
 * Creating surface polygonal mesh and return it in a MimmoObject.
 * Pidding the "M I O" texture with PID=1;
 *
 * \return unique ptr to the polygonal surface mesh
 */
mimmo::MimmoSharedPointer<mimmo::MimmoObject> createMIOMesh(){

    mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh(new mimmo::MimmoObject(1));

    int rank = -1;
#if MIMMO_ENABLE_MPI
    rank = mesh->getRank();
    if ( rank == 0)
#endif
    {
        //create the vertices set.
        mesh->addVertex({{0.0,0.0,0.0}}, 0);
        mesh->addVertex({{0.0,1.0,0.0}}, 1);
        mesh->addVertex({{0.3,0.0,0.0}}, 2);
        mesh->addVertex({{0.3,1.0,0.0}}, 3);
        mesh->addVertex({{0.4,0.0,0.0}}, 4);
        mesh->addVertex({{0.4,0.75,0.0}}, 5);
        mesh->addVertex({{0.4,1.0,0.0}}, 6);
        mesh->addVertex({{0.65,0.25,0.0}}, 7);
        mesh->addVertex({{0.65,0.5,0.0}}, 8);
        mesh->addVertex({{0.9,0.0,0.0}}, 9);
        mesh->addVertex({{0.9,0.75,0.0}}, 10);
        mesh->addVertex({{0.9,1.0,0.0}}, 11);
        mesh->addVertex({{1.0,0.0,0.0}}, 12);
        mesh->addVertex({{1.0,1.0,0.0}}, 13);
        mesh->addVertex({{1.3,0.0,0.0}}, 14);
        mesh->addVertex({{1.3,1.0,0.0}}, 15);
        mesh->addVertex({{1.4,0.0,0.0}}, 16);
        mesh->addVertex({{1.4,1.0,0.0}}, 17);
        mesh->addVertex({{1.7,0.0,0.0}}, 20);
        mesh->addVertex({{1.7,1.0,0.0}}, 21);
        mesh->addVertex({{1.8,0.0,0.0}}, 22);
        mesh->addVertex({{1.8,0.1,0.0}}, 23);
        mesh->addVertex({{1.8,0.9,0.0}}, 24);
        mesh->addVertex({{1.8,1.0,0.0}}, 25);
        mesh->addVertex({{2.3,0.0,0.0}}, 26);
        mesh->addVertex({{2.3,0.1,0.0}}, 27);
        mesh->addVertex({{2.3,0.9,0.0}}, 28);
        mesh->addVertex({{2.3,1.0,0.0}}, 29);
        mesh->addVertex({{2.4,0.0,0.0}}, 30);
        mesh->addVertex({{2.4,1.0,0.0}}, 31);
        mesh->addVertex({{2.7,0.0,0.0}}, 32);
        mesh->addVertex({{2.7,1.0,0.0}}, 33);

        //add and create polygonal cells
        mesh->addConnectedCell(livector1D({{0,2,3,1}}), bitpit::ElementType::QUAD, long(0), long(0), rank);
        mesh->addConnectedCell(livector1D({{5,2,4,5,6,3}}), bitpit::ElementType::POLYGON, long(1), long(1), rank);
        mesh->addConnectedCell(livector1D({{5,4,9,10,7,5}}), bitpit::ElementType::POLYGON, long(0), long(2), rank);
        mesh->addConnectedCell(livector1D({{5,7,8,6}}), bitpit::ElementType::QUAD, long(1), long(3), rank);
        mesh->addConnectedCell(livector1D({{7,10,11,8}}), bitpit::ElementType::QUAD, long(1), long(4), rank);
        mesh->addConnectedCell(livector1D({{6,8,11}}), bitpit::ElementType::TRIANGLE, long(0), long(5), rank);
        mesh->addConnectedCell(livector1D({{5,9,12,13,11,10}}), bitpit::ElementType::POLYGON, long(1), long(6), rank);
        mesh->addConnectedCell(livector1D({{12,14,15,13}}), bitpit::ElementType::QUAD, long(0), long(7), rank);
        mesh->addConnectedCell(livector1D({{14,16,17,15}}), bitpit::ElementType::QUAD, long(1), long(8), rank);
        mesh->addConnectedCell(livector1D({{16,20,21,17}}), bitpit::ElementType::QUAD, long(0), long(9), rank);
        mesh->addConnectedCell(livector1D({{6, 20,22, 23,24,25,21}}), bitpit::ElementType::POLYGON, long(1), long(10), rank);
        mesh->addConnectedCell(livector1D({{22,26,27,23}}), bitpit::ElementType::QUAD, long(1), long(11), rank);
        mesh->addConnectedCell(livector1D({{23,27,28,24}}), bitpit::ElementType::QUAD, long(0), long(12), rank);
        mesh->addConnectedCell(livector1D({{24,28,29,25}}), bitpit::ElementType::QUAD, long(1), long(13), rank);
        mesh->addConnectedCell(livector1D({{6,26,30,31,29,28,27}}), bitpit::ElementType::POLYGON, long(1), long(14), rank);
        mesh->addConnectedCell(livector1D({{30,32,33,31}}), bitpit::ElementType::QUAD, long(0), long(15), rank);

    }

	return mesh;
}



void test00001() {

    mimmo::setExpertMode(true);

    /*
     * Create a surface polygonal mesh with texture MIO pidded as PID=1
     */
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> geo = createMIOMesh();

#if MIMMO_ENABLE_MPI
    /* Instantiation of a Partition object with default partition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    std::unique_ptr<mimmo::Partition> partition(new mimmo::Partition());
    partition->setPlotInExecution(true);
    partition->setGeometry(geo);
#endif

    /*
     * extract texture MIO with SelectionByPID
     */
    std::unique_ptr<mimmo::SelectionByPID> sel(new mimmo::SelectionByPID());
    sel->setName("PIDExtraction");
#if MIMMO_ENABLE_MPI==0
    sel->setGeometry(geo);
#endif
    sel->setPID(1);
    sel->setDual(false);
    sel->setPlotInExecution(true);

    /*
     * isolate M from MIO using a plane clipping
     */
    std::unique_ptr<mimmo::ClipGeometry> clip(new mimmo::ClipGeometry());
    clip->setName("PlaneClipping");
    clip->setOrigin({{1.1,0.0,0.0}});
    clip->setNormal({{1.0,0.0,0.0}});
    clip->setInsideOut(true);
    clip->setPlotInExecution(true);

    /*
     * triangulate the M polygonal tessellation
     */
    std::unique_ptr<mimmo::SurfaceTriangulator> triang(new mimmo::SurfaceTriangulator());
    triang->setName("TriangulateSurface");
    triang->setPlotInExecution(true);


    /* Setup pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(partition.get(), sel.get(), M_GEOM, M_GEOM);
#endif
    mimmo::pin::addPin(sel.get(), clip.get(), M_GEOM, M_GEOM);
    mimmo::pin::addPin(clip.get(), triang.get(), M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition.get());
#endif
    ch0.addObject(sel.get());
    ch0.addObject(clip.get());
    ch0.addObject(triang.get());

    /* Execution of chain.
     * Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

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
            /**<Calling mimmo Test routine*/
            test00001();
        }
        catch(std::exception & e){
            std::cout<<"geohandlers_example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
