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
#include "MeshSelection.hpp"
#include "ReconstructFields.hpp"
#include "MRBF.hpp"
#include "Apply.hpp"
#if MIMMO_ENABLE_MPI
    #include "Partition.hpp"
#endif
// =================================================================================== //
/*!
	\example utils_example_00003.cpp

	\brief Deform a sphere using a set of mirrored RBF points. Check for collisions
           of the deformed object with the D=0.25 level-set isolevel of the original geometry.

	Using: MimmoGeometry, SpecularPoints, RBFBox, ControlDeformMaxDistance,
           MRBF, Applier, SelectionByBox, ReconstructVector

	<b>To run</b>: ./utils_example_00003 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00003() {

    mimmo::setExpertMode(true); // avoid control of unlinked ports.

    /* Creation of mimmo containers.
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setName("converterGeometry");
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::SURFVTU);
    mimmo0->setWriteFilename("utils_mesh_00003.0000");

#if MIMMO_ENABLE_MPI
    /* Instantiation of a Partition object for target mesh partitioning.
    */
    mimmo::Partition* partitioner = new mimmo::Partition();
    partitioner->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partitioner->setName("mimmo.Partitioner");
    partitioner->setPlotInExecution(true);
#endif

    //declare a support Radius for RBF cloud
    double suppR = 0.2;

    //create a RBF cloud with a set of displacements.
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> rbfPointCloud(new mimmo::MimmoObject(3));
    mimmo::MimmoPiercedVector<darray3E> rbfDispls(rbfPointCloud, mimmo::MPVLocation::POINT);

#if MIMMO_ENABLE_MPI
    if(rbfPointCloud->getRank() == 0)
    {
#endif

    //create 4 points on a y-plane at {{0.0,-0.6,0.0}}
    rbfPointCloud->addVertex({{0.0,-0.6,0.0}}, 0);
    rbfPointCloud->addVertex({{-0.1,-0.6,0.0}}, 1);
    rbfPointCloud->addVertex({{0.1,-0.6,0.12}}, 2);
    rbfPointCloud->addVertex({{-0.03,-0.6,-0.1}}, 3);

    rbfPointCloud->addConnectedCell(std::vector<long>(1,0), bitpit::ElementType::VERTEX, 0, 10, 0);
    rbfPointCloud->addConnectedCell(std::vector<long>(1,1), bitpit::ElementType::VERTEX, 0, 11, 0);
    rbfPointCloud->addConnectedCell(std::vector<long>(1,2), bitpit::ElementType::VERTEX, 0, 12, 0);
    rbfPointCloud->addConnectedCell(std::vector<long>(1,3), bitpit::ElementType::VERTEX, 0, 13, 0);

    // insert displacements for POINTS
    rbfDispls.insert(0, {{0.0,-0.2,0.0}});
    rbfDispls.insert(1, {{-0.6,-0.2,0.0}});
    rbfDispls.insert(2, {{0.16,-0.2,0.4}});
    rbfDispls.insert(3, {{-0.4,-0.2,-0.18}});

#if MIMMO_ENABLE_MPI
    }
    rbfPointCloud->cleanPatchInfo();
    rbfPointCloud->updateAdjacencies();
    rbfPointCloud->update();
#endif

    /*
        Project onto surface and mirror these rbf points with their displacements attached
        w.r.t. the y plane crossing the origin
    */
    mimmo::SpecularPoints * spec = new mimmo::SpecularPoints();
    spec->setName("MirroringRBF");
    spec->setPointCloud(rbfPointCloud);
    spec->setVectorData(&rbfDispls);
    spec->setOrigin({{0.0,0.0,0.0}});
    spec->setNormal({{0.0,1.0,0.0}});
    spec->setInsideOut(true);
    spec->setForce(false);
    spec->setPlotInExecution(true);


    /*
        Calculate AABB of the Mirrored RBF set accounting for the current support
        radius chosen. Mirrored points are passed through port connection.
    */
    mimmo::RBFBox * rbfbox = new mimmo::RBFBox();
    rbfbox->setName("AABB_RBFBox");
    rbfbox->setSupportRadius(suppR);
    rbfbox->setPlotInExecution(true);
    /*
        Use RBFBox to provide a sub-selection of the original geometry.
        Origin and span of RBFBox AABB are provided through ports.
     */
    mimmo::SelectionByBox * sel = new mimmo::SelectionByBox();
    sel->setName("RBFBoxSelection");
    sel->setDual(false);
    sel->setPlotInExecution(true);

    /*
        Get deformation on the box selected surface with the current RBF mirrored set,
        but not apply it yet.
    */
    mimmo::MRBF * manip = new mimmo::MRBF();
    manip->setName("RBFManipulator");
    manip->setSupportRadiusReal(suppR);

    /*
        Reconstruct the deformation field on the whole body
    */
    mimmo::ReconstructVector * recon = new mimmo::ReconstructVector();
    recon->setName("ReconstructDeformation");
    recon->setPlotInExecution(true);

    /*
        Verify actual geometry deformation vs the original geometry
        level-set isolevel at d=0.14. Get the max penetration value on file.
        Original geometry and its to-check deformation are passed through ports
    */
    mimmo::ControlDeformMaxDistance * isolevelCheck = new mimmo::ControlDeformMaxDistance();
    isolevelCheck->setName("IsoLevelCheckCollision");
    isolevelCheck->setLimitDistance(0.5);
    isolevelCheck->setPlotInExecution(true);

    /*
        Apply deformation to the original geometry
    */
    mimmo::Apply * applier = new mimmo::Apply();

#if MIMMO_ENABLE_MPI
    /* Instantiation of a Partition object for target mesh reserialization.
    */
    mimmo::Partition* serialize = new mimmo::Partition();
    serialize->setName("mimmo.Serialization");
    serialize->setPlotInExecution(false);
    serialize->setPartitionMethod(mimmo::PartitionMethod::SERIALIZE);
#endif




    /* Setup pin connections.
     */
    // original geometry passed to spec
#if MIMMO_ENABLE_MPI
    //read target geometry, partition it  and distribute it to the other subblocks
    mimmo::pin::addPin(mimmo0, partitioner, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partitioner, spec, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partitioner, sel, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partitioner, recon, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partitioner, isolevelCheck, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partitioner, applier, M_GEOM, M_GEOM);

#else
    //read target geometry and distribute it to the other subblocks
    mimmo::pin::addPin(mimmo0, spec, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, sel, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, recon, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, isolevelCheck, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);

#endif
    //passing mirrored pc to rbfbox
    mimmo::pin::addPin(spec, rbfbox, M_GEOM, M_GEOM);

    //passing box dimension from rbfbox to box selection
    mimmo::pin::addPin(rbfbox, sel, M_POINT, M_POINT);
    mimmo::pin::addPin(rbfbox, sel, M_SPAN, M_SPAN);

    // manip get as working geometry the selected patch of sel.
    // rbf point cloud and displacements are provided by spec
    mimmo::pin::addPin(sel, manip, M_GEOM, M_GEOM);
    mimmo::pin::addPin(spec, manip, M_GEOM, M_GEOM2);
    mimmo::pin::addPin(spec, manip, M_VECTORFIELD, M_VECTORFIELD);

    //pass the deformation field of manip to reconstructor
    mimmo::pin::addPin(manip, recon, M_GDISPLS, M_VECTORFIELD);

    //pass the reconstructed deformation field to isolevelCheck
    mimmo::pin::addPin(recon, isolevelCheck, M_VECTORFIELD, M_GDISPLS);

    // pass the reconstructed deformation to applier
    mimmo::pin::addPin(recon, applier, M_VECTORFIELD, M_GDISPLS);

#if MIMMO_ENABLE_MPI
    //add serialization
    mimmo::pin::addPin(applier, serialize, M_GEOM, M_GEOM);
#endif
    /* Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(manip);
    ch0.addObject(spec);
    ch0.addObject(rbfbox);
    ch0.addObject(sel);
    ch0.addObject(recon);
    ch0.addObject(isolevelCheck);
    ch0.addObject(applier);
#if MIMMO_ENABLE_MPI
    ch0.addObject(partitioner);
    ch0.addObject(serialize);
#endif

    /* Execution of chain.
     * Use debug flag false to avoid to print out the execution steps on console.
     */
    ch0.exec(true);

   //Last step write deformed original geometry in vtu
#if MIMMO_ENABLE_MPI
    mimmo0->getGeometry()->getPatch()->write("utils_mesh_00003.0001");
#else
    serialize->getGeometry()->getPatch()->write("utils_mesh_00003.0001");
#endif

    /* Clean up & exit;
     */
    delete mimmo0;
    delete spec;
    delete rbfbox;
    delete sel;
    delete manip;
    delete recon;
    delete isolevelCheck;
    delete applier;
#if MIMMO_ENABLE_MPI
    delete partitioner;
    delete serialize;
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

        /**<Change the name of mimmo logger file (default mimmo.log)
         * before initialization of BaseManipulation objects*/
        mimmo::setLogger("mimmo");

        /**<Calling mimmo Test routines*/
        try{
            test00003() ;
        }
        catch(std::exception & e){
            std::cout<<"utils_example_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;
}
