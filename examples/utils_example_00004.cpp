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
#include "FFDLattice.hpp"
#if MIMMO_ENABLE_MPI
    #include "Partition.hpp"
#endif

// =================================================================================== //
/*!
	\example utils_example_00004.cpp

	\brief Comparison of Oriented Bounding Boxes between a deformed and undeformed geometry.

	Using: MimmoGeometry, OBBox, FFDLattice, TranslationPoint, RotationAxes

	<b>To run</b>: ./utils_example_00004 \n
    <b>To run</b>: mpirun -np x utils_example_00004 \n
	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00004() {

    /* Creation of mimmo containers for target StanfordBunny.
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setName("StanfordBunnyReader");
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("stanfordBunny2");
    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::SURFVTU);
    mimmo0->setWriteFilename("utils_mesh_00004.0000");

#if MIMMO_ENABLE_MPI
    /* Partitionate  StanfordBunny.
     */
	mimmo::Partition * mimmo0part = new mimmo::Partition();
    mimmo0part->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    mimmo0part->setName("StanfordBunnyPartitioner");
    mimmo0part->setPlotInExecution(true);
#endif

    /*
        Calculate the OBB the original StanfordBunny
     */
    mimmo::OBBox * obb_original = new mimmo::OBBox();
    obb_original->setName("OBBOriginal");
    obb_original->setWriteInfo(true);
    obb_original->setPlotInExecution(true);

    /*
        Calculate the OBB the deformed StanfordBunny
     */
    mimmo::OBBox * obb_deformed = new mimmo::OBBox();
    obb_deformed->setName("OBBDeformed");
    obb_deformed->setWriteInfo(true);
    obb_deformed->setPlotInExecution(true);

    /* Translate point. This is meant as the new origin of the FFDLattice for deformation
     */
    mimmo::TranslationPoint * translp = new mimmo::TranslationPoint();
    translp->setName("TranslationOriginLattice");
    translp->setOrigin({{0.0,0.0,0.0}});
    translp->setDirection({{-0.714,0.0,1.0}});
    translp->setTranslation(0.9);

    /* Rotate axes reference system. This is meant as the new sdr axes of the FFDLattice for deformation
     */
    mimmo::RotationAxes * rot_axes = new mimmo::RotationAxes();
    rot_axes->setName("SDRAxesLattice");
    rot_axes->setOrigin({{0.0,0.0,0.0}});
    rot_axes->setDirection({{0.0,0.0,1.0}});
    rot_axes->setRotation(-30.0*BITPIT_PI/180.0);
    rot_axes->setAxes({{1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0}});
    rot_axes->setAxesOrigin({{0.0,0.0,0.0}});

    /*
        Create a cubic FFD Lattice to deform the head of the rabbit.
        Axes and Origin of the cube are provided by rot_axes and translp
        through ports.
    */
    mimmo::FFDLattice * latt = new mimmo::FFDLattice();
    latt->setName("FFDLattice");
    latt->setShape(mimmo::ShapeType::CUBE);
    latt->setSpan({{0.7,1.2,0.62}});
    latt->setDimension(iarray3E({{2,2,2}}));
    latt->setDegrees(iarray3E({{1,1,1}}));
    dvecarr3E displs(8,{{0.0,0.0,0.0}});
    displs[1][0] = -0.5;
    displs[3][0] = -0.5;
    latt->setDisplacements(displs);
    latt->setApply(true);
    latt->setPlotInExecution(true);

    // create connections
    // original geoemetry passed to obb and lattice
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, mimmo0part, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0part, obb_original, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0part, latt, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(mimmo0, obb_original, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, latt, M_GEOM, M_GEOM);
#endif

    // origin, axes passed to lattice
    mimmo::pin::addPin(translp, latt, M_POINT, M_POINT);
    mimmo::pin::addPin(rot_axes, latt, M_AXES, M_AXES);

    // pass deformed geoemetry to obb_deformed
    mimmo::pin::addPin(latt, obb_deformed, M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    mimmo::Chain ch0,ch1;
    ch0.addObject(mimmo0);
    ch0.addObject(obb_original);
#if MIMMO_ENABLE_MPI
    ch0.addObject(mimmo0part);
#endif

    ch1.addObject(translp);
    ch1.addObject(rot_axes);
    ch1.addObject(latt);
    ch1.addObject(obb_deformed);

    /* Execution of chain.
     * Use debug flag false to avoid to print out the execution steps on console.
     */
    ch0.exec(true);
    ch1.exec(true);


    /*
        Write deformed geometry;
    */
    mimmo0->getGeometry()->getPatch()->write("utils_mesh_00004.0001");

    /* Clean up & exit;
     */
    delete mimmo0;
    delete translp;
    delete rot_axes;
    delete obb_original;
    delete obb_deformed;
    delete latt;
#if MIMMO_ENABLE_MPI
    delete mimmo0part;
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

        /**<Calling mimmo Test routines*/
        try{
            test00004() ;
        }
        catch(std::exception & e){
            std::cout<<"utils_example_00004 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;
}
