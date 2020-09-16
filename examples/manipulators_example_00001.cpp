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

#include "mimmo_manipulators.hpp"
#include "mimmo_iogeneric.hpp"
#include <bitpit_common.hpp>
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

// =================================================================================== //
/*!
	\example manipulators_example_00001.cpp

	\brief Example of usage of global deformation blocks of an input geometry.

	Global geometry deformation blocks used: TranslationGeometry, ScaleGeometry, TwistGeometry, BendGeometry and RotationGeometry.

	<b>To run</b>: ./manipulators_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00001() {

    /* Creation of mimmo containers.
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("prism");

    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("manipulators_output_00001.0000");

    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("manipulators_output_00001.0001");

    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo2->setWriteDir(".");
    mimmo2->setWriteFileType(FileType::STL);
    mimmo2->setWriteFilename("manipulators_output_00001.0002");

    mimmo::MimmoGeometry * mimmo3 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo3->setWriteDir(".");
    mimmo3->setWriteFileType(FileType::STL);
    mimmo3->setWriteFilename("manipulators_output_00001.0003");

    mimmo::MimmoGeometry * mimmo4 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo4->setWriteDir(".");
    mimmo4->setWriteFileType(FileType::STL);
    mimmo4->setWriteFilename("manipulators_output_00001.0004");

    mimmo::MimmoGeometry * mimmo5 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo5->setWriteDir(".");
    mimmo5->setWriteFileType(FileType::STL);
    mimmo5->setWriteFilename("manipulators_output_00001.0005");

#if MIMMO_ENABLE_MPI
    /* Creation of a Partition object.
     */
    mimmo::Partition* partition= new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
#endif

    /* Creation of translation block.
     * Translation performed in z direction by -1.0 unit length.
     */
    mimmo::TranslationGeometry* translation = new mimmo::TranslationGeometry();
    translation->setDirection(darray3E{0.0, 0.0, -1.0});
    translation->setTranslation(1.0);

    /* Creation of applier block for translation.
     */
    mimmo::Apply* applierTranslation = new mimmo::Apply();


    /* Creation of scaling block.
     * Scaling performed in y & z direction by 0.5 factor.
     */
    mimmo::ScaleGeometry* scaling = new mimmo::ScaleGeometry();
    scaling->setScaling(darray3E{1.0, 0.5, 0.5});

    /* Creation of applier block for scaling.
     */
    mimmo::Apply* applierScaling = new mimmo::Apply();

    /* Creation of twisting block.
     * Twisting performed in x direction by a maximum angle at a distance 1.0
     * from the origin equal to pi/3 radiants.
     */
    mimmo::TwistGeometry* twist = new mimmo::TwistGeometry();
    twist->setDirection(darray3E{1.0, 0.0, 0.0});
    twist->setTwist((BITPIT_PI/3));
    twist->setMaxDistance(1.0);

    /* Creation of applier block for twisting.
     */
    mimmo::Apply* applierTwist = new mimmo::Apply();

    /* Creation of bending block.
     * Bending performed in x direction by a maximum angle at a distance 1.0
     * from the origin equal to pi/3 radiants.
     */
    mimmo::BendGeometry* bend = new mimmo::BendGeometry();
    umatrix33E degree = {0,0,0, 0,0,0, 2,0,0};
    bend->setDegree(degree);
    dmat33Evec coeffs;
    coeffs[2][0].resize(3);
    coeffs[2][0][0] = 0;
    coeffs[2][0][1] = 0;
    coeffs[2][0][2] = 0.5;
    bend->setCoeffs(&coeffs);

    /*
     * Bend directly applied during execution
     */
    bend->setApply();

    /* Creation of rotation block.
     * Rotation performed around an axis through the origin and
     * with direction (0.25,0.25,0.75) by pi/4 radiants.
     */
    mimmo::RotationGeometry* rotation = new mimmo::RotationGeometry();
    rotation->setDirection(darray3E{0.25,0.25,0.75});
    rotation->setRotation((BITPIT_PI/4));

    /* Creation of applier block for rotation.
     */
    mimmo::Apply* applierRotation = new mimmo::Apply();

    /* Setup pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, partition, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition, translation, M_GEOM, M_GEOM);
#else
    mimmo::pin::addPin(mimmo0, translation, M_GEOM, M_GEOM);
#endif
    mimmo::pin::addPin(mimmo0, applierTranslation, M_GEOM, M_GEOM);
    mimmo::pin::addPin(translation, applierTranslation, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applierTranslation, mimmo1, M_GEOM, M_GEOM);

    mimmo::pin::addPin(applierTranslation, scaling, M_GEOM, M_GEOM);
    mimmo::pin::addPin(applierTranslation, applierScaling, M_GEOM, M_GEOM);
    mimmo::pin::addPin(scaling, applierScaling, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applierScaling, mimmo2, M_GEOM, M_GEOM);

    mimmo::pin::addPin(applierScaling, twist, M_GEOM, M_GEOM);
    mimmo::pin::addPin(applierScaling, applierTwist, M_GEOM, M_GEOM);
    mimmo::pin::addPin(twist, applierTwist, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applierTwist, mimmo3, M_GEOM, M_GEOM);

    mimmo::pin::addPin(applierTwist, bend, M_GEOM, M_GEOM);
    mimmo::pin::addPin(bend, mimmo4, M_GEOM, M_GEOM);

    mimmo::pin::addPin(bend, rotation, M_GEOM, M_GEOM);
    mimmo::pin::addPin(bend, applierRotation, M_GEOM, M_GEOM);
    mimmo::pin::addPin(rotation, applierRotation, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(applierRotation, mimmo5, M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition);
#endif
    ch0.addObject(mimmo0);
    ch0.addObject(translation);
    ch0.addObject(applierTranslation);
    ch0.addObject(mimmo1);
    ch0.addObject(scaling);
    ch0.addObject(applierScaling);
    ch0.addObject(mimmo2);
    ch0.addObject(twist);
    ch0.addObject(applierTwist);
    ch0.addObject(mimmo3);
    ch0.addObject(bend);
    ch0.addObject(mimmo4);
    ch0.addObject(rotation);
    ch0.addObject(applierRotation);
    ch0.addObject(mimmo5);

    /* Execution of chain.
     * Use debug flag false to avoid to print out the execution steps on console.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
#if MIMMO_ENABLE_MPI
    delete partition;
#endif
    delete mimmo0;
    delete translation;
    delete applierTranslation;
    delete mimmo1;
    delete scaling;
    delete applierScaling;
    delete mimmo2;
    delete twist;
    delete applierTwist;
    delete mimmo3;
    delete bend;
    delete mimmo4;
    delete rotation;
    delete applierRotation;
    delete mimmo5;

#if MIMMO_ENABLE_MPI
    partition = nullptr;
#endif
    mimmo0 = nullptr;
    translation = nullptr;
    applierTranslation = nullptr;
    mimmo1 = nullptr;
    scaling = nullptr;
    applierScaling = nullptr;
    mimmo2 = nullptr;
    twist = nullptr;
    applierTwist = nullptr;
    mimmo3 = nullptr;
    bend = nullptr;
    mimmo4 = nullptr;
    rotation = nullptr;
    applierRotation = nullptr;
    mimmo5 = nullptr;

    return;

}

int main(int argc, char *argv[]) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI==1
    MPI_Init(&argc, &argv);

    {
#endif


        try{
            /**<Calling mimmo Test routines*/
            test00001() ;
        }

        catch(std::exception & e){
            std::cout<<"manipulators_example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;
}
