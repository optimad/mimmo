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
#include "bitpit.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmo::pin;

// =================================================================================== //
/*!
	\example manipulators_example_00001.cpp

	\brief Example of usage of global deformation blocks of an input geometry.

	Global geometry deformation blocks used: translation, scaling, twisting, bending and rotation.

	<b>To run</b>: ./manipulators_example_00001 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00001() {

    /* Creation of mimmo containers.
     */
    MimmoGeometry * mimmo0 = new MimmoGeometry();

    mimmo0->setIOMode(IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("prism");

    mimmo0->setWriteDir(".");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("manipulators_output_00001.0000");

    /**<Test control on change the logger name of mimmo after a base object declaration*/
    setLogger("manipulators_example_00002");

    MimmoGeometry * mimmo1 = new MimmoGeometry();
    mimmo1->setIOMode(IOMode::WRITE);
    mimmo1->setWriteDir(".");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("manipulators_output_00001.0001");

    MimmoGeometry * mimmo2 = new MimmoGeometry();
    mimmo2->setIOMode(IOMode::WRITE);
    mimmo2->setWriteDir(".");
    mimmo2->setWriteFileType(FileType::STL);
    mimmo2->setWriteFilename("manipulators_output_00001.0002");

    MimmoGeometry * mimmo3 = new MimmoGeometry();
    mimmo3->setIOMode(IOMode::WRITE);
    mimmo3->setWriteDir(".");
    mimmo3->setWriteFileType(FileType::STL);
    mimmo3->setWriteFilename("manipulators_output_00001.0003");

    MimmoGeometry * mimmo4 = new MimmoGeometry();
    mimmo4->setIOMode(IOMode::WRITE);
    mimmo4->setWriteDir(".");
    mimmo4->setWriteFileType(FileType::STL);
    mimmo4->setWriteFilename("manipulators_output_00001.0004");

    MimmoGeometry * mimmo5 = new MimmoGeometry();
    mimmo5->setIOMode(IOMode::WRITE);
    mimmo5->setWriteDir(".");
    mimmo5->setWriteFileType(FileType::STL);
    mimmo5->setWriteFilename("manipulators_output_00001.0005");


    /* Creation of translation block.
     * Translation performed in z direction by -1.0 unit length.
     */
    TranslationGeometry* translation = new TranslationGeometry();
    translation->setDirection(darray3E{0.0, 0.0, -1.0});
    translation->setTranslation(1.0);

    /* Creation of applier block for translation.
     */
    Apply* applierTranslation = new Apply();


    /* Creation of scaling block.
     * Scaling performed in y & z direction by 0.5 factor.
     */
    ScaleGeometry* scaling = new ScaleGeometry();
    scaling->setScaling(darray3E{1.0, 0.5, 0.5});

    /* Creation of applier block for scaling.
     */
    Apply* applierScaling = new Apply();


    /* Creation of twisting block.
     * Twisting performed in x direction by a maximum angle at a distance 1.0
     * from the origin equal to pi/3 radiants.
     */
    TwistGeometry* twist = new TwistGeometry();
    twist->setDirection(darray3E{1.0, 0.0, 0.0});
    twist->setTwist((M_PI/3));
    twist->setMaxDistance(1.0);

    /* Creation of applier block for twisting.
     */
    Apply* applierTwist = new Apply();


    /* Creation of bending block.
     * Bending performed in x direction by a maximum angle at a distance 1.0
     * from the origin equal to pi/3 radiants.
     */
    BendGeometry* bend = new BendGeometry();
    umatrix33E degree = {0,0,0, 0,0,0, 2,0,0};
    bend->setDegree(degree);
    dmat33Evec coeffs;
    coeffs[2][0].resize(3);
    coeffs[2][0][0] = 0;
    coeffs[2][0][1] = 0;
    coeffs[2][0][2] = 0.5;
    bend->setCoeffs(coeffs);
    /*
     * Bend directly applied during execution
     */
    bend->setApply();

    /* Creation of rotation block.
     * Rotation performed around an axis through the origin and
     * with direction (0.25,0.25,0.75) by pi/4 radiants.
     */
    RotationGeometry* rotation = new RotationGeometry();
    rotation->setDirection(darray3E{0.25,0.25,0.75});
    rotation->setRotation((M_PI/4));

    /* Creation of applier block for rotation.
     */
    Apply* applierRotation = new Apply();


    /* Setup pin connections.
     */
    addPin(mimmo0, translation, M_GEOM, M_GEOM);
    addPin(mimmo0, applierTranslation, M_GEOM, M_GEOM);
    addPin(translation, applierTranslation, M_GDISPLS, M_GDISPLS);
    addPin(applierTranslation, mimmo1, M_GEOM, M_GEOM);

    addPin(applierTranslation, scaling, M_GEOM, M_GEOM);
    addPin(applierTranslation, applierScaling, M_GEOM, M_GEOM);
    addPin(scaling, applierScaling, M_GDISPLS, M_GDISPLS);
    addPin(applierScaling, mimmo2, M_GEOM, M_GEOM);

    addPin(applierScaling, twist, M_GEOM, M_GEOM);
    addPin(applierScaling, applierTwist, M_GEOM, M_GEOM);
    addPin(twist, applierTwist, M_GDISPLS, M_GDISPLS);
    addPin(applierTwist, mimmo3, M_GEOM, M_GEOM);

    addPin(applierTwist, bend, M_GEOM, M_GEOM);
    addPin(bend, mimmo4, M_GEOM, M_GEOM);

    addPin(bend, rotation, M_GEOM, M_GEOM);
    addPin(bend, applierRotation, M_GEOM, M_GEOM);
    addPin(rotation, applierRotation, M_GDISPLS, M_GDISPLS);
    addPin(applierRotation, mimmo5, M_GEOM, M_GEOM);


    /* Setup execution chain.
     */
    Chain ch0;
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
    ch0.exec();

    /* Clean up & exit;
     */
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

    mimmo0 = NULL;
    translation = NULL;
    applierTranslation = NULL;
    mimmo1 = NULL;
    scaling = NULL;
    applierScaling = NULL;
    mimmo2 = NULL;
    twist = NULL;
    applierTwist = NULL;
    mimmo3 = NULL;
    bend = NULL;
    mimmo4 = NULL;
    rotation = NULL;
    applierRotation = NULL;
    mimmo5 = NULL;

    return;

}


int main(int argc, char *argv[]) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc, &argv);

    {
#endif

        /**<Change the name of mimmo logger file (default mimmo.log)
         * before initialization of BaseManipulation objects*/
        mimmo::setLogger("mimmo");

        /**<Calling mimmo Test routines*/
        test00001() ;

#if BITPIT_ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return (1);
}


