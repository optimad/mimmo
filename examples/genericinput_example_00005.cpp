/*----------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
 *  Copyright (c) 2015-2020 Optimad Engineering S.r.l., All Rights Reserved.
 *
 *  --------------------------------------------------------------------------
 *
 *  NOTICE:  All information contained herein is, and remains the property
 *  of COMPANY. The intellectual and technical concepts contained herein are
 *  proprietary to COMPANY and may be covered by Italian and Foreign Patents,
 *  patents in process, and are protected by trade secret or copyright law.
 *  Dissemination of this information or reproduction of this material is
 *  strictly forbidden unless prior written permission is obtained from
 *  COMPANY. Access to the source code contained herein is hereby forbidden
 *  to anyone except current COMPANY employees, managers or contractors who
 *  have executed Confidentiality and Non-disclosure agreements explicitly
 *  covering such access.
 *
 *  The copyright notice above does not evidence any actual or intended
 *  publication or disclosure of this source code, which includes information
 *  that is confidential and/or proprietary, and is a trade secret, of
 *  COMPANY. ANY REPRODUCTION, MODIFICATION, DISTRIBUTION, PUBLIC PERFORMANCE,
 *  OR PUBLIC DISPLAY OF OR THROUGH USE  OF THIS  SOURCE CODE  WITHOUT THE
 *  EXPRESS WRITTEN CONSENT OF COMPANY IS STRICTLY PROHIBITED, AND IN
 *  VIOLATION OF APPLICABLE LAWS AND INTERNATIONAL TREATIES. THE RECEIPT OR
 *  POSSESSION OF THIS SOURCE CODE AND/OR RELATED INFORMATION DOES NOT CONVEY
 *  OR IMPLY ANY RIGHTS TO REPRODUCE, DISCLOSE OR DISTRIBUTE ITS CONTENTS, OR
 *  TO MANUFACTURE, USE, OR SELL ANYTHING THAT IT  MAY DESCRIBE, IN WHOLE OR
 *  IN PART.
 *
\*----------------------------------------------------------------------------*/

#include <mimmo_iogeneric.hpp>
#include <mimmo_geohandlers.hpp>
#include <mimmo_manipulators.hpp>
#include <mimmo_utils.hpp>


/*!
    \example genericinput_example_00005.cpp

    \brief Example of usage of IOWavefrontOBJ and MRBF Style Object block

    mimic core blocks used: Proj3DCurveOnSurface, MRBFStyleObj
    mimmo main blocks used: IOWavefrontOBJ, Apply

    <b>To run</b>: ./genericinput_example_00005 \n

 */

int genericinput00005() {

    /*
     * Set expert mode true
     */
    mimmo::setExpertMode();

    /* Creation of mimmo containers.
     * Input and output MimmoGeometry are instantiated
     * as two different objects (no loop in chain are permitted).
     */
    mimmo::IOWavefrontOBJ * mimmo0 = new mimmo::IOWavefrontOBJ(mimmo::IOWavefrontOBJ::IOMode::READ);
    mimmo0->setDir("geodata");
    mimmo0->setFilename("eyeball");
    mimmo0->printResumeFile(false);

    mimmo::IOWavefrontOBJ * mimmo1 = new mimmo::IOWavefrontOBJ(mimmo::IOWavefrontOBJ::IOMode::WRITE);
    mimmo1->setFilename("iogeneric_output_00001.0001");
    mimmo1->printResumeFile(false);
    mimmo1->setTextureUVMode(true);

    /* Instantiation of MRBF object.
     * Set displacements.
     * Plot Optional results during execution active for MRBF block.
     */
    mimmo::MRBF * rbf = new mimmo::MRBF();
    dvecarr3E displ(1);
    displ[0] = {{0.0,0.0,1.0}};
    rbf->addNode(std::array<double,3>({{0.0,0.0,2.0}}));
    rbf->setDisplacements(displ);
    rbf->setFunction(bitpit::RBFBasisFunction::C0C2);
    rbf->setSupportRadiusReal(1.0);
    rbf->setPlotInExecution(true);

    /* Create applier block.
     * It applies the deformation displacements to the original input geometry.
     */
    mimmo::Apply * applier = new mimmo::Apply();


    /* Setup pin connections.
     */
    mimmo::pin::addPin(mimmo0, mimmo1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, mimmo1, M_WAVEFRONTDATA, M_WAVEFRONTDATA);

    mimmo::pin::addPin(mimmo0, rbf, M_GEOM, M_GEOM);

    mimmo::pin::addPin(rbf, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);

    mimmo::pin::addPin(applier, mimmo1, M_GEOM, M_GEOM);

    /* Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(mimmo1);
    ch0.addObject(rbf);
    ch0.addObject(applier);

    /* Execution of chain.
     * Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

    /* Clean up & exit;
     */
    delete mimmo0;
    delete mimmo1;
    delete rbf;
    delete applier;

    return 0;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
        /**<Calling mimic core routines*/
        int err = 1;
        try{
            err = genericinput00005() ;
        }catch(std::exception & e){
            std::cout<<"core_example_00011 exited with an error of type : "<<e.what()<<std::endl;
            return 0;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return err;
}
