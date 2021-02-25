/*----------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
 *  Copyright (c) 2015-2021 Optimad Engineering S.r.l., All Rights Reserved.
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

    \brief Example of usage of IOWavefrontOBJ to parse an OBJ mesh and deform it
           with a RBF manipulator

    Using: IOWavefrontOBJ, MRBF, Apply, Chain

    <b>To run</b>              : ./genericinput_example_00005 \n
    <b>To run (MPI version)</b>: mpirun -np X genericinput_example_00005 \n

    <b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


int genericinput00005() {

    /*
     * Set expert mode true
     */
    mimmo::setExpertMode();

    /*
        Create a block to parse an OBJ mesh from file.
     */
    mimmo::IOWavefrontOBJ * mimmo0 = new mimmo::IOWavefrontOBJ(mimmo::IOWavefrontOBJ::IOMode::READ);
    mimmo0->setDir("geodata");
    mimmo0->setFilename("eyeball");
    mimmo0->printResumeFile(false);

    /*
        Create a block to write an OBJ mesh to file.
     */
    mimmo::IOWavefrontOBJ * mimmo1 = new mimmo::IOWavefrontOBJ(mimmo::IOWavefrontOBJ::IOMode::WRITE);
    mimmo1->setFilename("iogeneric_output_00001.0001");
    mimmo1->printResumeFile(false);
    mimmo1->setTextureUVMode(true);

    /*
        Define a single node position with its displacement to
        build up the input needed by the RBF manipulator
     */
    std::array<double,3>node({{0.0,0.0,2.0}});
    dvecarr3E displ(1);
    displ[0] = {{0.0,0.0,1.0}};

    /*
        Create the RBF manipulator.
     */
    mimmo::MRBF * rbf = new mimmo::MRBF();
    rbf->addNode(node);
    rbf->setDisplacements(displ);
    rbf->setFunction(bitpit::RBFBasisFunction::C0C2);
    rbf->setSupportRadiusReal(1.0);
    rbf->setPlotInExecution(true);

    /*
       Create the applier block.
     * It applies the deformation displacements from rbf block to the original input geometry.
     */
    mimmo::Apply * applier = new mimmo::Apply();


    /*
        Define block pin connections.
     */
    mimmo::pin::addPin(mimmo0, mimmo1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, mimmo1, M_WAVEFRONTDATA, M_WAVEFRONTDATA);
    mimmo::pin::addPin(mimmo0, rbf, M_GEOM, M_GEOM);

    mimmo::pin::addPin(rbf, applier, M_GDISPLS, M_GDISPLS);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);

    mimmo::pin::addPin(applier, mimmo1, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(mimmo1);
    ch0.addObject(rbf);
    ch0.addObject(applier);

    /*
        Execute the chain.
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
        /**<Calling core function*/
        int err = 1;
        try{
            err = genericinput00005() ;
        }catch(std::exception & e){
            std::cout<<"genericinput_example_00005 exited with an error of type : "<<e.what()<<std::endl;
            return 0;
        }
#if MIMMO_ENABLE_MPI
    MPI_Finalize();
#endif

    return err;
}
