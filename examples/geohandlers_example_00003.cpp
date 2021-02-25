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

#include "mimmo_geohandlers.hpp"
#include "mimmo_manipulators.hpp"
#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

// =================================================================================== //
/*!
	\example geohandlers_example_00003.cpp

	\brief Example of usage of subpatch selections from an input geometry,
           geometry field manipulation and final deformation definition.

    Using: MimmoGeometry, SelectionByMap, MRBF, ReconstructVector,
           ExtractField, Apply Partition(MPI version)

	<b>To run</b>              : ./geohandlers_example_00003 \n
    <b>To run (MPI version)</b>: mpirun -np X geohandlers_example_00003 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */

// =================================================================================== //

void test00003() {

    /*
        Read a sphere geometry. Convert mode is to save the just read geometry in
        another file with name geohandlers_output_00003.0000.stl
     */
	mimmo::MimmoGeometry * mimmo0 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo0->setReadDir("geodata");
    mimmo0->setReadFileType(FileType::STL);
    mimmo0->setReadFilename("sphere2");
    mimmo0->setWriteDir("./");
    mimmo0->setWriteFileType(FileType::STL);
    mimmo0->setWriteFilename("geohandlers_output_00003.0000");

    /*
        Read a plane geometry. Convert mode is to save the just read geometry in
        another file with name geohandlers_output_00003p1.0000.stl
     */
    mimmo::MimmoGeometry * mimmo1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo1->setReadDir("geodata");
    mimmo1->setReadFileType(FileType::STL);
    mimmo1->setReadFilename("plane1");
    mimmo1->setWriteDir("./");
    mimmo1->setWriteFileType(FileType::STL);
    mimmo1->setWriteFilename("geohandlers_output_00003p1.0000");

    /*
        Read another plane geometry. Convert mode is to save the just read geometry in
        another file with name geohandlers_output_00003p2.0000.stl
     */
    mimmo::MimmoGeometry * mimmo2 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
    mimmo2->setReadDir("geodata");
    mimmo2->setReadFileType(FileType::STL);
    mimmo2->setReadFilename("plane2");
    mimmo2->setWriteDir("./");
    mimmo2->setWriteFileType(FileType::STL);
    mimmo2->setWriteFilename("geohandlers_output_00003p2.0000");

    /*
        Write the deformed sphere on file
     */
    mimmo::MimmoGeometry * mimmo3 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo3->setWriteDir(".");
    mimmo3->setWriteFileType(FileType::STL);
    mimmo3->setWriteFilename("geohandlers_output_00003.0001");

    /*
        Write the first deformed sub-patch selection on file
     */
    mimmo::MimmoGeometry * mimmo4 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo4->setWriteDir(".");
    mimmo4->setWriteFileType(FileType::STL);
    mimmo4->setWriteFilename("geohandlers_output_00003.0002");

    /*
        Write the second deformed sub-patch selection on file
     */
    mimmo::MimmoGeometry * mimmo5 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    mimmo5->setWriteDir(".");
    mimmo5->setWriteFileType(FileType::STL);
    mimmo5->setWriteFilename("geohandlers_output_00003.0003");

#if MIMMO_ENABLE_MPI
    /*
        Distribute the sphere geometry among processes.
        Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition* partition0 = new mimmo::Partition();
    partition0->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition0->setPlotInExecution(true);

    /*
        Distribute the plane1 geometry among processes.
        Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition* partition1 = new mimmo::Partition();
    partition1->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition1->setPlotInExecution(true);

    /*
        Please note MPI version of this example will not distribute plane2 geometry
        among processes. This is done on purpose, just to prove the code works
        even if the plane2 geometry is retained only on a single rank (that is namely the master rank)
    */
#endif

    /*
      Create two sub-patches of the sphere using proximity selection with the two external
      geometries plane1 and plane2. The sphere cells within a certain offset distance
      (tolerance) from the chosen external geometry will be selected to form the subpatch.
     */
    mimmo::SelectionByMapping  * mapSel1 = new mimmo::SelectionByMapping();
    mapSel1->setTolerance(1.0e-01);
    mapSel1->setPlotInExecution(true);

    mimmo::SelectionByMapping  * mapSel2 = new mimmo::SelectionByMapping();
    mapSel2->setTolerance(1.0e-01);
    mapSel2->setPlotInExecution(true);

    /*
        Creation of a two points with x coordinates -0.5 and 0.5 (y=z=0.0)
        to be used as RBF nodes
     */
    dvecarr3E rbfNodes1(1, {{-0.5, 0.0, 0.0}});
    dvecarr3E rbfNodes2(1, {{0.5, 0.0, 0.0}});

    /*
        Creation of a two displacement list associated to RBF nodes
     */
    dvecarr3E displ1(1, {{-0.25, 0.0, 0.0}});
    dvecarr3E displ2(1, {{0.25, 0.0, 0.0}});

    /*
       Create two RBF manipulators, with nodes and displacements defined before.
       Each manipulator will act on a specific subpatch selected before
     */
    mimmo::MRBF* mrbf1 = new mimmo::MRBF(mimmo::MRBFSol::NONE);
    mrbf1->setSupportRadiusLocal(0.4);
    mrbf1->setPlotInExecution(true);
    mrbf1->setNode(rbfNodes1);
    mrbf1->setDisplacements(displ1);
    mrbf1->setApply(false);

    mimmo::MRBF* mrbf2 = new mimmo::MRBF(mimmo::MRBFSol::NONE);
    mrbf2->setSupportRadiusLocal(0.4);
    mrbf2->setPlotInExecution(true);
    mrbf2->setNode(rbfNodes2);
    mrbf2->setDisplacements(displ2);
    mrbf2->setApply(false);

    //please note apply is set to false, because the two deformation fields onto
    //the subpatches will be combined to form a single deformation field on the
    //whole sphere

    /*
       Recombine the two deformation field calculated on subpatches into a unique
       field defined on the mother sphere mesh. The block will deal with possible
       overlaps and resolve them with its default overlapMethod engine. See doxy documentation.
     */
    mimmo::ReconstructVector* recon = new mimmo::ReconstructVector();
    recon->setPlotInExecution(true);

    /*
       Create an applier block.
     * It applies the reconstructed deformation field onto the sphere.
     */
    mimmo::Apply* applier = new mimmo::Apply();

    /*
       Extract the portion of the whole recostructed deformation field
       referring to the first subpatch. Method of extraction is the ID association
       between mother mesh and sub-patch vertices.
     */
    mimmo::ExtractVectorField* extr = new mimmo::ExtractVectorField();
    extr->setMode(mimmo::ExtractMode::ID);
    extr->setPlotInExecution(true);

    /*
        Create an applier block.
        It will deform the first subpatch with the result of the extr block
    */
    mimmo::Apply* applierextr = new mimmo::Apply();

    /*
       Extract the portion of the whole recostructed deformation field
       referring to the second subpatch. Method of extraction is the MAPPING that is
       data extracted will be imported from the portion of sphere mesh near to the
       second subpatch, within a prescribed proximity distance(tolerance).
     */
    mimmo::ExtractVectorField* extr2 = new mimmo::ExtractVectorField();
    extr2->setMode(mimmo::ExtractMode::MAPPING);
    extr2->setTolerance(1.0e-01);
    extr2->setPlotInExecution(true);

    /*
        Create an applier block.
        It will deform the second subpatch with the result of the extr2 block
    */
    mimmo::Apply* applierextr2 = new mimmo::Apply();

    /*
        Define block pin connections.
     */
#if MIMMO_ENABLE_MPI
    mimmo::pin::addPin(mimmo0, partition0, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition0, mapSel1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition0, mapSel2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition0, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, partition1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(partition1, mapSel1, M_GEOM, M_GEOM2);
#else
    mimmo::pin::addPin(mimmo0, mapSel1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, mapSel2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo0, applier, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mimmo1, mapSel1, M_GEOM, M_GEOM2);
#endif

    mimmo::pin::addPin(mimmo2, mapSel2, M_GEOM, M_GEOM2);

    mimmo::pin::addPin(mapSel1, mrbf1, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mapSel2, mrbf2, M_GEOM, M_GEOM);

    mimmo::pin::addPin(mimmo0, recon, M_GEOM, M_GEOM);
    mimmo::pin::addPin(mrbf1, recon, M_GDISPLS, M_VECTORFIELD);
    mimmo::pin::addPin(mrbf2, recon, M_GDISPLS, M_VECTORFIELD);

    mimmo::pin::addPin(recon, applier, M_VECTORFIELD, M_GDISPLS);
    mimmo::pin::addPin(applier, mimmo3, M_GEOM, M_GEOM);

    mimmo::pin::addPin(mapSel1, extr, M_GEOM, M_GEOM);
    mimmo::pin::addPin(recon, extr, M_VECTORFIELD, M_VECTORFIELD);
    mimmo::pin::addPin(mapSel1, applierextr, M_GEOM, M_GEOM);
    mimmo::pin::addPin(extr, applierextr, M_VECTORFIELD, M_GDISPLS);
    mimmo::pin::addPin(applierextr, mimmo4, M_GEOM, M_GEOM);

    mimmo::pin::addPin(mapSel2, extr2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(recon, extr2, M_VECTORFIELD, M_VECTORFIELD);
    mimmo::pin::addPin(mapSel2, applierextr2, M_GEOM, M_GEOM);
    mimmo::pin::addPin(extr2, applierextr2, M_VECTORFIELD, M_GDISPLS);
    mimmo::pin::addPin(applierextr2, mimmo5, M_GEOM, M_GEOM);

    /*
        Setup execution chain.
     */
    mimmo::Chain ch0;
    ch0.addObject(mimmo0);
    ch0.addObject(mimmo1);
    ch0.addObject(mimmo2);
#if MIMMO_ENABLE_MPI
    ch0.addObject(partition0);
    ch0.addObject(partition1);
#endif
    ch0.addObject(mapSel1);
    ch0.addObject(mapSel2);
    ch0.addObject(mrbf1);
    ch0.addObject(mrbf2);
    ch0.addObject(recon);
    ch0.addObject(applier);
    ch0.addObject(mimmo3);
    ch0.addObject(extr);
    ch0.addObject(applierextr);
    ch0.addObject(mimmo4);
    ch0.addObject(extr2);
    ch0.addObject(applierextr2);
    ch0.addObject(mimmo5);

    /*
        Execute the chain.
        Use debug flag true to to print out the execution steps.
     */
    ch0.exec(true);

    /*
        Clean up & exit;
     */
    delete mimmo0;
    delete mimmo1;
#if MIMMO_ENABLE_MPI
    delete partition0;
    delete partition1;
#endif
    delete mimmo2;
    delete mimmo3;
    delete mimmo4;
    delete mimmo5;
    delete mapSel1;
    delete mapSel2;
    delete mrbf1;
    delete mrbf2;
    delete recon;
    delete extr;
    delete extr2;
    delete applier;
    delete applierextr;
    delete applierextr2;

	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif
		/**<Calling core function*/
		try{
            test00003();
        }
        catch(std::exception & e){
            std::cout<<"geohandlers_example_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return 0;
}
