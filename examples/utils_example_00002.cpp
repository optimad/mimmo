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
#if MIMMO_ENABLE_MPI
    #include "Partition.hpp"
#endif

// =================================================================================== //
/*!
	\example utils_example_00002.cpp

	\brief Example of calculating penetrations with ControlDeformExtSurface. The
    Stanford Bunny fit the open box?

	Using: MimmoGeometry, ControlDeformExtSurface

	<b>To run</b>: ./utils_example_00002 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


void test00002() {

    /* reading target stanford Bunny
     */
    mimmo::MimmoGeometry * bunny = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    bitpit::Logger & log = bunny->getLog();
    bunny->setName("ReaderStanfordBunny");
    bunny->setReadDir("geodata");
    bunny->setReadFileType(FileType::STL);
    bunny->setReadFilename("stanfordBunny2");
    bunny->exec();

    log.setPriority(bitpit::log::Priority::NORMAL);
    log<<"Target Bunny read"<<std::endl;
    log.setPriority(bitpit::log::Priority::DEBUG);
    /* reading constraint inclinedPlane
     */
    mimmo::MimmoGeometry * incplane = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    incplane->setName("ReaderInclinedPlaneConstraint");
    incplane->setReadDir("geodata");
    incplane->setReadFileType(FileType::STL);
    incplane->setReadFilename("inclinedPlane");
    incplane->exec();

    log.setPriority(bitpit::log::Priority::NORMAL);
    log<<"Inclined plane constraint read"<<std::endl;
    log.setPriority(bitpit::log::Priority::DEBUG);

#if MIMMO_ENABLE_MPI
     /* Partitioner of the bunny.
     */
     mimmo::Partition* partbunny = new mimmo::Partition();
     partbunny->setName("PartitionerBunny");
     partbunny->setGeometry(bunny->getGeometry());
     partbunny->setPlotInExecution(true);
     partbunny->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
     partbunny->exec();

     log.setPriority(bitpit::log::Priority::NORMAL);
     log<<"Target Bunny partitioned"<<std::endl;
     log.setPriority(bitpit::log::Priority::DEBUG);

     /* Partitioner of the inclined plane.
     */
     mimmo::Partition* partincplane = new mimmo::Partition();
     partincplane->setName("PartitionerInclinedPlane");
     partincplane->setGeometry(incplane->getGeometry());
     partincplane->setPlotInExecution(true);
     partincplane->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
     partincplane->exec();

     log.setPriority(bitpit::log::Priority::NORMAL);
     log<<"Inclined plane constraint partitioned"<<std::endl;
     log.setPriority(bitpit::log::Priority::DEBUG);

#endif


    /*Creation of violation/penetration field calculator
     */
#if MIMMO_ENABLE_MPI
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> target = partbunny->getGeometry();
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> constraint = partincplane->getGeometry();
#else
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> target = bunny->getGeometry();
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> constraint = incplane->getGeometry();
#endif

    mimmo::dmpvecarr3E def;
    def.initialize(target, mimmo::MPVLocation::POINT, {{0.0,0.0,0.0}});

    mimmo::ControlDeformExtSurface * cdes = new mimmo::ControlDeformExtSurface();
    cdes->setGeometry(target);
    cdes->addConstraint(constraint);
    cdes->addConstraintFile("geodata/openBox.stl", FileType::STL);
    cdes->setDefField(&def);
    cdes->setTolerance(0.0);
    cdes->setPlotInExecution(true);
    cdes->exec();

    log.setPriority(bitpit::log::Priority::NORMAL);
    log<<"Violation field calculated"<<std::endl;
    log.setPriority(bitpit::log::Priority::DEBUG);


    /* Clean up & exit;
     */
    delete bunny;
    delete incplane;
    delete cdes;
#if MIMMO_ENABLE_MPI
    delete partbunny;
    delete partincplane;
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
            test00002() ;
        }
        catch(std::exception & e){
            std::cout<<"utils_example_00002 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;
}
