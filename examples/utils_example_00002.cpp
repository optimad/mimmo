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
#include "mimmo_iogeneric.hpp"

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
    mimmo::MimmoGeometry * bunny = new mimmo::MimmoGeometry();

    bunny->setIOMode(IOMode::READ);
    bunny->setReadDir("geodata");
    bunny->setReadFileType(FileType::STL);
    bunny->setReadFilename("stanfordBunny2");
    bunny->exec();

    /*Creation of violation/penetration field calculator
     */
    mimmo::ControlDeformExtSurface * cdes = new mimmo::ControlDeformExtSurface();
    cdes->setGeometry(bunny->getGeometry());

    mimmo::dmpvecarr3E def(bunny->getGeometry(), mimmo::MPVLocation::POINT);
    def.reserve(bunny->getGeometry()->getPatch()->getVertexCount());
    for(auto it=bunny->getGeometry()->getPatch()->vertexBegin(); it!=bunny->getGeometry()->getPatch()->vertexEnd(); ++it){
        def.insert(it.getId(), {{0.0,0.0,0.0}});
    }
    cdes->setDefField(&def);
    cdes->addFile("geodata/openBox.stl", 0.0, FileType::STL);
    cdes->setBackgroundDetails();
    cdes->setPlotInExecution(true);
    cdes->exec();

    /* Clean up & exit;
     */
    delete bunny;
    delete cdes;
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
