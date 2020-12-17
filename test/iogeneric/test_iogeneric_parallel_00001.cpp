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

#include "mimmo_iogeneric.hpp"
#include "Partition.hpp"
#include <exception>


// =================================================================================== //
/*!
 * //testing parallel :
   - reading vtu ascii with master rank
   - partition mesh and writing it in parallel vtu binary
   - reread parallel vtu binary and write it in parallel vtu ascii
 */
int test1() {

    mimmo::MimmoGeometry * reader = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    bitpit::Logger & log = reader->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);

    reader->setReadDir("geodata");
    reader->setReadFilename("curve");
    reader->setReadFileType(FileType::CURVEVTU);
    reader->setCodex(0);
    reader->exec();

    log<<"Read single Curve.vtu with master rank"<<std::endl;

    mimmo::Partition * part = new mimmo::Partition();
    part->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    part->setGeometry(reader->getGeometry());
    part->setPlotInExecution(false);
    part->exec();

    log<<"Partitioned curve"<<std::endl;

    mimmo::MimmoGeometry * writer = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    writer->setName("WriteCurve-Partitioned-binary");
    writer->setGeometry(part->getGeometry());
    writer->setWriteDir(".");
    writer->setWriteFilename("tip1_curveBinary");
    writer->setWriteFileType(FileType::CURVEVTU);
    writer->setCodex(1);
    writer->exec();

    log<<"Written curve in *pvtu binary format (if np > 1)"<<std::endl;

    mimmo::MimmoGeometry * readbin = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    readbin->setReadDir(".");
    readbin->setReadFilename("tip1_curveBinary");
    readbin->setReadFileType(FileType::CURVEVTU);
    readbin->setCodex(1);
    readbin->exec();

    mimmo::MimmoGeometry * writeAscii = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
    writeAscii->setGeometry(readbin->getGeometry());
    writeAscii->setWriteDir(".");
    writeAscii->setWriteFilename("tip1_curveASCII");
    writeAscii->setWriteFileType(FileType::CURVEVTU);
    writeAscii->setCodex(0);
    writeAscii->exec();

    log<<"Converted curve in *pvtu ascii format (if np > 1)"<<std::endl;

    bool check = writeAscii->getGeometry()->getNGlobalCells() == 412;
    check = check && writeAscii->getGeometry()->getNGlobalVertices() == 414;
    log<<"test passed : "<<check<<std::endl;


    delete reader;
    delete part;
    delete writer;
    delete readbin;
    delete writeAscii;
    return int(!check);

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
    int val = 1;
    try{
        /**<Calling mimmo Test routines*/
        val = test1() ;
    }
    catch(std::exception & e){
        std::cout<<"test_iogeneric_parallel_00001 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
