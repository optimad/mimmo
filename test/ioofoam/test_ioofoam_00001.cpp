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

#include "IOOFOAM.hpp"
#include "Chain.hpp"
#include <exception>

// =================================================================================== //
int test1() {

    mimmo::IOOFOAM * reader = new mimmo::IOOFOAM(false);
    reader->setDir("geodata/OFOAM");

    mimmo::IOOFOAMScalarField * fieldreader = new mimmo::IOOFOAMScalarField(false);
    fieldreader->setDir("geodata/OFOAM");
    fieldreader->setFieldName("p");

    mimmo::pin::addPin(reader, fieldreader, M_GEOMOFOAM, M_GEOMOFOAM);
    mimmo::pin::addPin(reader, fieldreader, M_GEOMOFOAM2, M_GEOMOFOAM2);
    mimmo::pin::addPin(reader, fieldreader, M_UMAPIDS, M_UMAPIDS);

    mimmo::Chain c0;
    c0.addObject(reader);
    c0.addObject(fieldreader);
    c0.exec(true);

    bool check = true;
    check = check && (reader->getGeometry()->getPatch()->getVertexCount() == 120450);
    check = check && (fieldreader->getBoundaryGeometry()->getPatch()->getVertexCount() == 120450);
    check = check && (reader->getGeometry()->getPatch()->getCellCount() == 59540);
    check = check && (fieldreader->getBoundaryGeometry()->getPatch()->getCellCount() == 120448);

    check = check && (reader->getGeometry() == fieldreader->getGeometry());
    check = check && (reader->getBoundaryGeometry() == fieldreader->getBoundaryGeometry());

    auto field = fieldreader->getBoundaryField();

    check = check && (field->getGeometry() == fieldreader->getBoundaryGeometry());

    double maxval = std::numeric_limits<double>::min();;
    double minval = std::numeric_limits<double>::max();
    for(auto it=field->begin(); it!=field->end(); ++it){
        maxval = std::max(maxval, *it);
        minval = std::min(minval, *it);
    }

    std::cout<<maxval<<"  "<<minval<<std::endl;
    check = check && std::abs(std::abs(maxval)-0.14221) <=1.0E-5;
    check = check && std::abs(std::abs(minval)-0.08285) <=1.0E-5;

    delete reader;
    delete fieldreader;

    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

    #if MIMMO_ENABLE_MPI
    	MPI_Init(&argc, &argv);
    #endif
    		/**<Calling mimmo Test routines*/
            int val = 1;
            try{
                val = test1() ;
            }
            catch(std::exception & e){
                std::cout<<"test_ioofoam_00001 exited with an error of type : "<<e.what()<<std::endl;
                return 1;
            }

    #if MIMMO_ENABLE_MPI
    	MPI_Finalize();
    #endif

    return val;
}
