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
#include "mimmo_utils.hpp"
#include "mimmo_iogeneric.hpp"
#include "bitpit.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmo::pin;

// =================================================================================== //
/*!
 */

void test0000() {

    int id = 0;
    int dim = 2;
    darray3E origin({-0.5, -0.5, 0.0});
    darray3E lengths({1.0, 1.0, 0.0});
    iarray3E nCells({51, 51, 0});

    bitpit::VolCartesian mesh(id, dim, origin, lengths, nCells);

    mesh.switchMemoryMode(VolCartesian::MEMORY_NORMAL);

    darray3E point = {0.0, 0.0, 0.0};
    long idCell = mesh.locateClosestCell(point);

    bitpit::Cell cell = mesh.getCell(idCell);

    std::unordered_map<long, bool> bounds;

    std::cout << std::endl;
    std::cout << "bounds" << std::endl;
    mimmo::MimmoPiercedVector<darray3E> bcfield;
//    mimmo::MimmoPiercedVector<double> bcfield;

    darray3E value({0.1, 0.1, 0.0});
//    double value(0.1);
    for (int iV = 0; iV < 4; iV++){
        long idV = cell.getVertex(iV);
        bounds[idV] = true;
        bcfield.insert(idV, value);
    }

    value = {0.0, 0.0, 0.0};
//    value = 0.0;
    for (auto inter : mesh.getInterfaces()){
        if (inter.isBorder()){
            for (int iV = 0; iV < 2; iV++){
                long idV = inter.getVertex(iV);
                bounds[idV] = true;
                if (!bcfield.exists(idV)){
                    bcfield.insert(idV, value);
                }
            }
        }
    }

    std::cout << "unstruct" << std::endl;

    bitpit::VolUnstructured meshV(id++, 2);
    bitpit::SurfUnstructured meshS(id++);
    long ID;
    for (auto vertex : mesh.getVertices()){
        ID = vertex.getId();
        meshV.addVertex(vertex, ID);
        if (bounds[ID]){
            meshS.addVertex(vertex, ID);
        }
    }

    for (auto inter : mesh.getInterfaces()){
        if (inter.isBorder() || inter.getOwner() == idCell || inter.getNeigh() == idCell){
            auto const conn = inter.getConnect();
            std::vector<long> vconn(2);
            for (int i=0; i<2; i++)
                vconn[i] = conn[i];
            meshS.addCell(bitpit::ElementInfo::LINE, true, vconn);
        }
    }

    for (auto cell : mesh.getCells()){
        auto const conn = cell.getConnect();
        std::vector<long> vconn(4);
        for (int i=0; i<4; i++)
            vconn[i] = conn[i];
        vconn[2] = conn[3];
        vconn[3] = conn[2];
       meshV.addCell(bitpit::ElementInfo::QUAD, true, vconn);
    }


    std::cout << "mimmo objects" << std::endl;
    meshV.buildAdjacencies();
    meshV.buildInterfaces();
    meshS.buildAdjacencies();
    meshS.buildInterfaces();

    meshS.write();

    MimmoObject* mimmoV = new MimmoObject(2, static_cast<bitpit::PatchKernel*>(&meshV));
    mimmoV->getPatch()->write("mimmo.0000");
    MimmoObject* mimmoS = new MimmoObject(1, static_cast<bitpit::PatchKernel*>(&meshS));
    bcfield.setGeometry(mimmoS);
    bcfield.setDataLocation(1);

    std::cout << "prop" << std::endl;

    PropagateVectorField* prop = new PropagateVectorField();
//    PropagateScalarField* prop = new PropagateScalarField();
    mimmo::setExpertMode(true);
    prop->setGeometry(mimmoV);
    prop->setBoundarySurface(mimmoS);
    prop->setBoundaryConditions(bcfield);
    prop->setWeightConstant(1.0);
    prop->setSmoothingSteps(1000);
    prop->setPlotInExecution(true);
    prop->setApply(true);
    prop->setConvergence(true);

    prop->setDumpingFactor(1.0);
    prop->setDumpingRadius(0.3);

    prop->exec();

}

int	main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
    MPI::Init(argc, argv);

    {
#endif
        /**<Calling mimmo Test routine*/
        test0000() ;

#if ENABLE_MPI==1
    }

    MPI::Finalize();
#endif

    return(1);
}
