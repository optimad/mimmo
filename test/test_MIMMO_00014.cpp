/*---------------------------------------------------------------------------*\
 * 
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

#include "bitpit.hpp"
#include "MiMMO.hpp"
#include "BvTree.hpp"
#include <functional>
#include <string>
#include <cstring>

using namespace std;
using namespace bitpit;
using namespace mimmo;

#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test00014(MimmoGeometry * mimmo0, MimmoGeometry * mimmo1 ) {

	//Creation of mimmo container.
	MimmoObject * object_ = mimmo0->getGeometry();


	//Initialize points
//	Lattice * mesh = new Lattice();
	{
//		darray3E origin =  {{0.0,0.0, 0.0}};
//		darray3E span = {{1.2,1.2,1.2}};

//		//spoiler
//		darray3E origin =  {{3.204,0.0, 1.0}};
//		darray3E span = {{0.2,1.1,0.2}};

		//		ahmed
//		darray3E origin =  {{-0.52, 0.0, 0.168}};
//		darray3E span = {{1.1, 0.5, 0.4}};

		//setting mesh
//		mesh->setShape(ShapeType::CUBE);
//		mesh->setOrigin(origin);
//		mesh->setSpan(span[0],span[1],span[2]);
//		mesh->setDimension(iarray3E{{10, 10, 10}});
	}

//	mesh->exec();
//
//	mesh->plotGrid(".", "Box_0006", 0, true);


	cout << "create bvtree" << endl;
	//create bvtree object
	bitpit::PatchKernel *patch_ = object_->getPatch();
	bitpit::SurfUnstructured* spatch_ = static_cast<bitpit::SurfUnstructured*>(patch_);
	BvTree tree(patch_);

	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();
	tree.setMaxLeafSize(1);
	tree.buildTree();
	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "mimmo build of bv-tree took me " << time_span.count() << " seconds.";
	std::cout << std::endl;

	cout << "bv-tree nodes		: " << tree.m_nnodes << endl;
	cout << "bv-tree elements	: " << tree.m_nelements << endl;
	cout << "bv-tree leaf nodes	: " << tree.m_nleaf << endl;

	cout << "execution start" << endl;
	t1 = steady_clock::now();

/*	int nP = mesh->getNNodes();
	dvecarr3E 	points = mesh->getGlobalCoords();
	dvector1D dF(nP);
	darray3E normal;
	long id;
	double r;
	for (int i=0; i<nP; i++){
		r = 0.01;
		dF[i] = signedDistance(&points[i], &tree, id, normal, r, spatch_);
	}

	r = 0.01;
	dvecarr3E Ppoints = projectPoint(&points, &tree, r);


	t2 = steady_clock::now();
	cout << "execution done" << endl;
	//Print execution time
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "mimmo distance and projection took me " << time_span.count() << " seconds.";
	std::cout << std::endl;




	{
		dvecarr3E * controlNodes = &points;
		ivector1D conn(controlNodes->size());
		for (int i=0; i<conn.size(); i++){
			conn[i] = i;
		}

		string nfield = "dF";
		bitpit::VTKUnstructuredGrid output( "./", "PointsDF_0006.0000", bitpit::VTKElementType::VERTEX, *controlNodes, conn);
		output.addData( nfield, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, dF ) ;
		output.write() ;

		controlNodes = &Ppoints;

		bitpit::VTKUnstructuredGrid output1( "./", "PointsDF_0006.0001", bitpit::VTKElementType::VERTEX, *controlNodes, conn);
		output1.addData( nfield, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, dF ) ;
		output1.write() ;
	}



	mesh->plotGridScalar("./", "GridDF_0006", 0, false, dF);
	mesh->plotGridScalar("./", "GridDF_0006", 1, false, dF);

*/


	//Selection by patch
	MimmoObject * objectsel_ = mimmo1->getGeometry();
	bitpit::PatchKernel *patchsel_ = objectsel_->getPatch();
	BvTree treesel(patchsel_);
	treesel.buildTree();

	cout << "execution start" << endl;
	t1 = steady_clock::now();

	double tol = 0.000001;
	std::vector<long> selected = mimmo::bvTreeUtils::selectByPatch(&treesel, &tree, tol);

	t2 = steady_clock::now();
	cout << "execution done" << endl;
	//Print execution time
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "mimmo selection took me " << time_span.count() << " seconds.";
	std::cout << std::endl;



	//create your subpatch.
	std::unique_ptr<MimmoObject> temp(new MimmoObject(1));
	bitpit::PatchKernel * tri = object_->getPatch();

	bitpit::PiercedVector<bitpit::Vertex> & mapV = temp->getPatch()->getVertices();

	livector1D TT;

	long idV;
	int sizeCC;
	int counter=0;
	bitpit::ElementInfo::Type eltype;

	for(auto && idCell : selected){

		bitpit::Cell & cell = tri->getCell(idCell);
		eltype = cell.getType();
		sizeCC = cell.getVertexCount();
		TT.resize(sizeCC);

		for(int i=0; i<sizeCC; ++i){
			idV = cell.getVertex(i);
			TT[i] = idV;

			if(!mapV.exists(idV))	temp->addVertex(tri->getVertexCoords(idV),idV);
		}
		temp->addConnectedCell(TT,eltype,idCell);
		TT.clear();
		counter++;
	}

	tri->deleteCells(selected, true, true);

	tri = NULL;


	temp->getPatch()->write("selection.0");
	object_->getPatch()->write("selection.1");

//	delete mesh;
//	mesh = NULL;


}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling mimmo Test routines*/

		MimmoGeometry * mimmo0 = new MimmoGeometry();
		mimmo0->setRead(true);
		mimmo0->setReadDir("geo_data");
		mimmo0->setReadFileType(FileType::STL);
//		mimmo0->setReadFilename("sphere2");
//		mimmo0->setReadFilename("cube");
//		mimmo0->setReadFilename("spoiler2");
//		mimmo0->setReadFilename("ahmed");
		mimmo0->setReadFilename("drivAerBin2");
		mimmo0->exec();

		MimmoGeometry * mimmo1 = new MimmoGeometry();
		mimmo1->setRead(true);
		mimmo1->setReadDir("geo_data");
		mimmo1->setReadFileType(FileType::STL);
//		mimmo1->setReadFilename("sphere");
//		mimmo1->setReadFilename("cube");
//		mimmo1->setReadFilename("spoiler2");
//		mimmo1->setReadFilename("ahmed");
		mimmo1->setReadFilename("prese");
		mimmo1->exec();

		test00014(mimmo0, mimmo1);

		delete mimmo0;
		mimmo0 = NULL;
		delete mimmo1;
		mimmo1 = NULL;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif

	return 0;
}
