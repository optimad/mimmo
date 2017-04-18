/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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
\*---------------------------------------------------------------------------*/
#include "IOVTKScalar.hpp"
#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkTriangleFilter.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!Default constructor of IOVTKScalar.
 */
IOVTKScalar::IOVTKScalar(){
	m_name 		= "mimmo.IOVTKScalar";
	m_read		= false;
	m_rfilename	= "mimmoVTKScalar";
	m_write		= false;
	m_wfilename	= "mimmoVTKScalar";
	m_rdir		= "./";
	m_wdir		= "./";
	m_polydata	= NULL;
	m_local		= false;
	m_normalize	= true;
	m_scaling	= 1.0;
}

/*!Default destructor of IOVTKScalar.
 */
IOVTKScalar::~IOVTKScalar(){
	if (m_local) delete getGeometry();
};

/*!Copy constructor of IOVTKScalar.
 */
IOVTKScalar::IOVTKScalar(const IOVTKScalar & other){
	*this = other;
};

/*!Assignement operator of IOVTKScalar.
 */
IOVTKScalar & IOVTKScalar::operator=(const IOVTKScalar & other){
	m_read = other.m_read;
	m_rfilename = other.m_rfilename;
	m_write = other.m_write;
	m_wfilename = other.m_wfilename;
	m_rdir = other.m_rdir;
	m_wdir = other.m_wdir;
	m_polydata = other.m_polydata;
	m_local = other.m_local;
	m_normalize = other.m_normalize;
	m_scaling = other.m_scaling;
	return *this;
};

/*! It builds the input/output ports of the object
 */
void
IOVTKScalar::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoObject*, IOVTKScalar>(this, &IOVTKScalar::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<double, IOVTKScalar>(this, &IOVTKScalar::setScaling, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<vtkPolyData*, IOVTKScalar>(this, &IOVTKScalar::setPolyData, PortType::M_POLYDATA_, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::POLYDATA_));
	built = (built && createPortIn<dvector1D, IOVTKScalar>(this, &IOVTKScalar::setField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));

	built = (built && createPortOut<MimmoObject*, IOVTKScalar>(this, &IOVTKScalar::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortOut<vtkPolyData*, IOVTKScalar>(this, &IOVTKScalar::getPolyData, PortType::M_POLYDATA_, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::POLYDATA_));
	built = (built && createPortOut<dvector1D, IOVTKScalar>(this, &IOVTKScalar::getField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

vtkPolyData*
IOVTKScalar::getPolyData(){
	return m_polydata;
}

double
IOVTKScalar::getScaling(){
	return m_scaling;
}

dvector1D
IOVTKScalar::getField(){
	return m_field;
}

void
IOVTKScalar::setPolyData(vtkPolyData* polydata){
	m_polydata = polydata;
}

void
IOVTKScalar::setScaling(double scaling){
	m_scaling = scaling;
}

void
IOVTKScalar::setField(dvector1D field){
	m_field = field;
}

void
IOVTKScalar::setNormalize(bool normalize){
	m_normalize = normalize;
}

/*!It sets the condition to read the geometry on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
IOVTKScalar::setRead(bool read){
	m_read = read;
}

/*!It sets the name of directory to read the geometry.
 * \param[in] dir Name of directory.
 */
void
IOVTKScalar::setReadDir(string dir){
	m_rdir = dir;
}

/*!It sets the name of file to read the geometry.
 * \param[in] filename Name of input file.
 */
void
IOVTKScalar::setReadFilename(string filename){
	m_rfilename = filename;
}

/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
IOVTKScalar::setWrite(bool write){
	m_write = write;
}

/*!It sets the name of directory to write the geometry.
 * \param[in] dir Name of directory.
 */
void
IOVTKScalar::setWriteDir(string dir){
	m_wdir = dir;
}

/*!It sets the name of file to write the geometry.
 * \param[in] filename Name of output file.
 */
void
IOVTKScalar::setWriteFilename(string filename){
	m_wfilename = filename;
}


/*!It reads the mesh geometry from an input file.
 * \return False if file doesn't exists or is not a polydata.
 */
bool
IOVTKScalar::read(){

	std::ifstream infile(m_rdir+"/"+m_rfilename+".vtk");
	bool check = infile.good();
	if (!check) return false;
	infile.close();

	if (getGeometry()==NULL){
		//Local Instantiation of mimmo Object.
		m_local = true;
		MimmoObject* mimmo0 = new MimmoObject(1);
		setGeometry(mimmo0);
	}

	string inputFilename = m_rdir+"/"+m_rfilename+".vtk";
	int	np = 0;
	int	nt = 0;
	darray3E point;

	// Get all data from the file
	vtkSmartPointer<vtkGenericDataObjectReader> reader =
			vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName(inputFilename.c_str());
	reader->Update();

	// All of the standard data types can be checked and obtained like this:
	if(reader->IsFilePolyData())
	{
		m_polydata = vtkPolyData::New();
		m_polydata->DeepCopy(reader->GetPolyDataOutput());
		vtkPolyData* output = reader->GetPolyDataOutput();

		vtkSmartPointer<vtkTriangleFilter> tri= vtkSmartPointer<vtkTriangleFilter>::New();
		tri->SetInputData(output);
		tri->Update();

		vtkPolyData* output2 = tri->GetOutput();

		vtkSmartPointer<vtkPoints> points = output2->GetPoints();
		vtkCellArray *cells = output2->GetPolys();
		vtkCellData *cdata = output2->GetCellData();
		vtkPointData *pdata = output2->GetPointData();

		double point_[3];
		for (vtkIdType id=0; id<points->GetNumberOfPoints(); id++ ){
			points->GetPoint(id, point_);
			for (int i=0; i<3; i++) point[i] = point_[i];
			getGeometry()->addVertex(point);
		}

		vtkIdType* conn_;
		vtkIdType npts;
		nt = cells->GetNumberOfCells();
		livector1D connectivity(3);
		bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
		cells->InitTraversal();
		for (vtkIdType id=0; id<nt; id++ ){
			cells->GetNextCell(npts, conn_);
			for (int i=0; i<npts; i++) connectivity[i] = conn_[i];
			getGeometry()->addConnectedCell(connectivity, eltype);
		}
		
		getGeometry()->cleanGeometry();
		np = getGeometry()->getNVertex();

		vtkDataArray* data = pdata->GetArray(0);
		m_field.resize(np);
		if (data != NULL){
			double maxf = 0.0;
			for (int i=0; i<np; i++){
				m_field[i] = (data->GetComponent(i,0));
				maxf = max(maxf, abs(m_field[i]));
			}
			if (m_normalize){
				for (int i=0; i<np; i++){
					m_field[i] /= maxf;
					m_field[i] *= m_scaling;
				}
			}
			else{
				for (int i=0; i<np; i++){
					m_field[i] *= m_scaling;
				}
			}
		}
	}
	else{
		std::cout << "mimmo : ERROR : polydata not found in : "<< m_rfilename << std::endl;
		return false;
	}

	return true;

}

//TODO It writes only if a polydata is linked!! No conversion from MimmoObject to VTK for this moment!!!
/*!It writes the mesh geometry on output .vtu file.
 * It modifies the linked polydata with the linked geometry if exists.
 *\return False if polydata is not linked.
 */
bool
IOVTKScalar::write(){

	if (m_polydata == NULL) return false;

	if (getGeometry() != NULL){
		vtkSmartPointer<vtkPoints> points = m_polydata->GetPoints();
		double point_[3];
		darray3E point;
		int np = points->GetNumberOfPoints();
		np = min(np, int(getGeometry()->getNVertex()));
		for (vtkIdType id=0; id<np; id++ ){
			point = getGeometry()->getVertexCoords(id);
			for (int i=0; i<3; i++) point_[i] = point[i];
			points->SetPoint(id, point_);
		}
	}

	string outputFilename = m_wdir+"/"+m_wfilename+".vtk";
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(outputFilename.c_str());
	writer->SetInputData(m_polydata);
	writer->Write();
	return true;
}

/*!Execution command.
 * It reads the geometry if the condition m_read is true.
 * It writes the geometry if the condition m_write is true.
 */
void
IOVTKScalar::execute(){
	bool check;
	if (m_read) check = read();
	if (!check){
		std::cout << "mimmo : ERROR : file not found : "<< m_rfilename << std::endl;
		std::cout << " " << std::endl;
		exit(10);
	}
	if (m_write) check = write();
	if (!check){
		std::cout << "mimmo : ERROR : write not done : geometry not linked " << std::endl;
		std::cout << " " << std::endl;
		exit(11);
	}
}

}
