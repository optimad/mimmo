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
\*---------------------------------------------------------------------------*/

#include "CreatePointCloud.hpp"

namespace mimmo {

/*!
 * Default constructor of CreatePointCloud.
 */
CreatePointCloud::CreatePointCloud(){
    m_name         = "mimmo.CreatePointCloud";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
CreatePointCloud::CreatePointCloud(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.CreatePointCloud";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);

    if(input == "mimmo.CreatePointCloud"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}


CreatePointCloud::~CreatePointCloud(){};

/*!
 * Copy constructor of CreatePointCloud. Only sensible inputs are copied.
 */
CreatePointCloud::CreatePointCloud(const CreatePointCloud & other):BaseManipulation(other){
    m_rawpoints    = other.m_rawpoints;
    m_rawvector    = other.m_rawvector;
    m_rawscalar    = other.m_rawscalar;
};

/*!Assignement operator of CreatePointCloud.
 */
CreatePointCloud & CreatePointCloud::operator=(CreatePointCloud other){
    swap(other);
    return *this;
};

/*!
 * Swap method
 * \param[in] x object to be swapped
 */
void CreatePointCloud::swap(CreatePointCloud & x) noexcept
{
    m_rawpoints.swap(x.m_rawpoints);
    m_rawscalar.swap(x.m_rawscalar);
    m_rawvector.swap(x.m_rawvector);
    m_scalarfield.swap(x.m_scalarfield);
    m_vectorfield.swap(x.m_vectorfield);

    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
CreatePointCloud::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dvecarr3E, CreatePointCloud>(this, &mimmo::CreatePointCloud::setRawPoints, M_COORDS, true,1));
    built = (built && createPortIn<dvecarr3E, CreatePointCloud>(this, &mimmo::CreatePointCloud::setRawVectorField, M_DISPLS));
    built = (built && createPortIn<dvector1D, CreatePointCloud>(this, &mimmo::CreatePointCloud::setRawScalarField, M_DATAFIELD));

    built = (built && createPortIn<dmpvecarr3E*, CreatePointCloud>(this, &mimmo::CreatePointCloud::setRawPoints, M_VECTORFIELD, true,1));
    built = (built && createPortIn<dmpvecarr3E*, CreatePointCloud>(this, &mimmo::CreatePointCloud::setRawVectorField, M_VECTORFIELD2));
    built = (built && createPortIn<dmpvector1D*, CreatePointCloud>(this, &mimmo::CreatePointCloud::setRawScalarField, M_SCALARFIELD));

    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, CreatePointCloud>(this, &mimmo::CreatePointCloud::getGeometry, M_GEOM));
    built = (built && createPortOut<dmpvector1D*, CreatePointCloud>(this, &mimmo::CreatePointCloud::getScalarField, M_SCALARFIELD));
    built = (built && createPortOut<dmpvecarr3E*, CreatePointCloud>(this, &mimmo::CreatePointCloud::getVectorField, M_VECTORFIELD));

    m_arePortsBuilt = built;
}

/*!
 * Return the scalar field stored in the class as pointer to MimmoPiercedVector object.
 * \return scalar field stored in the class
 */
dmpvector1D*
CreatePointCloud::getScalarField(){
    return &m_scalarfield;
};

/*!
 * Return the vector field stored in the class as pointer to MimmoPiercedVector object.
 * \return vector field stored in the class
 */
dmpvecarr3E*
CreatePointCloud::getVectorField(){
    return &m_vectorfield;
};


/*!
 * It sets the point cloud raw list
 * \param[in] rawPoints raw list of points.
 */
void
CreatePointCloud::setRawPoints(dmpvecarr3E * rawPoints){
    if(rawPoints)    m_rawpoints = *rawPoints;
};

/*!
 * It sets the point cloud raw list
 * \param[in] rawPoints raw list of points.
 */
void
CreatePointCloud::setRawPoints(dvecarr3E rawPoints){
    m_rawpoints.clear();
    m_rawpoints.reserve(rawPoints.size());
    long count(0);
    for(std::array<double,3> &val : rawPoints){
        m_rawpoints.insert(count, val);
        ++count;
    }
};

/*!
 * It sets the vector field optionally attached to point cloud raw list
 * \param[in] rawVectorField raw vector field.
 */
void
CreatePointCloud::setRawVectorField(dmpvecarr3E *rawVectorField){
    if(rawVectorField)  m_rawvector = *rawVectorField;
};

/*!
 * It sets the vector field optionally attached to point cloud raw list
 * \param[in] rawVectorField raw vector field.
 */
void
CreatePointCloud::setRawVectorField(dvecarr3E rawVectorField){
    m_rawvector.clear();
    m_rawvector.reserve(rawVectorField.size());
    long count(0);
    for(std::array<double,3> &val : rawVectorField){
        m_rawvector.insert(count, val);
        ++count;
    }
};

/*!
 * It sets the scalar field optionally attached to point cloud raw list
 * \param[in] rawScalarField raw scalar field.
 */
void
CreatePointCloud::setRawScalarField(dmpvector1D * rawScalarField){
    if(rawScalarField)  m_rawscalar = *rawScalarField;
};

/*!
 * It sets the scalar field optionally attached to point cloud raw list
 * \param[in] rawScalarField raw scalar field.
 */
void
CreatePointCloud::setRawScalarField(dvector1D rawScalarField){
    m_rawscalar.clear();
    m_rawscalar.reserve(rawScalarField.size());
    long count(0);
    for(double &val : rawScalarField){
        m_rawscalar.insert(count, val);
        ++count;
    }
};

/*!
 * Clear all data stored in the class
 */
void
CreatePointCloud::clear(){
    m_vectorfield.clear();
    m_scalarfield.clear();
    m_rawscalar.clear();
    m_rawvector.clear();
    m_rawpoints.clear();
    BaseManipulation::clear();
}

/*!
 * Execution command.
 * Read data from or Write data on linked filename
 */
void
CreatePointCloud::execute(){

    long np = m_rawpoints.size();

#if MIMMO_ENABLE_MPI
    // Rawpoints can be shared among all procs or retained by 0 only.
    // The condition is that rank 0 must have the input points, because
    // the geometry is filled only by master processor 0. The only useful
    // points are those owned by rank 0.
    // Check if rank 0 owns some points (not communicated to the other processors).
    MPI_Bcast(&np, 1, MPI_INT, 0, m_communicator);
#endif

    if(np == 0){
        (*m_log)<< "Warning in "<<m_name<<" : no raw points to work on"<<std::endl;
    }

    //fill the point cloud
    {
        MimmoSharedPointer<MimmoObject> dum(new MimmoObject(3));
        m_geometry = dum;
    }
#if MIMMO_ENABLE_MPI
    //leave the filling to the master rank
    if(getRank() == 0)
#endif
    {
        m_geometry->getPatch()->reserveVertices(np);
        m_geometry->getPatch()->reserveCells(np);

        long id;
        for(auto it= m_rawpoints.begin(); it != m_rawpoints.end(); ++it){
            id = it.getId();
            m_geometry->addVertex(*it, id);
            m_geometry->addConnectedCell(livector1D(1,id), bitpit::ElementType::VERTEX, 0, id, 0);
        }
    }

    m_scalarfield.initialize(m_geometry,MPVLocation::POINT, 0.);
    m_scalarfield.setName("PCScalar");

    m_vectorfield.initialize(m_geometry,MPVLocation::POINT, {{0.,0.,0.}});
    m_vectorfield.setName("PCVector");

    //fill data attached
#if MIMMO_ENABLE_MPI
    //leave the filling to the master rank
    if(getRank() == 0)
#endif
    {
        for(auto it=m_rawscalar.begin(); it!=m_rawscalar.end(); ++it){
            if(m_scalarfield.exists(it.getId()))    m_scalarfield[it.getId()] = *it;
        }
        for(auto it=m_rawvector.begin(); it!=m_rawvector.end(); ++it){
            if(m_vectorfield.exists(it.getId()))    m_vectorfield[it.getId()] = *it;
        }
    }

    m_geometry->update();
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
CreatePointCloud::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);
}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
CreatePointCloud::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

};


/*!
 * Plot cloud point and store it in *.vtu file
 */
void
CreatePointCloud::plotOptionalResults(){

    if (getGeometry() == nullptr) return;
    BaseManipulation::write(getGeometry(), m_scalarfield, m_vectorfield);
}

}
