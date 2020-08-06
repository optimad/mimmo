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
#include "SpecularPoints.hpp"
#include "SkdTreeUtils.hpp"

namespace mimmo{

/*!
 * Default constructor of SpecularPoints
 */
SpecularPoints::SpecularPoints(){
    m_name = "mimmo.SpecularPoints";
    m_plane.fill(0.0);
    m_insideout = false;
    m_force = true;
    m_origin.fill(0.0);
    m_normal.fill(0.0);
    m_implicit = false;
    m_scalar = nullptr;
    m_vector = nullptr;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SpecularPoints::SpecularPoints(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.SpecularPoints";
    m_plane.fill(0.0);
    m_insideout = false;
    m_force = true;
    m_origin.fill(0.0);
    m_normal.fill(0.0);
    m_implicit = false;
    m_scalar = nullptr;
    m_vector = nullptr;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.SpecularPoints"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of SpecularPoints
 */
SpecularPoints::~SpecularPoints(){};

/*!Copy constructor of SpecularPoints. No result are copied.
 */
SpecularPoints::SpecularPoints(const SpecularPoints & other):ProjPatchOnSurface(other){
    m_plane = other.m_plane;
    m_scalar = other.m_scalar;
    m_vector = other.m_vector;
    m_insideout = other.m_insideout;
    m_force = other.m_force;
    m_origin = other.m_origin;
    m_normal = other.m_normal;
    m_implicit = other.m_implicit;
    m_pc = other.m_pc;
};

/*!
 * Assignment operator. No result are copied
 */
SpecularPoints & SpecularPoints::operator=(SpecularPoints other){
    swap(other);
    return *this;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void SpecularPoints::swap(SpecularPoints & x) noexcept
{
    std::swap(m_plane, x.m_plane);
    std::swap(m_scalar, x.m_scalar);
    std::swap(m_vector, x.m_vector);
    std::swap(m_insideout, x.m_insideout);
    std::swap(m_force, x.m_force);
    std::swap(m_origin, x.m_origin);
    std::swap(m_normal, x.m_normal);
    std::swap(m_implicit, x.m_implicit);
    std::swap(m_scalarMirrored, x.m_scalarMirrored);
    std::swap(m_vectorMirrored, x.m_vectorMirrored);
    std::swap(m_pc,x.m_pc);

    ProjPatchOnSurface::swap(x);
}
/*! It builds the input/output ports of the object
 */
void
SpecularPoints::buildPorts(){
    bool built = true;

    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, SpecularPoints>(this,&mimmo::SpecularPoints::setPointCloud, M_GEOM2, true));
    built = (built && createPortIn<dmpvecarr3E*, SpecularPoints>(this, &mimmo::SpecularPoints::setVectorData, M_VECTORFIELD));
    built = (built && createPortIn<dmpvector1D*, SpecularPoints>(this, &mimmo::SpecularPoints::setScalarData, M_DATAFIELD));
    built = (built && createPortIn<darray4E, SpecularPoints>(this, &mimmo::SpecularPoints::setPlane, M_PLANE));
    built = (built && createPortIn<darray3E, SpecularPoints>(this, &mimmo::SpecularPoints::setOrigin, M_POINT));
    built = (built && createPortIn<darray3E, SpecularPoints>(this, &mimmo::SpecularPoints::setNormal, M_AXIS));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, SpecularPoints>(this, &mimmo::SpecularPoints::setGeometry, M_GEOM));

    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, SpecularPoints>(this, &mimmo::SpecularPoints::getMirroredPointCloud, M_GEOM));
    built = (built && createPortOut<dmpvecarr3E*, SpecularPoints>(this, &mimmo::SpecularPoints::getMirroredVectorData, M_VECTORFIELD));
    built = (built && createPortOut<dmpvector1D*, SpecularPoints>(this, &mimmo::SpecularPoints::getMirroredScalarData, M_SCALARFIELD));
    m_arePortsBuilt = built;
};

/*!
 * Returns pointer to the original scalar data field attached to cloud points
 * \returns input scalar field
 */
dmpvector1D *
SpecularPoints::getOriginalScalarData(){
    return m_scalar;
};

/*!
 * Returns pointer to the original vector data field attached to cloud points
 * \returns input vector field
 */
dmpvecarr3E *
SpecularPoints::getOriginalVectorData(){
    return m_vector;
};

/*!
 * Returns pointer to the scalar data field attached to mirrored cloud points
 * \returns output scalar field
 */
dmpvector1D *
SpecularPoints::getMirroredScalarData(){
    return &m_scalarMirrored;
};

/*!
 * Returns pointer to the vector data field attached to mirrored cloud points
 * \returns output vector field
 */
dmpvecarr3E *
SpecularPoints::getMirroredVectorData(){
    return &m_vectorMirrored;
};


/*!
    It gets the coordinates of the mirrored point cloud. Meaningful after execution.
 * \return raw list of coordinates of local point cloud.
 */
dvecarr3E
SpecularPoints::getMirroredRawCoords(){

    dvecarr3E points;
    if(m_patch != nullptr){
        points.reserve(m_patch->getNVertices());
        for(bitpit::Vertex & vert : m_patch->getVertices()){
            points.push_back(vert.getCoords());
        }
    }
    return points;
};

/*!
    It gets the labels associated to the list of getMirroredRawCoords. Meaningful after execution.
 * \return labels of mirrored raw point cloud coordinates.
 */
livector1D
SpecularPoints::getMirroredLabels(){

    livector1D labels;
    if(m_patch != nullptr){
        labels.reserve(m_patch->getNVertices());
        for(bitpit::Vertex & vert : m_patch->getVertices()){
            labels.push_back(vert.getId());
        }
    }
    return labels;
};

/*!
    \return final mirrored point cloud.
*/
MimmoSharedPointer<MimmoObject>
SpecularPoints::getMirroredPointCloud(){
    return m_patch;
}

/*!
 * Returns plane set up in the class, for mirroring
 * \returns plane coefficients
 */
darray4E
SpecularPoints::getPlane(){
    return  m_plane;
};

/*!
 * Returns which half-space intercepeted by the plane is interested by mirroring.
 * False represents the half-space where plane normal is directed, true the other one.
 * \returns insideout flag
 */
bool
SpecularPoints::isInsideOut(){
    return m_insideout;
}

/*!
 * Returns if even the points belonging to the symmetry
 * plane has to be mirrored (i.e. duplicated).
 * \returns true if points on symmetry plane has to be duplicated
 */
bool
SpecularPoints::isForced(){
    return m_force;
}


/*!
   It sets the target point cloud to be processed.
 * \param[in] targetpatch point cloud to be mirrored .
 */
void
SpecularPoints::setPointCloud(MimmoSharedPointer<MimmoObject> targetpatch){

    if(targetpatch == nullptr) return;
    if(targetpatch->getType() !=3) return;
    m_pc = targetpatch;
};


/*!
 * Set the original scalar data field attached to target point cloud
 * \param[in] scalardata scalar field pointer;
 */
void
SpecularPoints::setScalarData(dmpvector1D * scalardata){
    if(scalardata == nullptr) return;
    if(scalardata->getDataLocation() != MPVLocation::POINT)  return;
    m_scalar = scalardata;
};

/*!
 * Set the original vector data field attached to target point cloud
 * \param[in] vectordata vector field pointer;
 */
void
SpecularPoints::setVectorData(dmpvecarr3E * vectordata){
    if(vectordata == nullptr) return;
    if(vectordata->getDataLocation() != MPVLocation::POINT)  return;
    m_vector = vectordata;
};

/*!
 * Set Plane for mirroring cloud points. All points not belonging to plane will be mirrored
 * \param[in] plane coefficients a,b,c,d of plane in its implicit form a*x+b*y+c*z+d = 0
 */
void
SpecularPoints::setPlane(darray4E plane){
    m_plane = plane;
    m_implicit = true;
};

/*!
 * Set Plane for mirroring cloud points. All points not belonging to plane will be mirrored
 * \param[in] origin points belonging to plane
 * \param[in] normal plane normal
 *
 */
void
SpecularPoints::setPlane(darray3E origin, darray3E normal){

    double normx=norm2(normal);
    if(normx > std::numeric_limits<double>::min())  normal /= normx;

    m_plane[0] = normal[0];
    m_plane[1] = normal[1];
    m_plane[2] = normal[2];
    m_plane[3] = -1.0*dotProduct(origin, normal);

    m_implicit = true;

};

/*!
 * It sets origin of mirroring plane
 * \param[in] origin point belonging to plane
 *
 */
void
SpecularPoints::setOrigin(darray3E origin){
    m_origin = origin;
};

/*!
 * It sets normal of mirroring plane
 * \param[in] normal plane normal
 *
 */
void
SpecularPoints::setNormal(darray3E normal){
    m_normal = normal;
};

/*!
 * Returns which half-space intercepeted by the plane is interested by mirroring.
 * \param[in] flag false to select the half-space where plane normal is directed, true to select the other one.
 */
void
SpecularPoints::setInsideOut(bool flag){
    m_insideout = flag;
};

/*!
 * Set if even the points belonging to the symmetry
 * plane have to be mirrored (i.e. duplicated).
 * \param[in] flag true if the points on the symmetry plane have to be duplicated during mirroring.
 */
void
SpecularPoints::setForce(bool flag){
    m_force = flag;
}


/*!Execution command.Mirror the list of points linked, with data attached if any.
 * If a geometry is linked, project all resulting points on it.
 */
void
SpecularPoints::execute(){

    //check if there is a valid a point cloud
    if (m_pc == nullptr){
        (*m_log)<<"Error in "<<m_name << " : nullptr pointer to target point cloud"<<std::endl;
        throw std::runtime_error (m_name + " : nullptr pointer to target point cloud");
    }

    //clone it into the internal member m_cobj (mother class ProjPatchOnSurface)
    m_cobj = m_pc->clone();
    //initialize mirrored data;
    m_scalarMirrored.clear();
    m_vectorMirrored.clear();
    m_scalarMirrored.initialize(m_cobj, MPVLocation::POINT, 0.);
    m_scalarMirrored.setName("MirroredScalarData");
    m_vectorMirrored.initialize(m_cobj, MPVLocation::POINT, {{0.0,0.0,0.0}});
    m_vectorMirrored.setName("MirroredVectorData");


    /* Check the plane stuff and if an implicit definition is not present it has to be computed
     * by using origin and normal.
     */
    if (!m_implicit) setPlane(m_origin, m_normal);

    darray3E norm;
    double offsetPlane;
    for(int i=0; i<3; ++i)    norm[i] = m_plane[i];
    offsetPlane = m_plane[3];
    double normPlane = norm2(norm);

    //if normal is 0 there is nothing to do
    if(normPlane <= std::numeric_limits<double>::min()){
        (*m_log)<< "Error: "<<m_name <<" : no valid plane normal found"<<std::endl;
        m_patch = m_cobj;
        return;
    }

    //adjust the normal to be unitary in modulus
    norm /= normPlane;

    //see if need to project the cloud on a surface geometry: set bool project as active.
    bool project = false;
    int cellcount(0);
    if (getGeometry() != nullptr){
        cellcount = getGeometry()->getNCells();
#if MIMMO_ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &cellcount, 1, MPI_INT, MPI_SUM, m_communicator);
#endif
        project = (cellcount > 0);
    }


    //calculate a reasonable margin for defining the plane offset. This is
    // done to distinguish quickly points which lies on plane or not.
    double margin = 1.0E-12;
    if(project){
        double areaTot = 0.0;
        bitpit::SurfaceKernel * tri = static_cast<bitpit::SurfaceKernel * >(getGeometry()->getPatch());
        for(auto it = tri->internalBegin(); it != tri->internalEnd(); ++it ){
            areaTot += tri->evalCellArea(it.getId());
        }
#if MIMMO_ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &areaTot, 1, MPI_DOUBLE, MPI_SUM, m_communicator);
#endif
        if(areaTot > std::numeric_limits<double>::min()) {
            margin =  1.2*std::pow(areaTot/((double)cellcount),0.5);
        }
    }

    //take into account the right emiplane with insideout given by the User.
    double sig = (1.0  - 2.0*((int)m_insideout));

    //before mirroring, evaluate starting offsets for labeling new vertices to be added.

    livector1D localVertIds = m_cobj->getVerticesIds();
    int nPCvert = int(localVertIds.size());

    long offsetVertLabel(0);
    if(!localVertIds.empty())
        offsetVertLabel = *(std::max_element(localVertIds.begin(), localVertIds.end())) + 1;
    long offsetCellLabel(0);
    {
        livector1D localCellIds = m_cobj->getCellsIds();
        if(!localCellIds.empty())
            offsetCellLabel = *(std::max_element(localCellIds.begin(), localCellIds.end())) + 1;
    }


#if MIMMO_ENABLE_MPI
    std::vector<int> rankPCSize (m_nprocs, 0);
    rankPCSize[m_rank] = nPCvert;
    MPI_Allreduce(MPI_IN_PLACE, &offsetVertLabel, 1, MPI_LONG, MPI_MAX, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &offsetCellLabel, 1, MPI_LONG, MPI_MAX, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, rankPCSize.data(), m_nprocs, MPI_INT, MPI_MAX, m_communicator);

    for(int i=0; i<m_rank; ++i) {
        offsetVertLabel += long(rankPCSize[i]);
        offsetCellLabel += long(rankPCSize[i]);
    }
#endif

    //copy data of scalar and vector (if any) into mirrored twin structures
    //(that are already initialized to default 0 value) and reserve them
    //for new vertices insertion
    if(m_scalar){
        for(long id: localVertIds){
            if(m_scalar->exists(id)){
                m_scalarMirrored[id] = m_scalar->at(id);
            }
        }
        m_scalarMirrored.reserve(2*nPCvert);
    }

    if(m_vector){
        for(long id: localVertIds){
            if(m_vector->exists(id)){
                m_vectorMirrored[id] = m_vector->at(id);
            }
        }
        m_vectorMirrored.reserve(2*nPCvert);
    }

    m_cobj->getPatch()->reserveVertices(2*nPCvert);
    m_cobj->getPatch()->reserveCells(2*nPCvert);


    //FINALLY mirror the points, scalar and vector fields.
    darray3E vert;
    bitpit::ElementType type = bitpit::ElementType::VERTEX;
    for (long id: localVertIds){
        vert = m_cobj->getVertexCoords(id);
        double distance = sig*(dotProduct(norm, vert) + offsetPlane);

        if(distance > margin || m_force){

            m_cobj->addVertex( (vert - 2.0*distance*sig*norm), offsetVertLabel);
            m_cobj->addConnectedCell(std::vector<long>(1, offsetVertLabel), type, 0, offsetCellLabel, m_rank);

            m_scalarMirrored.insert(offsetVertLabel, m_scalarMirrored[id] );
            m_vectorMirrored.insert(offsetVertLabel, (m_vectorMirrored[id] -2.0*dotProduct(m_vectorMirrored[id], sig*norm)*sig*norm) );

            ++offsetVertLabel;
            ++offsetCellLabel;
        }
    }

    if(project){
        //this chunk will project the original cloud onto the target surface
        //m_cobj is filled up before.
        m_workingOnTarget = true;
        //final m_patch is m_cobj itself;
        projection();
        //for MPI version, this method already update the patch.
    }else{
        m_patch = m_cobj;

#if MIMMO_ENABLE_MPI

        if(m_patch->getPatch()->isPartitioned()){
            m_patch->updatePointGhostExchangeInfo();
        }
#endif

    }

    //squeeze stuff out
    m_patch->getPatch()->squeezeVertices();
    m_patch->getPatch()->squeezeCells();
    m_scalarMirrored.shrinkToFit();
    m_vectorMirrored.shrinkToFit();

    //evaluate kdtree if required by the class
    if(m_patch){
        if(m_buildKdTree)    m_patch->buildKdTree();
    }else{
        (*m_log)<<m_name << " : failed mirroring Point Cloud"<<std::endl;
    }
};

/*!
 * Clear all content of the class
 */
void
SpecularPoints::clear(){
    ProjPatchOnSurface::clear();
    m_scalar = nullptr;
    m_vector = nullptr;
    m_scalarMirrored.clear();
    m_vectorMirrored.clear();
    m_plane.fill(0.0);
    m_insideout = false;
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SpecularPoints::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    //start absorbing
    ProjPatchOnSurface::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Force")){
        std::string input = slotXML.get("Force");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setForce(value);
    }

    if(slotXML.hasOption("InsideOut")){
        std::string input = slotXML.get("InsideOut");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setInsideOut(value);
    }

    if(slotXML.hasSection("Plane")){
        const bitpit::Config::Section & planeXML = slotXML.getSection("Plane");

        std::string input1 = planeXML.get("Point");
        std::string input2 = planeXML.get("Normal");
        input1 = bitpit::utils::string::trim(input1);
        input2 = bitpit::utils::string::trim(input2);

        darray3E temp1 = {{0.0,0.0,0.0}};
        darray3E temp2 = {{0.0,0.0,0.0}};

        if(!input1.empty()){
            std::stringstream ss(input1);
            ss>>temp1[0]>>temp1[1]>>temp1[2];
        }
        if(!input2.empty()){
            std::stringstream ss(input2);
            ss>>temp2[0]>>temp2[1]>>temp2[2];
        }

        setPlane(temp1, temp2);
    };
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SpecularPoints::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    ProjPatchOnSurface::flushSectionXML(slotXML, name);

    int value = m_force;
    slotXML.set("Force", std::to_string(value));

    value = m_insideout;
    slotXML.set("InsideOut", std::to_string(value));


    {
        darray4E org = getPlane();
        darray3E normal;
        darray3E point = {{0.0,0.0,0.0}};
        int imax = -1;
        double dum = 0.0;
        for(int i=0; i<3; ++i)    {
            normal[i] =org[i];
            if(abs(normal[i]) > dum) {
                imax = i;
                dum = abs(normal[i]);
            }
        }
        if(imax != -1)    point[imax] = -1.0*org[3]/normal[imax];

        std::stringstream ss1, ss2;
        ss1<<std::scientific<<point[0]<<'\t'<<point[1]<<'\t'<<point[2];
        ss2<<std::scientific<<normal[0]<<'\t'<<normal[1]<<'\t'<<normal[2];

        bitpit::Config::Section & planeXML = slotXML.addSection("Plane");
        planeXML.set("Point",ss1.str());
        planeXML.set("Normal",ss2.str());
    }

};

/*!
 * Plot as optional results the mirrored list of points with the updated
 * data field associated to it
 */
void
SpecularPoints::plotOptionalResults(){
    write(m_patch, m_scalarMirrored, m_vectorMirrored);
};

}
