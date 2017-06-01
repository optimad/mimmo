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
#include "ScaleGeometry.hpp"

namespace mimmo{


/*!
 * Default constructor of ScaleGeometry
 */
ScaleGeometry::ScaleGeometry(darray3E scaling){
    m_scaling = scaling;
    m_origin.fill(0.0);
    m_meanP = false;
    m_name = "mimmo.ScaleGeometry";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ScaleGeometry::ScaleGeometry(const bitpit::Config::Section & rootXML){

    m_scaling.fill(1.0);
    m_origin.fill(0.0);
    m_meanP = false;
    m_name = "mimmo.ScaleGeometry";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.ScaleGeometry"){
        absorbSectionXML(rootXML);
    }else{
        (*m_log)<<"Warning in custom xml mimmo::ScaleGeometry constructor. No valid xml data found"<<std::endl;
    };
}

/*!Default destructor of ScaleGeometry
 */
ScaleGeometry::~ScaleGeometry(){};

/*!Copy constructor of ScaleGeometry.
 */
ScaleGeometry::ScaleGeometry(const ScaleGeometry & other):BaseManipulation(other){
    m_scaling = other.m_scaling;
    m_origin = other.m_origin;
    m_meanP = other.m_meanP;;
};

/*!Assignement operator of ScaleGeometry.
 */
ScaleGeometry & ScaleGeometry::operator=(const ScaleGeometry & other){
    *(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
    m_scaling = other.m_scaling;
    m_origin = other.m_origin;
    m_meanP = other.m_meanP;;
    return(*this);
};


/*! It builds the input/output ports of the object
 */
void
ScaleGeometry::buildPorts(){
    bool built = true;
    built = (built && createPortIn<darray3E, ScaleGeometry>(&m_origin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<darray3E, ScaleGeometry>(&m_scaling, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<dvector1D, ScaleGeometry>(this, &mimmo::ScaleGeometry::setFilter, PortType::M_FILTER, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<MimmoObject*, ScaleGeometry>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    built = (built && createPortOut<dvecarr3E, ScaleGeometry>(this, &mimmo::ScaleGeometry::getDisplacements, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<std::pair<MimmoObject*, dvecarr3E*> , ScaleGeometry>(this, &mimmo::ScaleGeometry::getDeformedField, PortType::M_PAIRVECFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
    built = (built && createPortOut<MimmoObject*, ScaleGeometry>(this, &BaseManipulation::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    m_arePortsBuilt = built;
};

/*!It sets the center point.
 * \param[in] origin origin for scaling transform
 */
void
ScaleGeometry::setOrigin(darray3E origin){
    m_origin = origin;
}

/*!It sets if the center point for scaling is the mean point of the geometry.
 * \param[in] meanP origin for scaling transform is the mean point?
 */
void
ScaleGeometry::setMeanPoint(bool meanP){
    m_meanP = meanP;
}

/*!It sets the scaling factors of each axis.
 * \param[in] scaling scaling factor values for x, y and z absolute axis.
 */
void
ScaleGeometry::setScaling(darray3E scaling){
    m_scaling = scaling;
}

/*!It sets the filter field to modulate the displacements of the vertices
 * of the target geometry.
 * \param[in] filter filter field defined on geometry vertices.
 */
void
ScaleGeometry::setFilter(dvector1D filter){
    m_filter = filter;
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * \return  deformation field
 */
dvecarr3E
ScaleGeometry::getDisplacements(){
    return m_displ;
};

/*!
 * Return actual computed deformation field (if any) for the geometry linked.
 * If no field is actually present, return null pointers;
 * \return     std::pair of pointers linking to actual geometry pointed by the class, and the computed deformation field on its vertices
 */
std::pair<MimmoObject * , dvecarr3E * >
ScaleGeometry::getDeformedField(){

    std::pair<MimmoObject *, dvecarr3E * > pairField;
    pairField.first = getGeometry();
    pairField.second = &m_displ;
    return pairField;
};

/*!Execution command. It perform the scaling by computing the displacements
 * of the points of the geometry. It applies a filter field eventually set as input.
 */
void
ScaleGeometry::execute(){

    if (getGeometry() == NULL) return;

    int nV = m_geometry->getNVertex();
    m_displ.resize(nV);
    m_filter.resize(nV, 1.0);


    long ID;
    int idx;
    liimap mapID = m_geometry->getMapDataInv();

    //computing centroid
    darray3E center = m_origin;
    if (m_meanP){
        center.fill(0.0);
        for (auto vertex : m_geometry->getVertices()){
            ID = vertex.getId();
            idx = mapID[ID];
            center += vertex.getCoords() / double(nV);
        }
    }
    for (auto vertex : m_geometry->getVertices()){
        ID = vertex.getId();
        idx = mapID[ID];

        darray3E coords = vertex.getCoords();
        m_displ[idx] = ( m_scaling*(coords - center) + center ) * m_filter[idx] - coords;
    }
};

/*!
 * Directly apply deformation field to target geometry.
 */
void
ScaleGeometry::apply(){

    if (getGeometry() == NULL) return;
    dvecarr3E vertex = getGeometry()->getVertexCoords();
    long nv = getGeometry()->getNVertex();
    nv = long(std::min(int(nv), int(m_displ.size())));
    livector1D & idmap = getGeometry()->getMapData();
    for (long i=0; i<nv; i++){
        vertex[i] += m_displ[i];
        getGeometry()->modifyVertex(vertex[i], idmap[i]);
    }

}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ScaleGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    if(slotXML.hasOption("Priority")){
        std::string input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

    if(slotXML.hasOption("MeanPoint")){
        std::string input = slotXML.get("MeanPoint");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setMeanPoint(value);
    };


    if(slotXML.hasOption("Origin")){
        std::string input = slotXML.get("Origin");
        input = bitpit::utils::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setOrigin(temp);
    }

    if(slotXML.hasOption("Scaling")){
        std::string input = slotXML.get("Scaling");
        input = bitpit::utils::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setScaling(temp);
    }

    if(slotXML.hasOption("Apply")){
        std::string input = slotXML.get("Apply");
        input = bitpit::utils::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setApply(value);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ScaleGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));

    {
        std::stringstream ss;
        ss<<std::scientific<<m_origin[0]<<'\t'<<m_origin[1]<<'\t'<<m_origin[2];
        slotXML.set("Origin", ss.str());
    }

    {
        bool value = m_meanP;
        std::string towrite = std::to_string(value);
        slotXML.set("MeanPoint", towrite);
    }

    {
        std::stringstream ss;
        ss<<std::scientific<<m_scaling[0]<<'\t'<<m_scaling[1]<<'\t'<<m_scaling[2];
        slotXML.set("Scaling", ss.str());
    }

    if(isApply()){
        slotXML.set("Apply", std::to_string(1));
    }

};

}
