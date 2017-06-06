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

#include "ReconstructFields.hpp"
namespace mimmo{

/*!
 * Constructor
 */
ReconstructScalar::ReconstructScalar(){
    m_name = "mimmo.ReconstructScalar";
    m_overlapCriterium = OverlapMethod::MAX;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ReconstructScalar::ReconstructScalar(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.ReconstructScalar";
    m_overlapCriterium = OverlapMethod::MAX;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.ReconstructScalar"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Destructor
 */
ReconstructScalar::~ReconstructScalar(){
    clear();
}

/*!
 * Copy Constructor
 */
ReconstructScalar::ReconstructScalar(const ReconstructScalar & other):BaseManipulation(){
    *this = other;
}

/*!
 * Copy Operator
 */
ReconstructScalar & ReconstructScalar::operator=(const ReconstructScalar & other){
    *(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
    m_overlapCriterium = other.m_overlapCriterium;
    m_subpatch = other.m_subpatch;
    m_result = other.m_result;
    return *this;
}

/*!
 * Get actually used criterium for overlap regions of given fields
 * \return criterium for overlap regions
 */
OverlapMethod
ReconstructScalar::getOverlapCriteriumENUM(){
    return m_overlapCriterium;
}

/*!
 * Get actually used criterium for overlap regions of given fields
 * \return criterium for overlap regions
 */
int    ReconstructScalar::getOverlapCriterium(){
    return static_cast<int>(m_overlapCriterium);
}

/*!
 * Return data pointed for a given subpatch mesh
 * \param[in] patch    pointer to a subpatch
 * \return Data of scalar field associated to the patch, if any.
 */
dvector1D
ReconstructScalar::getData(MimmoObject * patch ){

    std::unordered_map<MimmoObject *, dvector1D *>::iterator it = m_subpatch.find(patch);

    dvector1D result;
    if(it == m_subpatch.end())    return result;

    return *(it->second);
}

/*!
 * Return number of fields data actually set in your class
 * \return number of fields
 */
int
ReconstructScalar::getNData(){

    return m_subpatch.size();
}

/*!
 * Return your result field
 * \return result fields
 */
dvector1D
ReconstructScalar::getResultField(){
    return(m_result);
}; 

/*!
 * Return actual computed scalar field (if any) for the geometry linked.
 * If no field is actually present, return null pointers;
 * \return     std::pair of pointers linking to actual geometry pointed by the class, and the computed deformation field on its vertices
 */
std::pair<MimmoObject * , dvector1D * >
ReconstructScalar::getResultFieldPair(){

    std::pair<MimmoObject *, dvector1D * > pairField;
    pairField.first = getGeometry();
    pairField.second = &m_result;
    return pairField;
};

/*!
 * Return list of sub-patch meshes actually stored in the class as a vector of copied pointers.
 * \return list of sub-patch meshes
 */
std::vector<MimmoObject    *>
ReconstructScalar::whichSubMeshes(){
    std::vector<MimmoObject    *> result(getNData());
    int counter=0;
    for (auto && pairInd : m_subpatch){
        result[counter] = pairInd.first;
        ++counter;
    }
    return result;
};

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::MAX. Enum overloading
 * \param[in] funct    Type of overlap method
 */
void
ReconstructScalar::setOverlapCriteriumENUM( OverlapMethod funct){
    setOverlapCriterium(static_cast<int>(funct));
};

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::MAX. Enum overloading
 * \param[in] funct    Type of overlap method
 */
void
ReconstructScalar::setOverlapCriterium( int funct){
    if(funct <1 ||funct > 4)    return;
    m_overlapCriterium = static_cast<OverlapMethod>(funct);
};


/*!
 * Insert in the list data field of a sub-patch, as typedef pField
 * (pointer to the sub-patch mesh, pointer to the sub-patch field)
 * \param[in] field    Sub-patch to be inserted
 */
void
ReconstructScalar::setData( pField  field){
    m_subpatch.insert(field);
};

/*!
 * Insert in the list data field of a sub-patch, as typedef pField
 * (pointer to the sub-patch mesh, pointer to the sub-patch field)
 * \param[in] fieldList    Vector of sub-patch to be inserted
 */
void
ReconstructScalar::setData(std::vector<pField>  fieldList){
    for(auto && data : fieldList){
        setData(data);
    }
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 * \param[in] patch Pointer to sub-patch to be removed.
 */
void
ReconstructScalar::removeData(MimmoObject * patch){
    std::unordered_map<MimmoObject *, dvector1D *>::iterator it = m_subpatch.find(patch);
    if(it != m_subpatch.end()){
        m_subpatch.erase(it);
    }
};

/*!
 * Remove all data present in the list
 */
void
ReconstructScalar::removeAllData(){
    m_subpatch.clear();
    m_result.clear();
};

/*!
 * Clear your class data and restore defaults settings
 */
void
ReconstructScalar::clear(){
    BaseManipulation::clear();
    removeAllData();
    m_overlapCriterium = OverlapMethod::MAX;
}

/*!
 * Plot data (resulting field data) on vtu unstructured grid file
 * \param[in]    dir        Output directory
 * \param[in]    name    Output filename
 * \param[in]    flag    Writing codex flag, false ascii, binary true
 * \param[in]    counter Counter identifying your output name
 */
void
ReconstructScalar::plotData(std::string dir, std::string name, bool flag, int counter){

    if(getGeometry() == NULL || getGeometry()->isEmpty())    return;
    dvecarr3E points = getGeometry()->getVertexCoords();
    ivector2D connectivity;
    bitpit::VTKElementType cellType = getGeometry()->desumeElement();

    if (getGeometry()->getType() != 3){
        connectivity = getGeometry()->getCompactConnectivity();
    }
    else{
        int np = points.size();
        connectivity.resize(np);
        for (int i=0; i<np; i++){
            connectivity[i].resize(1);
            connectivity[i][0] = i;
        }
    }
    bitpit::VTKUnstructuredGrid output(dir,name,cellType);
    output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
    output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
    output.setDimensions(connectivity.size(), points.size());

    output.addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_result);

    std::vector<long> ids(points.size());
    long ID;
    for (auto vertex : getGeometry()->getVertices()){
        ID = vertex.getId();
        ids[getGeometry()->getMapDataInv(ID)] = ID;
    }

    output.addData("ID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, ids);

    output.setCounter(counter);
    output.setCodex(bitpit::VTKFormat::APPENDED);
    if(!flag) output.setCodex(bitpit::VTKFormat::ASCII);

    output.write();
};

/*!
 * Execution command.
 * Reconstruct fields and save result in m_results member.
 */
void
ReconstructScalar::execute(){

    if(getGeometry() == NULL)    return;

    liimap & vMotherMap = getGeometry()->getMapDataInv();

    m_result.resize(getGeometry()->getPatch()->getVertexCount(), 0.0);

    if(m_subpatch.empty()) return;

    std::unordered_map<long, dvector1D> map = checkOverlapping();

    for(auto && obj : map){
        m_result[vMotherMap[obj.first]] = overlapFields(obj.second);
        obj.second.clear();
    }

}

/*!
 * Plot Optional results in execution, that is the reconstructed scalar field .
 */
void     ReconstructScalar::plotOptionalResults(){
    std::string dir = m_outputPlot;
    std::string name = m_name;
    plotData(dir, name, true, getClassCounter());
}

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 *\param[in] locField List of value of concurrent field. If value is unique, simply assigns it
 *\return    assigned value
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
double
ReconstructScalar::overlapFields(dvector1D & locField){
    int size  = locField.size();
    if(size < 1) return 0.0;

    double value = locField[0];
    if(size ==1 )return value;

    switch(m_overlapCriterium){
    case OverlapMethod::MAX :
        for(auto && loc : locField){
            value = std::fmax(loc,value);
        }
        break;
    case OverlapMethod::MIN :
        for(auto && loc : locField){
            value = std::fmin(loc,value);
        }
        break;
    case OverlapMethod::AVERAGE :
        value = 0.0;
        for(auto && loc : locField){
            value += loc/double(size);
        }
        break;
    case OverlapMethod::SUM :
        value = 0.0;
        for(auto && loc : locField){
            value += loc;
        }
        break;
    default : //never been reached
        break;
    }

    return value;
};

/*!
 * Check your sub-patch fields and reorder them in an unordered map carrying as
 * key the bitpit::PatchKernel ID of the mother mesh vertex and as value a vector
 * of doubles with all different field values concurring in that vertex
 *\return reordered map
 */
std::unordered_map<long, dvector1D>
ReconstructScalar::checkOverlapping(){

    std::unordered_map<long, dvector1D> result;
    int counter;

    for(auto && pairInd : m_subpatch){

        livector1D & vMap = pairInd.first->getMapData();
        counter = 0;
        dvector1D field = *(pairInd.second);
        for(int i=0; i<(int)field.size(); ++i){
            result[vMap[counter]].push_back(field[i]);
            ++counter;
        }
    }

    return(result);
};


/*! 
 * It builds the input/output ports of the object
 */
void
ReconstructScalar::buildPorts(){

    bool built = true;

    //input
    built = (built && createPortIn<MimmoObject *, ReconstructScalar>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    built = (built && createPortIn<std::pair<MimmoObject *, dvector1D *>,ReconstructScalar>(this, &mimmo::ReconstructScalar::setData, PortType::M_PAIRSCAFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECFLOAT_));
    built = (built && createPortIn<std::vector<std::pair<MimmoObject *, dvector1D *>>,ReconstructScalar>(this, &mimmo::ReconstructScalar::setData, PortType::M_VECPAIRSF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECFLOAT_));

    //output
    built = (built && createPortOut<dvector1D, ReconstructScalar>(this, &ReconstructScalar::getResultField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<MimmoObject *, ReconstructScalar>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<std::pair<MimmoObject *, dvector1D *>,ReconstructScalar>(this, &mimmo::ReconstructScalar::getResultFieldPair, PortType::M_PAIRSCAFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECFLOAT_));
    m_arePortsBuilt = built;
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ReconstructScalar::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    //start absorbing
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("OverlapCriterium")){
        std::string input = slotXML.get("OverlapCriterium");
        input = bitpit::utils::trim(input);
        int value = 1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::min(std::max(1, value),4);
        }
        setOverlapCriterium(value);
    }
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ReconstructScalar::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    
    int value = static_cast<int>(m_overlapCriterium);
    slotXML.set("OverlapCriterium", std::to_string(value));

};


}
