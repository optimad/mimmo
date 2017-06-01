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
ReconstructVector::ReconstructVector(){
    m_name = "mimmo.ReconstructVector";
    m_overlapCriterium = OverlapMethod::SUM;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ReconstructVector::ReconstructVector(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.ReconstructVector";
    m_overlapCriterium = OverlapMethod::MAX;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.ReconstructVector"){
        absorbSectionXML(rootXML);
    }else{
        (*m_log)<<"Warning in custom xml mimmo::ReconstructVector constructor. No valid xml data found"<<std::endl;
    };
}

/*!
 * Destructor
 */
ReconstructVector::~ReconstructVector(){
    clear();
}

/*!
 * Copy Constructor
 */
ReconstructVector::ReconstructVector(const ReconstructVector & other):BaseManipulation(){
    *this = other;
}

/*!
 * Copy Operator
 */
ReconstructVector & ReconstructVector::operator=(const ReconstructVector & other){
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
ReconstructVector::getOverlapCriteriumENUM(){
    return m_overlapCriterium;
}

/*!
 * Get actually used criterium for overlap regions of given fields
 * \return criterium for overlap regions
 */
int
ReconstructVector::getOverlapCriterium(){
    return static_cast<int>(m_overlapCriterium);
}

/*!
 * Return data pointed for a given sub-patch mesh
 * \param[in] patch    Pointer to a sub-patch
 * \return data of scalar field associated to the patch, if any. 
 */
dvecarr3E
ReconstructVector::getData(MimmoObject * patch ){

    std::unordered_map<MimmoObject *, dvecarr3E *>::iterator it = m_subpatch.find(patch);

    dvecarr3E result;
    if(it == m_subpatch.end())    return result;

    return *(it->second);
}

/*!
 * Return number of fields data actually set in your class
 * \return number of fields
 */
int
ReconstructVector::getNData(){

    return m_subpatch.size();
}

/*!
 * Return your result field
 * \return result fields
 */
dvecarr3E
ReconstructVector::getResultField(){
    return(m_result);
}; 

/*!
 * Return actual computed scalar field (if any) for the geometry linked.
 * If no field is actually present, return null pointers;
 * \return     std::pair of pointers linking to actual geometry pointed by the class, and the computed deformation field on its vertices
 */
std::pair<MimmoObject * , dvecarr3E * >
ReconstructVector::getResultFieldPair(){

    std::pair<MimmoObject *, dvecarr3E * > pairField;
    pairField.first = getGeometry();
    pairField.second = &m_result;
    return pairField;
};

/*!
 * Return list of sub-patch meshes actually stored in the class as a vector of copied pointers.
 * \return list of sub-patch meshes
 */
std::vector<MimmoObject    *>
ReconstructVector::whichSubMeshes(){
    std::vector<MimmoObject    *> result(getNData());
    int counter=0;
    for (auto && pairInd : m_subpatch){
        result[counter] = pairInd.first;
        ++counter;
    }
    return result;
};

/*!
 * Set overlap criterium for multi-fields reconstruction. See OverlapMethod enum
 * Class default is OverlapMethod::MAX. Enum overloading
 * \param[in] funct    Type of overlap method
 */
void
ReconstructVector::setOverlapCriteriumENUM( OverlapMethod funct){
    setOverlapCriterium(static_cast<int>(funct));
};


/*!
 * Set overlap criterium for multi-fields reconstruction. See OverlapMethod enum
 * Class default is OverlapMethod::MAX. Enum overloading
 * \param[in] funct    Type of overlap method
 */
void
ReconstructVector::setOverlapCriterium( int funct){
    if(funct<1 ||funct>4)    return;
    m_overlapCriterium = static_cast<OverlapMethod>(funct);
};

/*!
 * Insert in the list data field of a sub-patch, as typedef pVector
 * (pointer to the sub-patch mesh, pointer to the sub-patch field)
 * \param[in] vfield    Sub-patch to be inserted
 */
void
ReconstructVector::setData( pVector  vfield){
    m_subpatch.insert(vfield);
};

/*!
 * Insert in the list data field of a sub-patch, as typedef pField
 * (pointer to the sub-patch mesh, pointer to the sub-patch field)
 * \param[in] vfieldList    Vector of sub-patch to be inserted
 */
void
ReconstructVector::setData(std::vector<pVector>  vfieldList){
    for(auto && data : vfieldList){
        setData(data);
    }
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 * \param[in] patch Pointer to sub-patch to be removed.
 */
void
ReconstructVector::removeData(MimmoObject * patch){
    std::unordered_map<MimmoObject *, dvecarr3E *>::iterator it = m_subpatch.find(patch);
    if(it != m_subpatch.end()){
        m_subpatch.erase(it);
    }
};

/*!
 * Remove all data present in the list
 */
void
ReconstructVector::removeAllData(){
    m_subpatch.clear();
    m_result.clear();
};

/*!
 * Clear your class data and restore defaults settings
 */
void
ReconstructVector::clear(){
    BaseManipulation::clear();
    removeAllData();
    m_overlapCriterium = OverlapMethod::AVERAGE;
}

/*!
 * Plot data (resulting field data) on vtu unstructured grid file
 * \param[in]    dir        Output directory
 * \param[in]    name    Output filename
 * \param[in]    flag    Writing codex flag, false ascii, binary true
 * \param[in]    counter Counter identifying your output name
 */
void
ReconstructVector::plotData(std::string dir, std::string name, bool flag, int counter){
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
    output.setGeomData(bitpit::VTKUnstructuredField::POINTS, points);
    output.setGeomData(bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
    output.setDimensions(connectivity.size(), points.size());

    output.addData("vectorfield", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, m_result);

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
ReconstructVector::execute(){
    if(getGeometry() == NULL)    return;

    liimap & vMotherMap = getGeometry()->getMapDataInv();

    m_result.resize(getGeometry()->getPatch()->getVertexCount(), {{0.0,0.0,0.0}});
    if(m_subpatch.empty())    return;
    std::unordered_map<long, dvecarr3E > map = checkOverlapping();

    for(auto && obj : map){
        m_result[vMotherMap[obj.first]] = overlapFields(obj.second);
        obj.second.clear();
    }

}

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 *\param[in] locField List of value of concurrent field. If value is unique, simply assigns it
 *\return    assigned values
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
darray3E
ReconstructVector::overlapFields(dvecarr3E & locField){
    int size  = locField.size();
    if(size < 1) return {{0.0,0.0,0.0}};

    if(size ==1 )return locField[0];
    darray3E value;
    darray3E dir; dir.fill(0.0);
    double match;

    switch(m_overlapCriterium){
    case OverlapMethod::MAX :
        value = {{0.0,0.0,0.0}};
        for(auto && loc : locField){
            dir += loc;
        }

        dir /= norm2(dir);

        match = 1.e-18;
        for(auto && loc : locField){
            double dummy = dotProduct(loc, dir);
            match = std::fmax(match,dummy);
        }

        value = match*dir;
        break;

    case OverlapMethod::MIN :

        value = {{0.0,0.0,0.0}};
        for(auto && loc : locField){
            dir += loc;
        }

        dir /= norm2(dir);

        match = 1.e18;
        for(auto && loc : locField){
            double dummy = dotProduct(loc, dir);
            match = std::fmin(match,dummy);
        }

        value = match*dir;
        break;

    case OverlapMethod::AVERAGE :
        value = {{0.0,0.0,0.0}};
        for(auto && loc : locField){
            value += loc/double(size);
        }
        break;

    case OverlapMethod::SUM :
        value = {{0.0,0.0,0.0}};
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
 * of double arrays with all different vector field values concurring in that vertex
 *\return reordered map
 */
std::unordered_map<long, dvecarr3E>
ReconstructVector::checkOverlapping(){

    std::unordered_map<long, dvecarr3E> result;
    int counter;

    for(auto && pairInd : m_subpatch){

        livector1D & vMap = pairInd.first->getMapData();
        counter = 0;
        dvecarr3E field = *(pairInd.second);
        for(auto & val : field){
            result[vMap[counter]].push_back(val);
            ++counter;
        }
    }

    return(result);
};


/*! 
 * It builds the input/output ports of the object
 */
void
ReconstructVector::buildPorts(){

    bool built = true;

    //input
    built = (built && createPortIn<MimmoObject *, ReconstructVector>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    built = (built && createPortIn<std::pair<MimmoObject *, dvecarr3E *>,ReconstructVector>(this, &mimmo::ReconstructVector::setData, PortType::M_PAIRVECFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_, true, 1));
    built = (built && createPortIn<std::vector<std::pair<MimmoObject *, dvecarr3E *> >,ReconstructVector>(this, &mimmo::ReconstructVector::setData, PortType::M_VECPAIRVF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECARR3FLOAT_, true, 1));

    //output
    built = (built && createPortOut<dvecarr3E, ReconstructVector>(this, &ReconstructVector::getResultField, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<MimmoObject *, ReconstructVector>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<std::pair<MimmoObject *, dvecarr3E *>,ReconstructVector>(this, &mimmo::ReconstructVector::getResultFieldPair, PortType::M_PAIRVECFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
    m_arePortsBuilt = built;
};

/*!
 * Plot Optional results in execution, that is the reconstructed vector field .
 */
void
ReconstructVector::plotOptionalResults(){
    std::string dir = m_outputPlot ;
    std::string name = m_name;
    plotData(dir, name, true, getClassCounter());
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ReconstructVector::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    //start absorbing
    if(slotXML.hasOption("Priority")){
        std::string input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

    if(slotXML.hasOption("OverlapCriterium")){
        std::string input = slotXML.get("OverlapCriterium");
        input = bitpit::utils::trim(input);
        int value = 4;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::min(std::max(1, value),4);
        }
        setOverlapCriterium(value);
    }


    if(slotXML.hasOption("PlotInExecution")){
        std::string input = slotXML.get("PlotInExecution");
        input = bitpit::utils::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setPlotInExecution(value);
    }

    if(slotXML.hasOption("OutputPlot")){
        std::string input = slotXML.get("OutputPlot");
        input = bitpit::utils::trim(input);
        std::string temp = ".";
        if(!input.empty())    setOutputPlot(input);
        else                  setOutputPlot(temp);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ReconstructVector::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));

    int value = static_cast<int>(m_overlapCriterium);
    slotXML.set("OverlapCriterium", std::to_string(value));

    if(isPlotInExecution()){
        slotXML.set("PlotInExecution", std::to_string(1));
    }

    if(m_outputPlot != "."){
        slotXML.set("OutputPlot", m_outputPlot);
    }

};

}
