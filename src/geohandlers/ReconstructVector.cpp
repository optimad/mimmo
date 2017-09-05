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
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.ReconstructVector"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
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
ReconstructVector::ReconstructVector(const ReconstructVector & other):BaseManipulation(other){
    m_overlapCriterium = other.m_overlapCriterium;
    m_subpatch = other.m_subpatch;
    m_result = other.m_result;
    m_subresults = other.m_subresults;
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
 * Return number of fields data actually set in your class
 * \return number of fields
 */
int
ReconstructVector::getNData(){
    return m_subpatch.size();
}

/*!
 * Return your result field
 * \return result field
 */
dmpvecarr3E
ReconstructVector::getResultField(){
    return(m_result);
};

/*!
 * Return your result fields
 * \return result fields
 */
std::vector<dmpvecarr3E>
ReconstructVector::getResultFields(){
    return(m_subresults);
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
ReconstructVector::addData(dmpvecarr3E field){
    m_subpatch.push_back(field);
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 * \param[in] patch Pointer to sub-patch to be removed.
 */
void
ReconstructVector::removeData(MimmoObject * patch){
    std::vector<dmpvecarr3E>::iterator it = m_subpatch.begin();
    while(it != m_subpatch.end()){
        if (it->getGeometry() == patch){
            m_subpatch.erase(it);
        }
        else{
            ++it;
        }
    }
};

/*!
 * Remove all data present in the list
 */
void
ReconstructVector::removeAllData(){
    m_subpatch.clear();
    m_result.clear();
    m_subresults.clear();
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
    if(getGeometry() == NULL || getGeometry()->isEmpty()) return;

    liimap mapData;
    dvecarr3E points = getGeometry()->getVertexCoords(&mapData);
    ivector2D connectivity;
    bitpit::VTKElementType cellType = getGeometry()->desumeElement();

    if (getGeometry()->getType() != 3){
        connectivity = getGeometry()->getCompactConnectivity(mapData);
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

    dvecarr3E field = m_result.getDataAsVector();
    std::vector<long> ids = m_result.getIds();

    output.addData("vectorfield", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, field);
    output.addData("ID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, ids);

    output.setCounter(counter);
    output.setCodex(bitpit::VTKFormat::APPENDED);
    if(!flag) output.setCodex(bitpit::VTKFormat::ASCII);
    output.write();
};

/*!
 * Plot sub data (resulting field data) on vtu unstructured grid file
 * \param[in]    dir        Output directory
 * \param[in]    name    Output filename (the function will add SubPatch-i to this name)
 * \param[in]    i       index of the sub-patch
 * \param[in]    flag    Writing codex flag, false ascii, binary true
 * \param[in]    counter Counter identifying your output name
 */
void
ReconstructVector::plotSubData(std::string dir, std::string name, int i, bool flag, int counter){
    if(m_subresults[i].getGeometry() == NULL || m_subresults[i].getGeometry()->isEmpty()) return;

    name = name+"SubPatch"+to_string(i);

    liimap mapData;
    dvecarr3E points = m_subresults[i].getGeometry()->getVertexCoords(&mapData);
    ivector2D connectivity;
    bitpit::VTKElementType cellType = m_subresults[i].getGeometry()->desumeElement();

    if (m_subresults[i].getGeometry()->getType() != 3){
        connectivity = m_subresults[i].getGeometry()->getCompactConnectivity(mapData);
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

    dvecarr3E field = m_subresults[i].getDataAsVector();
    std::vector<long> ids = m_subresults[i].getIds();

    output.addData("vectorfield", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, field);
    output.addData("ID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, ids);

    output.setCounter(counter);
    output.setCodex(bitpit::VTKFormat::APPENDED);
    if(!flag) output.setCodex(bitpit::VTKFormat::ASCII);
    output.write();
};

/*!
 * Execution command.
 * Reconstruct fields and save result in results member.
 */
void
ReconstructVector::execute(){

    //Overlap fields
    m_result.clear();
    m_result.reserve(getGeometry()->getNVertex());
    m_subresults.clear();
    bitpit::PiercedVector<int> counter;
    for (int i=0; i<getNData(); i++){
        dmpvecarr3E* pv = &m_subpatch[i];
        long int ID;
        for (const auto & vertex : pv->getGeometry()->getVertices()){
            ID = vertex.getId();
            if (!m_result.exists(ID)){
                m_result.insert(ID, (*pv)[ID]);
                counter.insert(ID, 1);
            }
            else{
                overlapFields(ID, (*pv)[ID]);
                counter[ID]++;
            }
        }
    }
    if (m_overlapCriterium == OverlapMethod::AVERAGE){
        long int ID;
        auto itend = m_result.end();
        for (auto it=m_result.begin(); it!=itend; ++it){
            ID = it.getId();
            for (int j=0; j<3; j++)
                (*it)[j] = (*it)[j] / counter[ID];
        }
    }

    //Create subresults
    m_subresults.resize(getNData());
    for (int i=0; i<getNData(); i++){
        dmpvecarr3E* pv = &m_subpatch[i];
        m_subresults[i].setGeometry(pv->getGeometry());
//         m_subresults[i].setName(pv->getName());
        m_subresults[i].reserve(pv->getGeometry()->getNVertex());
        long int ID;
        for (const auto & vertex : pv->getGeometry()->getVertices()){
            ID = vertex.getId();
            m_subresults[i].insert(ID, m_result[ID]);
        }
    }

    //Update field on whole geometry
    if(getGeometry() != NULL){
        m_result.setGeometry(getGeometry());
//         if (m_subresults.size() != 0)
//             m_result.setName(m_subresults[0].getName());
        darray3E zero;
        zero.fill(0.0);
        long int ID;
        for (const auto & vertex : getGeometry()->getVertices()){
            ID = vertex.getId();
            if (!m_result.exists(ID)){
                m_result.insert(ID, zero);
            }
        }
    }
    else{
        m_result.clear();
    }
}

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 * \param[in] ID of the vertex to be checked.
 *\param[in] value Value of concurrent field. It updates the value in result field by using the input value of actual concurrent field.
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
void
ReconstructVector::overlapFields(long int ID, darray3E & locField){

    switch(m_overlapCriterium){
    case OverlapMethod::MAX :

    {
        double actual = norm2(m_result[ID]);
        double dummy = norm2(locField);

        if (actual < dummy)
            m_result[ID] = locField;
    }
    // TODO ??? WHAT'S IT ???
    //        value = {{0.0,0.0,0.0}};
    //        for(auto && loc : locField){
    //            dir += loc;
    //        }
    //
    //        dir /= norm2(dir);
    //
    //        match = 1.e-18;
    //        for(auto && loc : locField){
    //            double dummy = dotProduct(loc, dir);
    //            match = std::fmax(match,dummy);
    //        }
    //
    //        value = match*dir;
    break;

    case OverlapMethod::MIN :
    {
        double actual = norm2(m_result[ID]);
        double dummy = norm2(locField);

        if (actual > dummy)
            m_result[ID] = locField;
    }
    break;

    //        value = {{0.0,0.0,0.0}};
    //        for(auto && loc : locField){
    //            dir += loc;
    //        }
    //
    //        dir /= norm2(dir);
    //
    //        match = 1.e18;
    //        for(auto && loc : locField){
    //            double dummy = dotProduct(loc, dir);
    //            match = std::fmin(match,dummy);
    //        }
    //
    //        value = match*dir;
    //        break;

    case OverlapMethod::AVERAGE :

        m_result[ID] += locField;
        break;

    case OverlapMethod::SUM :

        m_result[ID] += locField;
        break;

    default : //never been reached
        break;
    }

};

/*! 
 * It builds the input/output ports of the object
 */
void
ReconstructVector::buildPorts(){

    bool built = true;

    //input
    built = (built && createPortIn<MimmoObject *, ReconstructVector>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortIn<dmpvecarr3E, ReconstructVector>(this, &mimmo::ReconstructVector::addData, PortType::M_VECTORFIELD, mimmo::pin::containerTAG::MPVECARR3, mimmo::pin::dataTAG::FLOAT));

    //output
    built = (built && createPortOut<dmpvecarr3E, ReconstructVector>(this, &ReconstructVector::getResultField, PortType::M_VECTORFIELD, mimmo::pin::containerTAG::MPVECARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<MimmoObject *, ReconstructVector>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<std::vector<dmpvecarr3E>, ReconstructVector>(this, &mimmo::ReconstructVector::getResultFields, PortType::M_VECVFIELDS, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MPVECARR3FLOAT));
    m_arePortsBuilt = built;
};

/*!
 * Plot Optional results in execution, that is the reconstructed vector field .
 */
void
ReconstructVector::plotOptionalResults(){
    std::string dir = m_outputPlot;
    std::string name = m_name;
    plotData(dir, name, true, getClassCounter());
    for (int i=0; i<getNData(); i++){
        plotSubData(dir, name, i, true, getClassCounter());
    }
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ReconstructVector::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    //start absorbing
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("OverlapCriterium")){
        std::string input = slotXML.get("OverlapCriterium");
        input = bitpit::utils::string::trim(input);
        int value = 4;
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
void ReconstructVector::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    int value = static_cast<int>(m_overlapCriterium);
    slotXML.set("OverlapCriterium", std::to_string(value));
};

}
