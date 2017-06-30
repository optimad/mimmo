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
    input = bitpit::utils::string::trim(input);
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
ReconstructScalar::ReconstructScalar(const ReconstructScalar & other):BaseManipulation(other){
    m_overlapCriterium = other.m_overlapCriterium;
    m_subpatch = other.m_subpatch;
    m_result = other.m_result;
    m_subresults = other.m_subresults;
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
 * Return number of fields data actually set in your class
 * \return number of fields
 */
int
ReconstructScalar::getNData(){
    return m_subpatch.size();
}

/*!
 * Return your result field
 * \return result field
 */
dmpvector1D
ReconstructScalar::getResultField(){
    return(m_result);
}; 

/*!
 * Return your result fields
 * \return result fields
 */
std::vector<dmpvector1D>
ReconstructScalar::getResultFields(){
    return(m_subresults);
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
ReconstructScalar::addData( dmpvector1D  field){
    m_subpatch.push_back(field);
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 * \param[in] patch Pointer to sub-patch to be removed.
 */
void
ReconstructScalar::removeData(MimmoObject * patch){
    std::vector<dmpvector1D>::iterator it = m_subpatch.begin();
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
ReconstructScalar::removeAllData(){
    m_subpatch.clear();
    m_result.clear();
    m_subresults.clear();
};

/*!
 * Clear your class data and restore defaults settings
 */
void
ReconstructScalar::clear(){
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
ReconstructScalar::plotData(std::string dir, std::string name, bool flag, int counter){

    if(getGeometry() == NULL || getGeometry()->isEmpty())    return;
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
    output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
    output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
    output.setDimensions(connectivity.size(), points.size());


    dvector1D field(points.size());
    std::vector<long> ids(points.size());
    long ID;
    for (const auto vertex : getGeometry()->getVertices()){
        ID = vertex.getId();
        ids[mapData[ID]] = ID;
        field[mapData[ID]] = m_result[ID];
    }

    output.addData("scalarfield", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, field);
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
ReconstructScalar::plotSubData(std::string dir, std::string name, int i, bool flag, int counter){
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

    dvector1D field(points.size());
    std::vector<long> ids(points.size());
    long ID;
    for (const auto vertex : m_subresults[i].getGeometry()->getVertices()){
        ID = vertex.getId();
        ids[mapData[ID]] = ID;
        field[mapData[ID]] = m_subresults[i][ID];
    }

    output.addData("scalarfield", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, field);
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


    //Overlap fields
    m_result.clear();
    m_subresults.clear();
    bitpit::PiercedVector<int> counter;
    for (int i=0; i<getNData(); i++){
        dmpvector1D* pv = &m_subpatch[i];
        long int ID;
        for (const auto vertex : pv->getGeometry()->getVertices()){
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
        MimmoPiercedVector<double>::iterator it;
        MimmoPiercedVector<double>::iterator itend = m_result.end();
        for (it=m_result.begin(); it!=itend; ++it){
            ID = it.getId();
            (*it) = (*it) / counter[ID];
        }
    }

    //Create subresults
    m_subresults.resize(getNData());
    for (int i=0; i<getNData(); i++){
        dmpvector1D* pv = &m_subpatch[i];
        m_subresults[i].setGeometry(pv->getGeometry());
        m_subresults[i].setName(pv->getName());
        long int ID;
        for (const auto vertex : pv->getGeometry()->getVertices()){
            ID = vertex.getId();
            m_subresults[i].insert(ID, m_result[ID]);
        }
    }

    //Update field on whole geometry
    if(getGeometry() != NULL){
        m_result.setGeometry(getGeometry());
        if (m_subresults.size() != 0)
            m_result.setName(m_subresults[0].getName());
        double zero = 0.0;
        long int ID;
        for (const auto vertex : getGeometry()->getVertices()){
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
 * Plot Optional results in execution, that is the reconstructed scalar field .
 */
void     ReconstructScalar::plotOptionalResults(){
    std::string dir = m_outputPlot;
    std::string name = m_name;
    plotData(dir, name, true, getClassCounter());
    for (int i=0; i<getNData(); i++){
        plotSubData(dir, name, i, true, getClassCounter());
    }
}

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 *\param[in] locField List of value of concurrent field. If value is unique, simply assigns it
 *\return    assigned value
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
void
ReconstructScalar::overlapFields(long int ID, double & locField){

    switch(m_overlapCriterium){
    case OverlapMethod::MAX :
    {
        double actual = m_result[ID];
        double dummy = locField;

        if (actual < dummy)
            m_result[ID] = locField;
    }
    break;
    case OverlapMethod::MIN :
    {
        double actual = m_result[ID];
        double dummy = locField;

        if (actual > dummy)
            m_result[ID] = locField;
    }
    break;
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
ReconstructScalar::buildPorts(){

    bool built = true;

    //input
    built = (built && createPortIn<MimmoObject *, ReconstructScalar>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    built = (built && createPortIn<dmpvector1D, ReconstructScalar>(this, &mimmo::ReconstructScalar::addData, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::MPVECTOR, mimmo::pin::dataTAG::FLOAT));

    //output
    built = (built && createPortOut<dmpvector1D, ReconstructScalar>(this, &ReconstructScalar::getResultField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::MPVECTOR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<MimmoObject *, ReconstructScalar>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<std::vector<dmpvector1D>, ReconstructScalar>(this, &mimmo::ReconstructScalar::getResultFields, PortType::M_VECSFIELDS, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MPVECFLOAT));
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
        input = bitpit::utils::string::trim(input);
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
