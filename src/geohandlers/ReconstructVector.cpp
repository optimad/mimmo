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
ReconstructVector::ReconstructVector(MPVLocation loc){
    m_name = "mimmo.ReconstructVector";
    m_overlapCriterium = OverlapMethod::SUM;
    m_loc = loc;
    if(m_loc == MPVLocation::UNDEFINED || m_loc == MPVLocation::INTERFACE) m_loc= MPVLocation::POINT;
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
    std::string fallback_loc  = "-1";
    std::string input_loc = rootXML.get("DataLocation", fallback_loc);
    
    input = bitpit::utils::string::trim(input);
    input_loc = bitpit::utils::string::trim(input_loc);
    
    int loc = std::stoi(input_loc);
    if(loc > 0 && loc < 3){
        m_loc  =static_cast<MPVLocation>(loc);
    }else{
        m_loc = MPVLocation::POINT;
    }

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
    m_loc = other.m_loc;
    m_overlapCriterium = other.m_overlapCriterium;
    m_subpatch = other.m_subpatch;
    m_result = other.m_result;
    m_subresults = other.m_subresults;
}

/*!
 * Swap function of ReconstructVector
 * \param[in] x object to be swapped.
 */
void ReconstructVector::swap(ReconstructVector & x ) noexcept
{
    std::swap(m_loc, x.m_loc);
    std::swap(m_overlapCriterium,x.m_overlapCriterium);
    std::swap(m_subpatch, x.m_subpatch);
    std::swap(m_subresults, x.m_subresults);
    //std::swap(m_result, x.m_result);
    m_result.swap(x.m_result);
    BaseManipulation::swap(x);
};

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
 * Insert in the list data field of a sub-patch.
 * \param[in] field    Sub-patch to be inserted
 */
void
ReconstructVector::addData(dmpvecarr3E field){
    if(field.getGeometry()== NULL && field.getDataLocation() != m_loc) return;
    if(field.getGeometry()->getType()==3 && m_loc==MPVLocation::CELL){
        (*m_log)<<"warning in "<<m_name<<" : trying to add field referred to a Point Cloud, while Class has Data Location referred to CELLS. Do Nothing."<<std::endl;
        return;
    }
    m_subpatch.push_back(field);
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 * \param[in] patch Pointer to sub-patch to be removed.
 */
void
ReconstructVector::removeData(MimmoObject * patch){
    //in m_subpatch remove progressively all pierced vector elements in last position
    //which links towards a geometry of type patch.
    while(m_subpatch.back().getGeometry() == patch){
        m_subpatch.pop_back();
    }

    //start searching from the begin on element linking patch
    std::vector<dmpvecarr3E>::iterator it = m_subpatch.begin();
    
    while(it != m_subpatch.end()){
        if (it->getGeometry() == patch){
            *it = m_subpatch.back();
            m_subpatch.pop_back();
            //be always sure the last element is clean and does not link
            // the patch
            while(m_subpatch.back().getGeometry() == patch){
                m_subpatch.pop_back();
            }
        }
        ++it; //forward the iterator.
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
    m_overlapCriterium = OverlapMethod::SUM;
}

/*!
 * Plot data (resulting field data) on vtu unstructured grid file
 * \param[in]    dir        Output directory
 * \param[in]    name    Output filename
 * \param[in]    flag    Writing codex flag, false ascii, binary true
 */
void
ReconstructVector::plotData(std::string dir, std::string name, bool flag){
    
    if(getGeometry() == NULL) return;
    if(getGeometry()->isEmpty())    return;
    if(!m_result.completeMissingData({{0.0,0.0,0.0}}))   return;

    bitpit::VTKLocation loc = bitpit::VTKLocation::POINT;
    if(m_loc == MPVLocation::CELL){
        loc = bitpit::VTKLocation::CELL;
    }

    dvecarr3E field = m_result.getDataAsVector();

    if (getGeometry()->getType() != 3){
        getGeometry()->getPatch()->getVTK().addData("vectorfield", bitpit::VTKFieldType::VECTOR, loc, field);
        getGeometry()->getPatch()->getVTK().setDataCodex(bitpit::VTKFormat::APPENDED);
        if(!flag) getGeometry()->getPatch()->getVTK().setDataCodex(bitpit::VTKFormat::ASCII);
        getGeometry()->getPatch()->write(dir+"/"+name);
        getGeometry()->getPatch()->getVTK().removeData("vectorfield");
    }else{
        liimap mapData;
        dvecarr3E points = getGeometry()->getVertexCoords(&mapData);
        ivector2D connectivity;
        int np = points.size();
        connectivity.resize(np);
        for (int i=0; i<np; i++){
            connectivity[i].resize(1);
            connectivity[i][0] = i;
        }
        bitpit::VTKUnstructuredGrid output(dir,name,bitpit::VTKElementType::VERTEX);
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
        output.setDimensions(connectivity.size(), points.size());
        output.addData("vectorfield", bitpit::VTKFieldType::VECTOR, loc, field);
        output.setDataCodex(bitpit::VTKFormat::APPENDED);
        if(!flag) output.setDataCodex(bitpit::VTKFormat::ASCII);
        output.write();
    }
};

/*!
 * Plot sub data (resulting field data) on vtu unstructured grid file
 * \param[in]    dir        Output directory
 * \param[in]    name    Output filename (the function will add SubPatch-i to this name)
 * \param[in]    i       index of the sub-patch
 * \param[in]    flag    Writing codex flag, false ascii, binary true
 */
void
ReconstructVector::plotSubData(std::string dir, std::string name, int i, bool flag){
    if(m_subresults[i].getGeometry() == NULL) return;
    if(m_subresults[i].getGeometry()->isEmpty()) return;
    
    std::string nameX = name+"SubPatch"+std::to_string(i);

    bitpit::VTKLocation loc = bitpit::VTKLocation::POINT;
    if(m_loc == MPVLocation::CELL){
        loc = bitpit::VTKLocation::CELL;
    }

    //check size of field and adjust missing values to zero for writing purposes only.
    dmpvecarr3E field_supp = m_subresults[i];
    if(!field_supp.completeMissingData({{0.0,0.0,0.0}}))    return;
    dvecarr3E field = field_supp.getDataAsVector();

    if (m_subresults[i].getGeometry()->getType() != 3){
        m_subresults[i].getGeometry()->getPatch()->getVTK().addData("vectorfield", bitpit::VTKFieldType::VECTOR,loc, field);
        m_subresults[i].getGeometry()->getPatch()->getVTK().setDataCodex(bitpit::VTKFormat::APPENDED);
        if(!flag) m_subresults[i].getGeometry()->getPatch()->getVTK().setDataCodex(bitpit::VTKFormat::ASCII);
        m_subresults[i].getGeometry()->getPatch()->write(dir+"/"+nameX);
        m_subresults[i].getGeometry()->getPatch()->getVTK().removeData("vectorfield");
    }else{
        liimap mapData;
        dvecarr3E points = m_subresults[i].getGeometry()->getVertexCoords(&mapData);
        ivector2D connectivity;
        int np = points.size();
        connectivity.resize(np);
        for (int i=0; i<np; i++){
            connectivity[i].resize(1);
            connectivity[i][0] = i;
        }
        bitpit::VTKUnstructuredGrid output(dir,nameX,bitpit::VTKElementType::VERTEX);
        output.setGeomData(bitpit::VTKUnstructuredField::POINTS, points);
        output.setGeomData(bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
        output.setDimensions(connectivity.size(), points.size());
        output.addData("vectorfield", bitpit::VTKFieldType::VECTOR,loc, field);
        output.setDataCodex(bitpit::VTKFormat::APPENDED);
        if(!flag) output.setDataCodex(bitpit::VTKFormat::ASCII);
        output.write();
    }
};

/*!
 * Execution command.
 * Reconstruct fields and save result in results member.
 */
void
ReconstructVector::execute(){

    if(getGeometry() == NULL){
        throw std::runtime_error(m_name + "NULL pointer to linked geometry found");
    } 
    
    if(getGeometry()->isEmpty()){
        throw std::runtime_error(m_name + "empty linked geometry found");
    } 
    
    //Overlap fields
    m_result.clear();
    m_result.setGeometry(getGeometry());
    m_result.setDataLocation(m_loc);
    
    m_subresults.clear();
    
    std::unordered_set<long> idsTarget;
    {
        livector1D ids = idsGeoDataLocation(getGeometry());
        idsTarget.insert(ids.begin(), ids.end());
    }
    
    bitpit::PiercedVector<int> counter;
    for (int i=0; i<getNData(); i++){
        dmpvecarr3E* pv = &m_subpatch[i];
        livector1D ids = idsGeoDataLocation(pv->getGeometry());
        for (auto ID: ids){
            if(idsTarget.count(ID)==0) continue;
            
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
    
    if (m_result.isEmpty()){
        (*m_log)<<"Warning in "<<m_name<<". Resulting reconstructed field is empty.This is could be caused by unrelated fields linked geometry and target geometry"<<std::endl;
        throw std::runtime_error(m_name + "empty field reconstructed in class execution.");
    }
    
    if (m_overlapCriterium == OverlapMethod::AVERAGE){
        long int ID;
        auto itend = m_result.end();
        for (auto it=m_result.begin(); it!=itend; ++it){
            ID = it.getId();
            for (int j=0; j<3; j++)
                (*it)[j] = (*it)[j] / double(counter[ID]);
        }
    }

    //Update field on whole geometry
    darray3E zero = {{0.0,0.0,0.0}};
    m_result.completeMissingData(zero);
    
    //Create subresults
    m_subresults.resize(getNData());
    for (int i=0; i<getNData(); i++){
        dmpvecarr3E* pv = &m_subpatch[i];
        m_subresults[i].setGeometry(pv->getGeometry());
        m_subresults[i].setDataLocation(m_loc);
        livector1D ids = idsGeoDataLocation(pv->getGeometry());
        for (auto ID : ids){
            if(idsTarget.count(ID)>0){
                m_subresults[i].insert(ID, m_result[ID]);
            }
        }
    }
}

/*!
 * Plot Optional results in execution, that is the reconstructed vector field .
 */
void
ReconstructVector::plotOptionalResults(){
    std::string dir = m_outputPlot;
    std::string name = m_name;
    plotData(dir, name, true);
    for (int i=0; i<getNData(); i++){
        plotSubData(dir, name, i, true);
    }
}

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 * \param[in] ID of the vertex to be checked.
 *\param[in] locField Value of concurrent field. It updates the value in result field by using the input value of actual concurrent field.
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
    break;

    case OverlapMethod::MIN :
    {
        double actual = norm2(m_result[ID]);
        double dummy = norm2(locField);

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
ReconstructVector::buildPorts(){

    bool built = true;

    //input
    built = (built && createPortIn<MimmoObject *, ReconstructVector>(&m_geometry, M_GEOM, true,1));
    built = (built && createPortIn<dmpvecarr3E, ReconstructVector>(this, &mimmo::ReconstructVector::addData, M_VECTORFIELD));

    //output
    built = (built && createPortOut<dmpvecarr3E, ReconstructVector>(this, &ReconstructVector::getResultField, M_VECTORFIELD));
    built = (built && createPortOut<MimmoObject *, ReconstructVector>(&m_geometry, M_GEOM));
    built = (built && createPortOut<std::vector<dmpvecarr3E>, ReconstructVector>(this, &mimmo::ReconstructVector::getResultFields, M_VECVFIELDS));
    m_arePortsBuilt = built;
};




/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ReconstructVector::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    if(slotXML.hasOption("DataLocation")){
        std::string input = slotXML.get("DataLocation");
        input = bitpit::utils::string::trim(input);
        int temp = -1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        if(int(m_loc) != temp){
            (*m_log)<<"Error absorbing DataLocation in "<<m_name<<". Class and read locations mismatch"<<std::endl;
            if (temp == 0) (*m_log)<<"XML DataLocation in "<<m_name<<" is set to 0-UNDEFINED"<<std::endl;
            if (temp == 3) (*m_log)<<"XML DataLocation in "<<m_name<<" is set to 3-INTERFACE, not supported for now."<<std::endl;            
            throw std::runtime_error (m_name + " : xml absorbing failed.");
        }
    }
    
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

    slotXML.set("DataLocation", std::to_string(int(m_loc)));
    int value = static_cast<int>(m_overlapCriterium);
    slotXML.set("OverlapCriterium", std::to_string(value));
};

/*!
 * Given a reference geometry, return list of ids relative to geometry vertices or cells
 * according to data location parameter m_loc of the class
 *\param[in] geo valid pointer to a MimmoObject geometry
 *\return list of ids relative to vertices or cells according to class m_loc.
 */
livector1D ReconstructVector::idsGeoDataLocation(MimmoObject * geo){
    
    if (m_loc == MPVLocation::POINT) return geo->getVertices().getIds();
    if (m_loc == MPVLocation::CELL)  return geo->getCells().getIds();
    return livector1D(0);
}


}
