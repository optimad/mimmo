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
#include "OverlappingFields.hpp"
namespace mimmo{

/*!
 * Constructor
 */
OverlapVectorFields::OverlapVectorFields(){
    m_name = "mimmo.OverlapVectorFields";
    m_overlapCriterium = OverlapMethod::SUM;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
OverlapVectorFields::OverlapVectorFields(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.OverlapVectorFields";
    m_overlapCriterium = OverlapMethod::SUM;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.OverlapVectorFields"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Destructor
 */
OverlapVectorFields::~OverlapVectorFields(){
    clear();
}

/*!
 * Copy Constructor
 */
OverlapVectorFields::OverlapVectorFields(const OverlapVectorFields & other):BaseManipulation(){
    m_overlapCriterium = other.m_overlapCriterium;
    m_originals = other.m_originals;
}

/*!
 * Get actually used criterium for overlap regions of given fields
 * \return criterium for overlap regions
 */
OverlapMethod
OverlapVectorFields::getOverlapCriteriumENUM(){
    return m_overlapCriterium;
}

/*!
 * Get actually used criterium for overlap regions of given fields
 * \return criterium for overlap regions
 */
int
OverlapVectorFields::getOverlapCriterium(){
    return static_cast<int>(m_overlapCriterium);
}

/*!
 * Return overlapped data pointed for a given subpatch mesh
 * \param[in] patch    Pointer to a geometry subpatch
 * \return data of processed vector field associated to the patch, if any. 
 */
dvecarr3E
OverlapVectorFields::getResultData(MimmoObject * patch ){
    dvecarr3E result;
    if(m_results.count(patch) == 0) return result;
    return m_results[patch];
}

/*!
 * Return effective number of fields that would count after overlapping.
 * \return number of fields
 */
int
OverlapVectorFields::getNEffectiveFields(){
    return m_originals.size();
}

/*!
 * Return the number of all linked fields in your class. Multiple fields associated to a unique
 * geometry that need to be overlapped are counted too.
 * \return number of fields
 */
int
OverlapVectorFields::getNLinkedFields(){

    int sum = 0;

    for(auto & val : m_originals){
        sum +=val.second.size();
    }
    return sum;
}

/*!
 * Return list of geometry pointers actually linked
 * \return vector of geometry pointers
 */
std::vector<MimmoObject * > 
OverlapVectorFields::whichGeometriesLinked(){

    std::vector<MimmoObject * > res(m_originals.size());
    int counter = 0;
    for(auto & val : m_originals){
        res[counter] = val.first;
        ++counter;
    }
    return res;
}

/*!
 * Return overlapped fields after class execution as map of geometry pointers 
 * and data field pointers which refers to
 * \return map of geometry pointers and related fields
 */
std::unordered_map<MimmoObject*, dvecarr3E* > 
OverlapVectorFields::getDataFieldMap(){

    std::unordered_map<MimmoObject*, dvecarr3E* > res;
    for(auto & val : m_results){
        res[val.first] = &(val.second);
    }
    return res;
}

/*!
 * Return overlapped fields after class execution as a list of pairs having geometry pointers 
 * as key and data field pointers which refers to as argument
 * \return vector of pair of geometry pointers and related fields
 */
std::vector<std::pair<MimmoObject*, dvecarr3E* > >
OverlapVectorFields::getDataFieldList(){

    std::vector<std::pair<MimmoObject*, dvecarr3E* > > res(m_results.size());
    int counter = 0;
    for(auto & val : m_results){
        res[counter] = std::make_pair(val.first, &(val.second));
        ++counter;
    }
    return res;
}

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::SUM. Enum overloading
 * \param[in] funct    Type of overlap method
 */
void
OverlapVectorFields::setOverlapCriteriumENUM( OverlapMethod funct){
    setOverlapCriterium(static_cast<int>(funct));
    m_results.clear();
};

/*!
 * Set overlap criterium for multi-fields overlapping. See OverlapMethod enum
 * Class default is OverlapMethod::SUM. 
 * \param[in] funct    Type of overlap method
 */
void
OverlapVectorFields::setOverlapCriterium( int funct){
    if(funct <1 ||funct > 4)    return;
    m_overlapCriterium = static_cast<OverlapMethod>(funct);
    m_results.clear();
};

/*!
 * Append data fields of a geometry in a list , as a pair
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] field    Pair with geometry pointer and related field to be inserted
 */
void
OverlapVectorFields::setAddDataField( std::pair<MimmoObject*, dvecarr3E*> field){

    if(field.first == NULL || field.second == NULL) return;
    if(field.first->isEmpty() || field.second->empty()) return;

    m_originals[field.first].push_back(field.second);
    m_results.clear();
};

/*!
 * Append data fields of a geometry in a list , as a map of
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] fieldMap    map of geometry pointers and related fields to be inserted
 */
void
OverlapVectorFields::setDataFieldMap( std::unordered_map<MimmoObject*, dvecarr3E*> fieldMap){

    for(auto & val : fieldMap){
        setAddDataField(val);
    }
};

/*!
 * Append data fields of a geometry in a list , as a vector of pairs of
 * (pointer to the geometry, pointer to the field it refers to).
 * \param[in] fieldList    Vector of pair with geometry pointers and related fields to be inserted
 */
void
OverlapVectorFields::setDataFieldList(std::vector<std::pair<MimmoObject*, dvecarr3E*> > fieldList){

    for(auto & val : fieldList){
        setAddDataField(val);
    }
};

/*!
 * Remove a data field on the list by passing as key its pointer to geometry mesh
 * \param[in] patch Pointer to geometry to be removed.
 */
void
OverlapVectorFields::removeData(MimmoObject * patch){
    std::unordered_map<MimmoObject *, std::vector<dvecarr3E *> >::iterator it = m_originals.find(patch);
    if(it != m_originals.end()){
        m_originals.erase(it);
    }
    m_results.clear();
};

/*!
 * Remove all data present in the list
 */
void
OverlapVectorFields::removeAllData(){
    m_originals.clear();
    m_results.clear();
};

/*!
 * Clear your class data and restore defaults settings
 */
void
OverlapVectorFields::clear(){
    BaseManipulation::clear();
    removeAllData();
    m_overlapCriterium = OverlapMethod::SUM;
}

/*! 
 * It builds the input/output ports of the object
 */
void
OverlapVectorFields::buildPorts(){

    bool built = true;

    //input
    built = (built && createPortIn<std::pair<MimmoObject *, dvecarr3E *>, OverlapVectorFields>(this, &mimmo::OverlapVectorFields::setAddDataField, PortType::M_PAIRVECFIELD, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
    built = (built && createPortIn<std::unordered_map<MimmoObject *, dvecarr3E *>,OverlapVectorFields>(this, &mimmo::OverlapVectorFields::setDataFieldMap, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
    built = (built && createPortIn<std::vector<std::pair<MimmoObject *, dvecarr3E *>>,OverlapVectorFields>(this, &mimmo::OverlapVectorFields::setDataFieldList, PortType::M_VECPAIRVF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECARR3FLOAT_));

    //output
    built = (built && createPortOut<std::unordered_map<MimmoObject *, dvecarr3E *>,OverlapVectorFields>(this, &mimmo::OverlapVectorFields::getDataFieldMap, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
    built = (built && createPortOut<std::vector<std::pair<MimmoObject *, dvecarr3E *>>,OverlapVectorFields>(this, &mimmo::OverlapVectorFields::getDataFieldList, PortType::M_VECPAIRVF, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::PAIRMIMMO_VECARR3FLOAT_));
    m_arePortsBuilt = built;
};

/*!
 * Plot overlapped data field of a target geometry. If geometry is not listed into the class, plot nothing.
 * If overlapped data field is not available into results (because of plotting before proper execution), returns
 * a vtu file with 0 field.
 * \param[in]    dir        Output directory
 * \param[in]    name    Output filename
 * \param[in]    flag    Writing codex flag, false ascii, binary true
 * \param[in]    counter Counter identifying your output name
 * \param[in]    geo     Pointer to a geometry listed to the class.
 */
void
OverlapVectorFields::plotData(std::string dir, std::string name, bool flag, int counter, mimmo::MimmoObject * geo){

    if(m_originals.count(geo) == 0) return;
    dvecarr3E field;
    bool flagCalc = (m_results.count(geo) > 0);


    dvecarr3E points = geo->getVertexCoords();
    ivector2D connectivity;
    bitpit::VTKElementType cellType = geo->desumeElement();

    if (geo->getType() != 3){
        connectivity = geo->getCompactConnectivity();
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

    if(flagCalc)    field = m_results[geo];
    else            field.resize(points.size(), {{0.0,0.0,0.0}});

    output.addData("vectorfield", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, field);

    std::vector<long> ids(points.size());
    long ID;
    auto convMap = geo->getMapDataInv();
    for (auto vertex : geo->getVertices()){
        ID = vertex.getId();
        ids[convMap[ID]] = ID;
    }

    output.addData("ID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, ids);

    output.setCounter(counter);
    output.setCodex(bitpit::VTKFormat::APPENDED);
    if(!flag) output.setCodex(bitpit::VTKFormat::ASCII);

    output.write();
};

/*!
 * Plot all overlapped data field of all target geometries. Results will returned with the target filename and 
 * a consecutive numeration. If class is not executed and overlapped fields are not available plot nothing
 * \param[in]    dir        Output directory
 * \param[in]    name    Output unique filename
 * \param[in]    flag    Writing codex flag, false ascii, binary true
 */
void
OverlapVectorFields::plotAllData(std::string dir, std::string name, bool flag){

    int counter = 0;
    for(auto & val : m_results){
        plotData(dir, name, flag, counter, val.first);
        ++counter;
    }
};

/*
 * Execution command.
 * Overlap fields and save result in m_results member.
 */
void
OverlapVectorFields::execute(){

    if(m_originals.empty())    return;
    m_results.clear();
    int size, listsize;
    dvecarr3E temp;
    dvecarr3E results;
    for(auto & obj : m_originals){

        size = obj.first->getPatch()->getVertexCount();
        //careful resizing
        results.resize(size, {{0.0,0.0,0.0}});
        for(auto &vec : obj.second)    vec->resize(size,{{0.0,0.0,0.0}});

        listsize = obj.second.size();
        temp.resize(listsize);

        for(int i=0; i<size; ++i){
            for(int j=0; j<listsize; j++)    temp[j] = (*(obj.second[j]))[i];
            results[i] = overlapFields(temp);
        }

        m_results[obj.first] = results;
    }
}

/*!
 * Plot Optional results in execution, that is the series of overlapped field .
 */
void
OverlapVectorFields::plotOptionalResults(){
    std::string dir = m_outputPlot;
    std::string name = m_name;
    plotAllData(dir, name, true);
}

/*!
 * Overlap concurrent value of different fields in the same node. Overlap Method is specified
 * in the class set
 *\param[in] locField List of value of concurrent field. If value is unique, simply assigns it
 *\return    assigned value
 */
//DEVELOPERS REMIND if more overlap methods are added refer to this method to implement them
darray3E
OverlapVectorFields::overlapFields(dvecarr3E & locField){
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
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
OverlapVectorFields::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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
void OverlapVectorFields::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    int value = static_cast<int>(m_overlapCriterium);
    slotXML.set("OverlapCriterium", std::to_string(value));

};

}
