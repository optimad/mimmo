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
#include "ControlDeformMaxDistance.hpp"

namespace mimmo{


/*!Default constructor of ControlDeformMaxDistance
 */
ControlDeformMaxDistance::ControlDeformMaxDistance(){
    m_name = "mimmo.ControlDeformMaxDistance";
    m_maxDist= 0.0 ;

};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ControlDeformMaxDistance::ControlDeformMaxDistance(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.ControlDeformMaxDistance";
    m_maxDist= 0.0 ;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.ControlDeformMaxDistance"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of ControlDeformMaxDistance
 */
ControlDeformMaxDistance::~ControlDeformMaxDistance(){};

/*!Copy constructor of ControlDeformMaxDistance.
 */
ControlDeformMaxDistance::ControlDeformMaxDistance(const ControlDeformMaxDistance & other):BaseManipulation(){
    *this = other;
};

/*!
 * Assignement operator of ControlDeformMaxDistance. Create an exact copy of the class,
 * except for the deformation field referred to the target geometry.
 */
ControlDeformMaxDistance & ControlDeformMaxDistance::operator=(const ControlDeformMaxDistance & other){
    *(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
    m_maxDist = other.m_maxDist;
    //deformation field is not copied
    return(*this);
};

/*! It builds the input/output ports of the object
 */
void
ControlDeformMaxDistance::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dvecarr3E, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setDefField, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT, true));
    built = (built && createPortIn<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setLimitDistance, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<MimmoObject*, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));

    built = (built && createPortOut<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolation, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<dvector1D, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolationField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<std::pair<BaseManipulation*, double>, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolationPair, PortType::M_VIOLATION, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::PAIRMIMMO_OBJFLOAT_));
    m_arePortsBuilt = built;
};

/*! 
 * Return the value of violation of deformed geometry, after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation value
 */
double 
ControlDeformMaxDistance::getViolation(){

    double result = -1.0E+18;
    for(auto & val : m_violationField){
        result = std::fmax(result, val);
    }

    return    result;
};


/*! 
 *  Return the value of violation of deformed geometry, after class execution. 
 *  A BaseManipulation object pointer, representing the sender of geometry to which violation is referred, 
 *  is attached to violation value; 
 *  See getViolation method doc for further details. 
 * \return std::pair structure reporting geometry sender pointer and violation value.
 */
std::pair<BaseManipulation*, double> 
ControlDeformMaxDistance::getViolationPair(){

    //get map of Input ports of the class.
    std::map<short int, mimmo::PortIn*> mapPorts = getPortsIn();

    //get class who send geometry here - portID = 99 -> M_GEOM

    std::vector<BaseManipulation*> senders = mapPorts[99]->getLink();

    std::string etiq;
    if(senders.size() == 0){
        return    std::make_pair(this, getViolation());
    }else{
        return    std::make_pair(senders[0], getViolation());
    }

};

/*! 
 * Return the violation distances of each point of deformed geometry, on the geometry itself. The info is available after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation field values
 */
dvector1D 
ControlDeformMaxDistance::getViolationField(){
    return(m_violationField);
};

/*!
 * Set the deformative field associated to each point of the target geometry. 
 * Field resize occurs in execution, if point dimension between field and geoemetry does not match.
 * \param[in]    field of deformation
 */
void
ControlDeformMaxDistance::setDefField(dvecarr3E field){
    m_defField.clear();
    m_violationField.clear();
    m_defField = field;
    m_violationField.resize(field.size(),-1.E+18);
};

/*! Set limit distance d of the constraint surface. Must be a positive definite value (>= 0).
 *  Given a target geometry surface (open or closed), its constraint surface is intended 
 *  as the iso-level surface where every points are at distance d from the target surface.
 * \param[in] dist limit distance value
 */
void 
ControlDeformMaxDistance::setLimitDistance(double dist){
    m_maxDist = std::fmax(1.0E-12,dist);
};

/*!Execution command. Calculate violation value and store it in the class member m_violation
 */
void
ControlDeformMaxDistance::execute(){

    MimmoObject * geo = getGeometry();
    if(geo->isEmpty()) return;
    if(!(geo->isBvTreeSupported())) return;

    m_defField.resize(getGeometry()->getNVertex(),darray3E{{0.0,0.0,0.0}});
    m_violationField.resize(m_defField.size());

    if(!(geo->isBvTreeBuilt()))    geo->buildBvTree();

    dvecarr3E points = geo->getVertexCoords();
    points+= m_defField;
    dvector1D normDef(m_defField.size());

    int count=0;
    for(auto & def : m_defField){
        normDef[count] = norm2(def);
        ++count;
    }

    double dist;
    double radius ;
    double rate = 0.05;
    int kmax = 200;
    int kiter;
    bool flag;
    long id;
    count = 0;

    for(auto &p : points){
        dist =1.0E+18;
        kiter = 0;
        flag = true;
        radius = std::fmax(1.0E-8, normDef[count]);
        while(flag && kiter < kmax){
            dist = bvTreeUtils::distance(&p, geo->getBvTree(), id, radius);
            flag = (dist == 1.0E+18);
            if(flag)    radius *= (1.0+ rate*((double)flag));
            kiter++;
        }
        if(kiter == kmax)    dist = m_maxDist - dist;
        m_violationField[count] =  (dist - m_maxDist);
        count++;
    }
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ControlDeformMaxDistance::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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

    if(slotXML.hasOption("LimitDistance")){
        std::string input = slotXML.get("LimitDistance");
        input = bitpit::utils::trim(input);
        double value = 0.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setLimitDistance(value);
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
void
ControlDeformMaxDistance::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));

    slotXML.set("LimitDistance", std::to_string(m_maxDist));

    if(isPlotInExecution()){
        slotXML.set("PlotInExecution", std::to_string(1));
    }

    if(m_outputPlot != "."){
        slotXML.set("OutputPlot", m_outputPlot);
    }

};

/*!
 * Plot optional results in execution, that is the violation distance field on current deformed geometry.
 * Reimeplemented from BaseManipulation::plotOptionalResults;
 */
void
ControlDeformMaxDistance::plotOptionalResults(){
    if(getGeometry()->isEmpty())    return;

    dvecarr3E  points = getGeometry()->getVertexCoords();
    m_defField.resize(points.size());
    points+=m_defField;
    ivector2D connectivity = getGeometry()->getCompactConnectivity();

    bitpit::VTKElementType  elDM = bitpit::VTKElementType::TRIANGLE;

    std::string name = m_name +std::to_string(getClassCounter())+ "_ViolationField";
    bitpit::VTKUnstructuredGrid output(m_outputPlot, name, elDM);
    output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points) ;
    output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity) ;
    output.setDimensions(connectivity.size(), points.size());

    m_violationField.resize(points.size(), -1.e+18);
    std::string sdfstr = "Violation Distance Field";
    output.addData(sdfstr, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_violationField);
    output.write();
}


}
