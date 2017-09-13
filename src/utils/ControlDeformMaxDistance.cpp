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
#include "SkdTreeUtils.hpp"

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
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.ControlDeformMaxDistance"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of ControlDeformMaxDistance
 */
ControlDeformMaxDistance::~ControlDeformMaxDistance(){};

/*!Copy constructor of ControlDeformMaxDistance.Deformation field referred to geometry 
 * and result violation field are not copied.
 */
ControlDeformMaxDistance::ControlDeformMaxDistance(const ControlDeformMaxDistance & other):BaseManipulation(other){
    m_maxDist = other.m_maxDist;
};

/*!
 * Assignment operator of ControlDeformMaxDistance. Deformation field referred to geometry 
 * and result violation field are not copied.
 */
ControlDeformMaxDistance & ControlDeformMaxDistance::operator=(ControlDeformMaxDistance other){
    swap(other);
    return(*this);
};

/*!
 * Swap Function.
 * \param[in] x object ot be swapped
 */
void ControlDeformMaxDistance::swap(ControlDeformMaxDistance & x) noexcept
{
    std::swap(m_maxDist, x.m_maxDist);
    std::swap(m_violationField, x.m_violationField);
    std::swap(m_defField, x.m_defField);
    BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
ControlDeformMaxDistance::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dmpvecarr3E, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setDefField, PortType::M_GDISPLS, mimmo::pin::containerTAG::MPVECARR3, mimmo::pin::dataTAG::FLOAT, true));
    built = (built && createPortIn<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setLimitDistance, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<MimmoObject*, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));

    built = (built && createPortOut<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolation, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<dmpvector1D, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolationField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::MPVECTOR, mimmo::pin::dataTAG::FLOAT));
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
    for(const auto & val : m_violationField){
        result = std::fmax(result, val);
    }

    return    result;
};


/*! 
 * Return the violation distances of each point of deformed geometry, on the geometry itself. The info is available after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation field values
 */
dmpvector1D
ControlDeformMaxDistance::getViolationField(){
    return(m_violationField);
};

/*!
 * Set the deformative field associated to each point of the target geometry. 
 * Field resize occurs in execution, if point dimension between field and geoemetry does not match.
 * \param[in]    field of deformation
 */
void
ControlDeformMaxDistance::setDefField(dmpvecarr3E field){
    m_defField.clear();
    m_violationField.clear();
    m_defField = field;
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
    if (geo == NULL) return;
    if(geo->isEmpty() || m_defField.getGeometry() != geo) return;
    if(!(geo->isBvTreeSupported())) return;

    m_violationField.clear();


    if(!(geo->isBvTreeBuilt()))    geo->buildBvTree();

    dmpvecarr3E points;
    dmpvector1D normDef;
    long int ID;
    for (const auto & v : geo->getVertices()){
        ID = v.getId();
        points.insert(ID, geo->getVertexCoords(ID) + m_defField[ID] );
        normDef.insert(ID, norm2(m_defField[ID]));
        m_violationField.insert(ID, -1.0e+18);
    }

    double dist;
    double radius ;
    double rate = 0.05;
    int kmax = 200;
    int kiter;
    bool flag;
    long int IDC;
    for(const auto &ID : points.getIds()){
        dist =1.0E+18;
        kiter = 0;
        flag = true;
        radius = std::fmax(1.0E-8, normDef[ID]);
        while(flag && kiter < kmax){
            dist = skdTreeUtils::distance(&points[ID], geo->getBvTree(), IDC, radius);
            flag = (dist == 1.0E+18);
            if(flag)    radius *= (1.0+ rate*((double)flag));
            kiter++;
        }
        if(kiter == kmax)    dist = m_maxDist - dist;
        m_violationField[ID] =  (dist - m_maxDist);
    }

    m_violationField.setGeometry(getGeometry());
//     m_violationField.setName("violation");

    //write log
    std::string logname = m_name+"_violation";
    bitpit::Logger* log = &bitpit::log::cout(logname);
    log->setPriority(bitpit::log::DEBUG);
    (*log)<<" violation value : " << getViolation() << std::endl;

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
    BaseManipulation::absorbSectionXML(slotXML, name);
    
    if(slotXML.hasOption("LimitDistance")){
        std::string input = slotXML.get("LimitDistance");
        input = bitpit::utils::string::trim(input);
        double value = 0.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setLimitDistance(value);
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

    BaseManipulation::flushSectionXML(slotXML, name);
    
    slotXML.set("LimitDistance", std::to_string(m_maxDist));

};

/*!
 * Plot optional results in execution, that is the violation distance field on current deformed geometry.
 * Reimeplemented from BaseManipulation::plotOptionalResults;
 */
void
ControlDeformMaxDistance::plotOptionalResults(){
    if(getGeometry() == NULL) return;
    if(getGeometry()->isEmpty())    return;

    liimap map;
    dvecarr3E  points = getGeometry()->getVertexCoords(&map);
    dvecarr3E deff = m_defField.getDataAsVector();
    points+=deff;
    ivector2D connectivity = getGeometry()->getCompactConnectivity(map);

    bitpit::VTKElementType  elDM = bitpit::VTKElementType::TRIANGLE;

    std::string name = m_name +std::to_string(getClassCounter())+ "_ViolationField";
    bitpit::VTKUnstructuredGrid output(m_outputPlot, name, elDM);
    output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points) ;
    output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity) ;
    output.setDimensions(connectivity.size(), points.size());

    dvector1D viol = m_violationField.getDataAsVector();
    std::string sdfstr = "Violation Distance Field";
    output.addData(sdfstr, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, viol);
    output.write();
}


}
