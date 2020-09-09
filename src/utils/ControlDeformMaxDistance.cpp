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
    m_violationField.swap(x.m_violationField);
    m_defField.swap(x.m_defField);
    BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
ControlDeformMaxDistance::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dmpvecarr3E*, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setDefField, M_GDISPLS));
    built = (built && createPortIn<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setLimitDistance, M_VALUED));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setGeometry, M_GEOM, true));

    built = (built && createPortOut<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolation, M_VALUED));
    built = (built && createPortOut<dmpvector1D*, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolationField, M_SCALARFIELD));
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

    double result = -1.0*std::numeric_limits<double>::max();
    for(const auto & val : m_violationField){
        result = std::fmax(result, val);
    }
#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX, m_communicator);
#endif
    return    result;
};


/*!
 * Return the violation distances of each point of deformed geometry, on the geometry itself. The info is available after class execution.
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation field values
 */
dmpvector1D *
ControlDeformMaxDistance::getViolationField(){
    return  &m_violationField;
};

/*!
 * Set the deformative field associated to each point of the target geometry.
 * Field resize occurs in execution, if point dimension between field and geoemetry does not match.
 * \param[in]    field of deformation
 */
void
ControlDeformMaxDistance::setDefField(dmpvecarr3E *field){
    if(!field) return;
    m_defField.clear();
    m_violationField.clear();
    m_defField = *field;
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

/*!
 * Set the link to MimmoObject target geometry. Geometry must be a 3D surface
 * (MimmoObject of type 1). Reimplemented from BaseManipulation::setGeometry().
 * \param[in] geo pointer to target geometry
 */
void
ControlDeformMaxDistance::setGeometry(MimmoSharedPointer<MimmoObject> geo){
    if (geo == nullptr) return;
    if (geo->getType() != 1 ) return;

    BaseManipulation::setGeometry(geo);
}

/*!Execution command. Calculate violation value and store it in the class member m_violation
 */
void
ControlDeformMaxDistance::execute(){

    MimmoSharedPointer<MimmoObject> geo = getGeometry();
    if(geo == nullptr){
        (*m_log)<<"Error "+m_name + " : null pointer to linked geometry found"<<std::endl;
        throw std::runtime_error("Error "+m_name + " : null pointer to linked geometry found");
    }

    if(geo->isEmpty()){
        (*m_log)<<"Warning in " + m_name + " : empty linked geometry found"<<std::endl;
    }

    bool check = m_defField.getGeometry() == geo;
    check = check && m_defField.getDataLocation() == MPVLocation::POINT;
    if(!check){
        (*m_log)<<"Warning in "<<m_name<<": Unsuitable deformation field linked"<<std::endl;
    }
    m_defField.completeMissingData({{0.0,0.0,0.0}});

    geo->buildSkdTree();

    m_violationField.clear();
    m_violationField.initialize(geo, MPVLocation::POINT, -1.0*std::numeric_limits<double>::max());

    std::vector<long> mapIDV = geo->getVerticesIds(); //take all internals and ghosts if partitioned mesh.
    //set points structure
    std::vector<std::array<double,3> > points;
    std::vector<double> normDef;
    points.reserve(mapIDV.size());
    normDef.resize(mapIDV.size(), 1.0E-08);

    int counter = 0;
    for (long idV : mapIDV){
        points.push_back(geo->getVertexCoords(idV) + m_defField[idV] );
        normDef[counter] = std::max(normDef[counter], 1.2*norm2(m_defField[idV]) );
        ++counter;
    }

    std::vector<double> distances(mapIDV.size());
    std::vector<long> suppCellIds(mapIDV.size());

#if MIMMO_ENABLE_MPI
    if (geo->getPatch()->isCommunicatorSet()){
        std::vector<int> suppCellRanks(mapIDV.size());
        skdTreeUtils::globalDistance(points.size(), points.data(), geo->getSkdTree(), suppCellIds.data(), suppCellRanks.data(), distances.data(), normDef.data(), false);
    }
    else
#endif
    {
        skdTreeUtils::distance(points.size(), points.data(), geo->getSkdTree(), suppCellIds.data(), distances.data(), normDef.data());
    }

    //transfer distance value inside m_violation field.(parallel case, ghost are already in)
    //Final value of violation is local distance of deformed point minus the offset m_maxDist
    // fixed by the user
    counter = 0;
    for( long id : mapIDV){
        m_violationField[id] =  (distances[counter] - m_maxDist);
        ++counter;
    }


    double violationMax = getViolation();
    //write resume file: in case of parallel version, only rank 0 writes.
#if MIMMO_ENABLE_MPI
    if(m_rank == 0){
#endif
        std::string logname = m_name+std::to_string(getId())+"_violation.log";
        std::ofstream log;
        log.open(logname);
        log<<"mimmo "<<m_name<<" resume file"<<std::endl;
        log<<std::endl;
        log<<std::endl;
        log<<" violation value : " << violationMax << std::endl;
        log.close();

#if MIMMO_ENABLE_MPI
    }
    MPI_Barrier(m_communicator); // other ranks stalled here, waiting 0 to finish writing.
#endif

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

    m_defField.setName("Deformation Field");
    m_violationField.setName("Violation Distance Field");
    write(getGeometry(), m_defField, m_violationField);

}



}
