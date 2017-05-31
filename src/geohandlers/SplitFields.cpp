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

#include "SplitFields.hpp"

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!Default constructor of SplitField.
 * Format admissible are linked to your choice of topology. See FileType enum.
 * No argument passed, give default topology 1.
 * \param[in] topo    Topology of your geometries. 1-surface, 2-volume, 3-pointcloud, 4-3D curve
 */
SplitField::SplitField(int topo){
    m_topo = std::max(1,topo);
    if (m_topo >4) m_topo = 1;
}

/*!
 * Default destructor of SplitField.
 */
SplitField::~SplitField(){
    clear();
};

/*!Copy constructor of SplitField.Soft Copy of MimmoObject;
 */
SplitField::SplitField(const SplitField & other):BaseManipulation(){
    *this = other;
};

/*!
 * Assignement operator of SplitField. Soft copy of MimmoObject
 */
SplitField & SplitField::operator=(const SplitField & other){
    clear();
    *(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
    m_topo = other.m_topo;
    m_originals = other.m_originals;
    m_mapCellDivision = other.m_mapCellDivision;
    m_mapVertDivision = other.m_mapVertDivision;

    return *this;
};

/*!
 * Build the ports of the class;
 */
void
SplitField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, SplitField>(this, &mimmo::SplitField::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    built = (built && createPortIn<std::vector<MimmoObject*>, SplitField>(this, &mimmo::SplitField::setSplittedGeometries, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_, true));

    built = (built && createPortIn<std::unordered_map<long,std::pair<int, long> >, SplitField>(this, &mimmo::SplitField::setCellDivisionMap, PortType::M_MAPDCELL, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGPAIRINTLONG));
    built = (built && createPortIn<std::unordered_map<long,std::pair<int, long> >, SplitField>(this, &mimmo::SplitField::setVertDivisionMap, PortType::M_MAPDVERT, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGPAIRINTLONG));
    m_arePortsBuilt = built;
}


/*!
 * Set the target geometry where your field is defined.Topology of the geometry must be coeherent
 * with topology of the class; Reimplementation of BaseManipulation::setGeometry
 * \param[in] geo  Pointer to MimmoObject
 */
void
SplitField::setGeometry(MimmoObject* geo){
        if(geo->isEmpty()) return;
        if(geo->getType() != m_topo)    return;

        m_geometry = geo;
};

/*!
 * Set original split geometries which refers to target geometry. List will be saved as is, replacing any other saved list.
 * \param[in] originals  Vector of pointers to MimmoObjects (split geometries)
 */
void
SplitField::setSplittedGeometries(std::vector<MimmoObject *> originals){

    if(originals.empty()) return;
    m_originals.resize(originals.size());

    int counter = 0;
    for(auto & obj: originals){
        if(obj->isEmpty()) continue;
        if(obj->getType() == m_topo){
            m_originals[counter]= obj;
            ++counter;
        }
    }
    m_originals.resize(counter);
};

/*!
 * It sets the Cell division map relative to the target geometry w.r.t the original split geometries
 * The class interprets this command as an explicit declaration that your current field to be split is
 * referred to geometry-cells. If any previous compiled Vertex Division map is set, it will be erased.
 * If class topology is a point cloud, the method do nothing. Please use setVertDivisionMap instead;
 * \param[in] map Cell division map of the target geometry w.r.t the original split geometries
 */
void
SplitField::setCellDivisionMap(std::unordered_map<long, std::pair<int,long > > map){
    if(m_topo == 3) return;
    m_mapCellDivision = map;
    m_mapVertDivision.clear();
}

/*!
 * It sets the Vertex division map relative to the target geometry w.r.t the original split geeometries
 * The class interprets this command as an explicit declaration that your current field to be split is
 * referred to geometry-vertices/points. If any previous compiled Cell Division map is set, it will be erased.
 * \param[in] map vertex division map of the target geometry w.r.t the original split geometries
 */
void
SplitField::setVertDivisionMap(std::unordered_map<long, std::pair<int,long > > map){

    m_mapVertDivision = map;
    m_mapCellDivision.clear();
}
/*!
 * Check if target geometry and its split originals are present or not.
 * \return true - no geometry present, false otherwise.
 */
bool 
SplitField::isEmpty(){
    return (m_geometry == NULL || m_originals.empty());
}

/*!
 * Return current topology type set for your class geometries.
 * \return topology (1-surface, 2-volume, 3-points cloud)
 */
int 
SplitField::getTopo(){
    return m_topo;
}

/*!
 * Clear all stuffs in your class
 */
void
SplitField::clear(){
    m_originals.clear();
    m_mapCellDivision.clear();
    m_mapVertDivision.clear();
    BaseManipulation::clear();
};

/*!Execution command.
 * It stitches together multiple geoemetries in the same object.
 */
void
SplitField::execute(){

    bool check = split();
    if(!check){
        std::cout<<"Error in class "<<m_name<<". Field cannot be split"<<std::endl;
        std::cout<<"This could be due to not correct setting of geometries or division maps"<<std::endl;
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SplitField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    //checking topology
    if(slotXML.hasOption("Topology")){
        std::string input = slotXML.get("Topology");
        input = bitpit::utils::trim(input);
        int temptop = -1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temptop;
        }
        if(m_topo != temptop)    return;
    }

    if(slotXML.hasOption("Priority")){
        input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

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
void SplitField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));
    slotXML.set("Topology", m_topo);

    std::string output;

    if(isPlotInExecution()){
        slotXML.set("PlotInExecution", std::to_string(1));
    }

    if(m_outputPlot != "."){
        slotXML.set("OutputPlot", m_outputPlot);
    }
};

}
