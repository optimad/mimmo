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

/*!
 * Default constructor. Requires topo flag. 1-surface, 2-volume, 3-pointcloud.
 */
SplitVectorField::SplitVectorField(int topo):SplitField(topo){
    m_name = "mimmo.SplitVectorField";
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SplitVectorField::SplitVectorField(const bitpit::Config::Section & rootXML){

    std::string fallback_name = "ClassNONE";
    std::string fallback_topo = "-1";
    std::string input_name = rootXML.get("ClassName", fallback_name);
    std::string input_topo = rootXML.get("Topology", fallback_topo);
    input_name = bitpit::utils::trim(input_name);
    input_topo = bitpit::utils::trim(input_topo);

    int topo = std::stoi(input_topo);
    m_topo = std::max(1,topo);
    if (m_topo >3) m_topo = 1;

    m_name = "mimmo.SplitVectorField";

    if(input_name == "mimmo.SplitVectorField"){
        absorbSectionXML(rootXML);
    }else{
        (*m_log)<<"Warning in custom xml mimmo::SplitVectorField constructor. No valid xml data found"<<std::endl;
    };
}

/*!
 * Default destructor
 */
SplitVectorField::~SplitVectorField(){
    m_field.clear();
    m_result.clear();
}

/*!
 * Build the ports of the class;
 */
void
SplitVectorField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dvecarr3E, SplitVectorField>(this, &mimmo::SplitVectorField::setField, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT, true));
    built = (built && createPortOut<std::unordered_map<MimmoObject*, dvecarr3E* >, SplitVectorField>(this, &mimmo::SplitVectorField::getSplittedData, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
    SplitField::buildPorts();
    m_arePortsBuilt = built;
}

/*!
 * Get splitted field as a map having pointers to splitted original geometries 
 * and the portions of field relative to them. 
 * \return map with pointers to split geometries and fileds
 */
std::unordered_map<MimmoObject*, dvecarr3E *> 
SplitVectorField::getSplittedData(){

    std::unordered_map<MimmoObject*, dvecarr3E *> map;
    if(m_originals.empty() || m_geometry==NULL) return map;
    if(m_result.size() != m_originals.size())    return map;

    int counter = 0;
    for(auto & val : m_originals){
        map[val] = &(m_result[counter]);
        ++counter;
    }

    return map;
}

/*!
 * Set Field associated to the target geometry and that need to splitted.
 * If the field is associated to the cells or to points of the target geometry,
 * please set this info, choosing the correct division map between setCellDivisionMap or 
 * setVertDivisionMap methods.  
 * \param[in]    field vector field of array at 3 double elements
 */
void
SplitVectorField::setField(dvecarr3E field){
    m_field = field;
}

/*!
 * Clear content of the class
 */
void
SplitVectorField::clear(){
    m_field.clear();
    m_result.clear();
    SplitField::clear();
}

/*!
 * Plot splitted field alongside its geometries ;
 */
void 
SplitVectorField::plotOptionalResults(){
    if(isEmpty()) return;
    if(m_originals.size() != m_result.size())    return;

    bitpit::VTKLocation loc = bitpit::VTKLocation::POINT;
    if(m_mapVertDivision.empty())    loc = bitpit::VTKLocation::CELL;


    int counter=0;
    for(auto & geo : m_originals){

        bitpit::VTKElementType cellType = geo->desumeElement();
        dvecarr3E points = geo->getVertexCoords();

        if (cellType == bitpit::VTKElementType::UNDEFINED) continue;

        if(cellType != bitpit::VTKElementType::VERTEX){
            ivector2D connectivity = geo->getCompactConnectivity();
            bitpit::VTKUnstructuredGrid output(".",m_name+std::to_string(getClassCounter()),cellType);
            output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
            output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
            output.setDimensions(connectivity.size(), points.size());
            output.addData("field", bitpit::VTKFieldType::VECTOR, loc, m_result[counter]);
            output.setCounter(counter);
            output.setCodex(bitpit::VTKFormat::APPENDED);
            output.write();
        }else{
            int size = points.size();
            ivector1D connectivity(size);
            for(int i=0; i<size; ++i)    connectivity[i]=i;
            bitpit::VTKUnstructuredGrid output(".",m_name+std::to_string(getClassCounter()),    cellType);
            output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
            output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
            output.setDimensions(connectivity.size(), points.size());
            output.addData("field", bitpit::VTKFieldType::VECTOR, loc, m_result[counter]);
            output.setCounter(counter);
            output.setCodex(bitpit::VTKFormat::APPENDED);
            output.write();
        }

        ++counter;
    }
}

/*!
 * Split your original field along the splitted original geometries provided
 * \return true if split without errors
 */
bool
SplitVectorField::split(){
    if(m_mapCellDivision.empty() && m_mapVertDivision.empty())    return false;
    if(isEmpty())    return false;

    //check original field;
    bool loc = m_mapVertDivision.empty(); //true by cells, false by vert.

    int nGeo = m_originals.size();
    m_result.resize(nGeo);

    std::unordered_map<long, std::pair<int, long>> * mapp;
    livector1D * mapIdTarget;
    std::vector< std::map<long, int> * > invIdLoc(nGeo);

    if(loc)    {
        m_field.resize(getGeometry()->getNCells(), {{0.0,0.0,0.0}});
        mapp = &m_mapCellDivision;
        mapIdTarget = &(getGeometry()->getMapCell());
        for(int i=0; i<nGeo; ++i){
            invIdLoc[i] = &(m_originals[i]->getMapCellInv());
        }
    }
    else    {
        m_field.resize(getGeometry()->getNVertex(), {{0.0,0.0,0.0}});
        mapp = &m_mapVertDivision;
        mapIdTarget = &(getGeometry()->getMapData());
        for(int i=0; i<nGeo; ++i){
            invIdLoc[i] = &(m_originals[i]->getMapDataInv());
        }
    }

    if(m_field.size() != mapp->size())    return false;

    //allocate memory for results;
    for(int i=0; i<nGeo; ++i){

        if(loc)    m_result[i].resize(m_originals[i]->getNCells(),{{0.0,0.0,0.0}});
        else    m_result[i].resize(m_originals[i]->getNVertex(),{{0.0,0.0,0.0}});
    }

    int counter = 0;
    long IDtarget;
    std::pair<int,long> douple;
    int locali;
    for(auto & val: m_field){

        IDtarget = (*mapIdTarget)[counter];
        if(mapp->count(IDtarget)> 0 ){
            douple = (*mapp)[IDtarget];
            if(douple.first > nGeo)    continue;
            if((invIdLoc[douple.first])->count(douple.second) > 0){
                locali = (*invIdLoc[douple.first])[douple.second];
                m_result[douple.first][locali] = val;
            }
        }
        ++counter;
    }

    return true;
}

}
