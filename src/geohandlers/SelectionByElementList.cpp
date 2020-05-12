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

#include "MeshSelection.hpp"

namespace mimmo{

/*!
 * Basic Constructor
 */
SelectionByElementList::SelectionByElementList(){
    m_name = "mimmo.SelectionByElementList";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByElementList::SelectionByElementList(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.SelectionByElementList";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.SelectionByElementList"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Destructor
 */
SelectionByElementList::~SelectionByElementList(){};


/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x SelectionByElementList object
 */
void SelectionByElementList::swap(SelectionByElementList & x) noexcept
{
   GenericSelection::swap(x);
   std::swap(m_annotatedcells, x.m_annotatedcells);
   std::swap(m_annotatedvertices, x.m_annotatedvertices);
   std::swap(m_rawcells, x.m_rawcells);
   std::swap(m_rawvertices, x.m_rawvertices);
}

/*!
 * It builds the input/output ports of the object.
 */
void
SelectionByElementList::buildPorts(){

    GenericSelection::buildPorts();
    bool built = m_arePortsBuilt;

    built = (built && createPortIn<MimmoPiercedVector<long> *, SelectionByElementList>(this, &SelectionByElementList::addAnnotatedCellList, M_LONGFIELD2,true,2 ));
    built = (built && createPortIn<MimmoPiercedVector<long> *, SelectionByElementList>(this, &SelectionByElementList::addAnnotatedVertexList, M_LONGFIELD,true,2 ));
    built = (built && createPortIn<std::vector<long>, SelectionByElementList>(this, &SelectionByElementList::addRawCellList, M_VECTORLI3,true,2 ));
    built = (built && createPortIn<std::vector<long>, SelectionByElementList>(this, &SelectionByElementList::addRawVertexList, M_VECTORLI2,true,2 ));

    m_arePortsBuilt = built;
};

/*!
    Add a list of mesh cells as candidates for selection.
    Multiple list can be added.
    \param[in] celldata list of mesh cell ids.
*/
void
SelectionByElementList::addAnnotatedCellList(MimmoPiercedVector<long> * celldata){
    if(!celldata) return;
    if(celldata->getDataLocation() != MPVLocation::CELL)   return;
    m_annotatedcells.push_back(celldata);
}

/*!
    Add a list of mesh vertices. Mesh cells strictly sharing those vertices will be
    marked as candidates for selection.
    Multiple list can be added.
    \param[in] vertexdata list of mesh vertex ids.
*/
void
SelectionByElementList::addAnnotatedVertexList(MimmoPiercedVector<long> * vertexdata){
    if(!vertexdata) return;
    if(vertexdata->getDataLocation() != MPVLocation::POINT)   return;
    m_annotatedvertices.push_back(vertexdata);
}

/*!
    Add a list of mesh cells as candidates for selection.
    Multiple list can be added.
    \param[in] celldata list of mesh cell ids.
*/
void
SelectionByElementList::addRawCellList(std::vector<long> celldata){
    if(celldata.empty())    return;
    m_rawcells.push_back(celldata);
}

/*!
    Add a list of mesh vertices. Mesh cells strictly sharing those vertices will be
    marked as candidates for selection.
    Multiple list can be added.
    \param[in] vertexdata list of mesh vertex ids.
*/
void
SelectionByElementList::addRawVertexList(std::vector<long> vertexdata){
    if(vertexdata.empty())    return;
    m_rawvertices.push_back(vertexdata);
}


/*!
 * Clear content of the class
 */
void
SelectionByElementList::clear(){
    m_subpatch.reset();
    BaseManipulation::clear();
    m_annotatedcells.clear();
    m_annotatedvertices.clear();
    m_rawcells.clear();
    m_rawvertices.clear();
};



/*!
 * Extract portion of target geometry which are enough near
 * to the external geometry provided
 * \return ids of cell of target tessellation extracted
 */
livector1D
SelectionByElementList::extractSelection(){

    mimmo::MimmoSharedPointer<MimmoObject> geo = getGeometry();
    std::unordered_set<long> setstruct;
    //fill a non repeated list of vertices ids
    for(MimmoPiercedVector<long> * ptr : m_annotatedvertices){
        std::vector<long> temp = ptr->getIds();
        setstruct.insert(temp.begin(), temp.end());
    }
    for(std::vector<long> list: m_rawvertices){
        setstruct.insert(list.begin(), list.end());
    }

    //check them w.r.t the geometry ids, then get their intersection vertexSureList;
    long id;
    livector1D vertexSureList;
    vertexSureList.reserve(std::min(setstruct.size(), std::size_t(geo->getNVertices())));
    for(auto it=geo->getPatch()->vertexBegin(); it!=geo->getPatch()->vertexEnd(); ++it){
        id = it.getId();
        if(setstruct.count(id) > 0){
            vertexSureList.push_back(id);
        }
    }
    vertexSureList.shrink_to_fit();

    //check if you are dealing with point clouds. if it's the case, return the vertex list;
    if(m_topo == 3)     return vertexSureList;

    //otherwise you need to do the same for cells.
    setstruct.clear();

    // push the cells identified by vertices in the unordered struct.
    livector1D addOtherCells = geo->getCellFromVertexList(vertexSureList, true);
    setstruct.insert(addOtherCells.begin(), addOtherCells.end());
    vertexSureList.clear();
    addOtherCells.clear();

    //fill the list with the cell ids from inputs.
    for(MimmoPiercedVector<long> * ptr : m_annotatedcells){
        std::vector<long> temp = ptr->getIds();
        setstruct.insert(temp.begin(), temp.end());
    }
    for(std::vector<long> list: m_rawcells){
        setstruct.insert(list.begin(), list.end());
    }

    //check them w.r.t the geometry ids, then get their intersection cellSureList;
    livector1D cellSureList;
    cellSureList.reserve(std::min(setstruct.size(), std::size_t(geo->getNCells())));
    for(auto it=geo->getPatch()->cellBegin(); it!=geo->getPatch()->cellEnd(); ++it){
        id = it.getId();
        if(setstruct.count(id) > 0){
            cellSureList.push_back(id);
        }
    }

    cellSureList.shrink_to_fit();

    return cellSureList;
};


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByElementList::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Dual")){
        std::string input = slotXML.get("Dual");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setDual(value);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByElementList::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML,name);

    slotXML.set("Dual", std::to_string(m_dual));

};

}
