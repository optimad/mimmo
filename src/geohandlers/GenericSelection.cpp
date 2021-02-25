/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

namespace mimmo {

/*!
 * Basic Constructor
 */
GenericSelection::GenericSelection(){
    m_type = SelectionType::UNDEFINED;
    m_topo = 1; /*default to surface geometry*/
    m_dual = false; /*default to exact selection*/
};

/*!
 * Basic Destructor
 */
GenericSelection::~GenericSelection(){
    m_subpatch.reset();
};

/*!
 * Copy Constructor, any already calculated selection is not copied.
 */
GenericSelection::GenericSelection(const GenericSelection & other):BaseManipulation(other){
    m_type = other.m_type;
    m_topo = other.m_topo;
    m_dual = other.m_dual;
};

/*!
 * Copy operator, any already calculated selection is not copied.
 */
GenericSelection & GenericSelection::operator=(const GenericSelection & other){
    *(static_cast<BaseManipulation *>(this)) = *(static_cast<const BaseManipulation * >(&other));
    m_type = other.m_type;
    m_topo = other.m_topo;
    m_dual = other.m_dual;
    /*m_subpatch is not copied and it is obtained in execution*/
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x GenericSelection object
 */
void GenericSelection::swap(GenericSelection & x) noexcept
{
    std::swap(m_type, x.m_type);
    std::swap(m_topo, x.m_topo);
    std::swap(m_dual, x.m_dual);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
GenericSelection::buildPorts(){

    bool built = true;

    built = (built && createPortIn<mimmo::MimmoSharedPointer<MimmoObject>, GenericSelection>(this, &GenericSelection::setGeometry, M_GEOM, true));
    built = (built && createPortIn<bool, GenericSelection>(this, &GenericSelection::setDual,M_VALUEB));

    built = (built && createPortOut<mimmo::MimmoSharedPointer<MimmoObject>, GenericSelection>(this, &GenericSelection::getPatch,M_GEOM));
    built = (built && createPortOut<livector1D, GenericSelection>(this, &GenericSelection::constrainedBoundary, M_VECTORLI));
    m_arePortsBuilt = built;
};


/*!
 * Return type of method for selection.
 * See SelectionType enum for available methods.
 * \return type of selection method
 */
SelectionType
GenericSelection::whichMethod(){
    return    m_type;
};

/*!
 * Return pointer by copy to sub-patch extracted by the class
 * \return pointer to MimmoObject extracted sub-patch
 */
mimmo::MimmoSharedPointer<MimmoObject>
GenericSelection::getPatch(){
    return    m_subpatch;
};

/*!
 * Return pointer by copy to subpatch extracted by the class [Const overloading].
 * \return pointer to MimmoObject extracted sub-patch
 */
const mimmo::MimmoSharedPointer<MimmoObject>
GenericSelection::getPatch() const{
    return    m_subpatch;
};

/*!
 * Set link to target geometry for your selection.
 * Reimplementation of mimmo::BaseManipulation::setGeometry();
 *  \param[in] target Pointer to MimmoObject with target geometry.
 */
void
GenericSelection::setGeometry( mimmo::MimmoSharedPointer<MimmoObject> target){
    if(target == nullptr)  return;
    m_geometry = target;
    /*set topology informations*/
    m_topo = target->getType();

};


/*!
 * Set your class behavior selecting a portion of a target geoemetry.
 * Given a initial set up, gets the dual
 * result (its negative) of current selection.
 * For instance, in a class extracting geometry inside the volume of an
 * elemental shape, gets all other parts of it not included in the shape.
 * For a class extracting portions of geometry
 * matching an external tessellation, gets all the other parts not matching it.
 * For a class extracting portion of
 * geometry identified by a PID list, get all other parts not marked by such list.
 * \param[in] flag Active/Inactive "dual" feature true/false.
 */
void
GenericSelection::setDual(bool flag ){
    m_dual = flag;
}

/*!
 * Return actual status of "dual" feature of the class. See setDual method.
 * \return  true/false for "dual" feature activated or not
 */
bool
GenericSelection::isDual(){
    return m_dual;
};

/*!
 * Return list of constrained boundary nodes (all those boundary nodes of
 * the subpatch extracted which are not part of the boundary of the mother
 * geometry)
  MPI version return a list of all subpatch physical border vertices not belonging to
  the physical border of the mother, ghost included (internal rank border faces are not accounted).
  This list is computed locally, and vary from rank to rank, according to mesh distribution
  among ranks.
 * \return list of physical boundary vertex ID
 */
livector1D
GenericSelection::constrainedBoundary(){

    if(getGeometry() == nullptr || getPatch() == nullptr)    return livector1D(0);
    std::unordered_map<long, std::set<int> > survivors;

    auto daughterBCells  = getPatch()->extractBoundaryFaceCellID(true);
    auto motherBCells = getGeometry()->extractBoundaryFaceCellID(true);

    //save all those daughter boundary faces not included in mother pot
    for(const auto & sT : daughterBCells){
        if(motherBCells.count(sT.first) > 0){
            for(const auto & val: sT.second){
                if(motherBCells[sT.first].count(val) > 0)   continue;
                survivors[sT.first].insert(val);
            }
        }else{
            survivors[sT.first] = daughterBCells[sT.first];
        }
    }
    motherBCells.clear();
    daughterBCells.clear();

    //get vertex of the cleaned daughter boundary.
    std::unordered_set<long> containerVert;
    containerVert.reserve(getPatch()->getPatch()->getVertexCount());
    for(const auto & sT  :survivors){
        const bitpit::Cell & cell = getPatch()->getPatch()->getCell(sT.first);
        for(const auto & val: sT.second){
            bitpit::ConstProxyVector<long> faceVertIds = cell.getFaceVertexIds(val);
            containerVert.insert(faceVertIds.begin(), faceVertIds.end());
        }
    }
    //move up in livector1D container.
    livector1D result;
    result.reserve(containerVert.size());
    result.insert(result.end(), containerVert.begin(), containerVert.end());

    return result;
};


/*!
 * Execute your object. A selection is extracted and trasferred in
 * an indipendent MimmoObject structure pointed by m_subpatch member
 */
void
GenericSelection::execute(){
    if(getGeometry() == nullptr) {
        (*m_log)<<m_name + " : nullptr pointer to target geometry found"<<std::endl;
        throw std::runtime_error (m_name + " : nullptr pointer to target geometry found");
    }

    m_subpatch.reset();

// extract all the interior cell satisfying the extraction criterium.
    livector1D extracted = extractSelection();

    /*Create subpatch.*/
    mimmo::MimmoSharedPointer<MimmoObject> temp(new MimmoObject(m_topo));
	
	/* Set the same tolerance of the mother geoemtry*/
    temp->setTolerance(getGeometry()->getTolerance());

    if (m_topo != 3){

        livector1D vertExtracted = getGeometry()->getVertexFromCellList(extracted);
        for(const auto & idV : vertExtracted){
            temp->addVertex(getGeometry()->getVertexCoords(idV), idV);
        }

        int rank;
        for(const auto & idCell : extracted){
            bitpit::Cell & cell = getGeometry()->getPatch()->getCell(idCell);
            rank = -1;
#if MIMMO_ENABLE_MPI
            rank = getGeometry()->getPatch()->getCellRank(idCell);
#endif
            temp->addCell(cell, idCell, rank);
        }

    }else{

        for(const auto & idV : extracted){
            temp->addVertex(getGeometry()->getVertexCoords(idV),idV);
        }

    }

    auto originalmap = getGeometry()->getPIDTypeListWNames();
    auto currentPIDmap = temp->getPIDTypeList();
    for(const auto & val: currentPIDmap){
        temp->setPIDName(val, originalmap[val]);
    }
    m_subpatch = temp;

    // Clean and Update selection patch. This will update even parallel structures if needed.
    m_subpatch->cleanGeometry();
    m_subpatch->update();

};

/*!
 * Plot optional result of the class in execution. It plots the selected patch
 * as standard vtk unstructured grid.
 */
void
GenericSelection::plotOptionalResults(){

	 write(getPatch());

}

}
