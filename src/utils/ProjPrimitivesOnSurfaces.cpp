/*----------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
 *  Copyright (c) 2015-2021 Optimad Engineering S.r.l., All Rights Reserved.
 *
 *  --------------------------------------------------------------------------
 *
 *  NOTICE:  All information contained herein is, and remains the property
 *  of COMPANY. The intellectual and technical concepts contained herein are
 *  proprietary to COMPANY and may be covered by Italian and Foreign Patents,
 *  patents in process, and are protected by trade secret or copyright law.
 *  Dissemination of this information or reproduction of this material is
 *  strictly forbidden unless prior written permission is obtained from
 *  COMPANY. Access to the source code contained herein is hereby forbidden
 *  to anyone except current COMPANY employees, managers or contractors who
 *  have executed Confidentiality and Non-disclosure agreements explicitly
 *  covering such access.
 *
 *  The copyright notice above does not evidence any actual or intended
 *  publication or disclosure of this source code, which includes information
 *  that is confidential and/or proprietary, and is a trade secret, of
 *  COMPANY. ANY REPRODUCTION, MODIFICATION, DISTRIBUTION, PUBLIC PERFORMANCE,
 *  OR PUBLIC DISPLAY OF OR THROUGH USE  OF THIS  SOURCE CODE  WITHOUT THE
 *  EXPRESS WRITTEN CONSENT OF COMPANY IS STRICTLY PROHIBITED, AND IN
 *  VIOLATION OF APPLICABLE LAWS AND INTERNATIONAL TREATIES. THE RECEIPT OR
 *  POSSESSION OF THIS SOURCE CODE AND/OR RELATED INFORMATION DOES NOT CONVEY
 *  OR IMPLY ANY RIGHTS TO REPRODUCE, DISCLOSE OR DISTRIBUTE ITS CONTENTS, OR
 *  TO MANUFACTURE, USE, OR SELL ANYTHING THAT IT  MAY DESCRIBE, IN WHOLE OR
 *  IN PART.
 *
\*----------------------------------------------------------------------------*/

#include "ProjPrimitivesOnSurfaces.hpp"

namespace mimmo{

/*!
 * Default constructor of ProjPrimitivesOnSurfaces.
 */
ProjPrimitivesOnSurfaces::ProjPrimitivesOnSurfaces(){
    m_name         = "";
    m_nC = 1000;
    m_topo = 0;
    m_buildSkdTree = false;
    m_buildKdTree = false;

}

/*!
 * Default destructor of ProjPrimitivesOnSurfaces.
 */
ProjPrimitivesOnSurfaces::~ProjPrimitivesOnSurfaces(){};

/*!
 * Copy constructor of ProjPrimitivesOnSurfaces. No resulting projected data structure is copied.
 */
ProjPrimitivesOnSurfaces::ProjPrimitivesOnSurfaces(const ProjPrimitivesOnSurfaces & other):BaseManipulation(other){
    m_topo = other.m_topo;
    m_nC = other.m_nC;
    m_buildSkdTree = other.m_buildSkdTree;
    m_buildKdTree = other.m_buildKdTree;
};

/*!
 * Assignement operator of ProjPrimitivesOnSurfaces. No resulting projected data structure is copied.
 */
ProjPrimitivesOnSurfaces & ProjPrimitivesOnSurfaces::operator=(const ProjPrimitivesOnSurfaces &other){
    m_topo = other.m_topo;
    m_nC = other.m_nC;
    m_buildSkdTree = other.m_buildSkdTree;
    m_buildKdTree = other.m_buildKdTree;
    *(static_cast<BaseManipulation *>(this)) = *(static_cast<const BaseManipulation *>(&other));
    return *this;
};

/*!
 * Swap function.
 * \param[in] x object ot be swapped
 */
void ProjPrimitivesOnSurfaces::swap(ProjPrimitivesOnSurfaces & x)   noexcept
{
    std::swap(m_topo, x.m_topo);
    std::swap(m_nC, x.m_nC);
    std::swap(m_buildSkdTree, x.m_buildSkdTree);
    std::swap(m_buildKdTree, x.m_buildKdTree);
    std::swap(m_patch, x.m_patch);
    BaseManipulation::swap(x);
}

/*!
 * Building ports of the class
 */
void
ProjPrimitivesOnSurfaces::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, ProjPrimitivesOnSurfaces>(this, &mimmo::ProjPrimitivesOnSurfaces::setGeometry,M_GEOM, true));

    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, ProjPrimitivesOnSurfaces>(this, &mimmo::ProjPrimitivesOnSurfaces::getProjectedElement, M_GEOM));
    m_arePortsBuilt = built;
}


/*!
 * Return the topology of the primitive element support. 0-0D(point clouds), 1-1D(3D curve), 2-2D (3D surfaces)
 * \return topology of the primitive element
 */
int
ProjPrimitivesOnSurfaces::getTopology(){
    return m_topo;
}

/*!
 * Get the actual number of the cells that will be used to represent a discrete mesh of the primitive projection.
 * \return number of cells
 */
int
ProjPrimitivesOnSurfaces::getProjElementTargetNCells(){
    return m_nC;
}


/*!
 * Get your current projected primitive as a 3D mesh
 * \return pointer to projected primitive as MimmoObject
 */
MimmoSharedPointer<MimmoObject>
ProjPrimitivesOnSurfaces::getProjectedElement(){
    return m_patch;
}



/*!
 * Set a target external surface mesh where primitive projection need to be performed.
 * Topology of the geometry must be of superficial type, so that MimmoObject::getType() must return 1, otherwise
 * nothing will be set.
 * \param[in] geo  pointer to MimmoObject
 */
void
ProjPrimitivesOnSurfaces::setGeometry(MimmoSharedPointer<MimmoObject> geo){
    if(geo == nullptr) return;
    if(geo->getType() != 1)    return;
    m_geometry = geo;
};


/*!It sets if the SkdTree of the projected primitive element data structure has to be built during execution.
 * \param[in] build If true the skdTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
ProjPrimitivesOnSurfaces::setBuildSkdTree(bool build){
    m_buildSkdTree = build;
}

/*!It sets if the KdTree the projected primitive element data structure has to be built during execution.
 * \param[in] build If true the KdTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
ProjPrimitivesOnSurfaces::setBuildKdTree(bool build){
    m_buildKdTree = build;
}

/*!It sets the number of the cells to be used to represent in a discrete mesh the primitive projection.
 * \param[in] nC number of cells
 */
void
ProjPrimitivesOnSurfaces::setProjElementTargetNCells(int nC){
    m_nC = nC;
}

/*!
 * Check if resulting primitive element is present or not.
 * True - no geometry evaluated, False otherwise.
 * \return empty flag
 */
bool
ProjPrimitivesOnSurfaces::isEmpty(){
    if(m_patch == nullptr) return true;
    return  m_patch->isEmpty();
}

/*!
 * Clear all stuffs in your class
 */
void
ProjPrimitivesOnSurfaces::clear(){
    m_patch.reset(nullptr);
    m_nC = 1000;
    BaseManipulation::clear();
};

/*!
 * Execution command.
 * Project the primitive on the target 3D surface.
 */
void
ProjPrimitivesOnSurfaces::execute(){

    if (m_geometry == nullptr){
        (*m_log)<<"Error in "<<m_name << " : nullptr pointer to linked geometry"<<std::endl;
        throw std::runtime_error (m_name + " : nullptr pointer to linked geometry");
    }

    projection();

    if(m_patch){
        if(m_buildSkdTree)   m_patch->buildSkdTree();
        if(m_buildKdTree)    m_patch->buildKdTree();
    }else{
        (*m_log)<<m_name << " : failed projecting object on target surface"<<std::endl;
    }
};

/*!
 * Plot resulting projected element in a vtu mesh file;
 */
void
ProjPrimitivesOnSurfaces::plotOptionalResults(){

	write(m_patch);

};


}
