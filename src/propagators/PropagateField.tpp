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

#include "customOperators.hpp"
#include "SkdTreeUtils.hpp"

namespace mimmo{

/*!
 * Constructor
 */
template<std::size_t NCOMP>
PropagateField<NCOMP>::PropagateField(){
	m_name = "mimmo.PropagateField";
	setDefaults();
};

/*!
 * Set most significant parameters to constructor defaults
 */
template<std::size_t NCOMP>
void PropagateField<NCOMP>::setDefaults(){
    this->m_thres = 1.0E-8;
    this->m_tol   = 1.0E-12;
    this->m_print = false;

    this->m_dampingActive = false;
    this->m_dampingType = 0;
	this->m_decayFactor = 1.0;
    this->m_radius = 0.0;
    this->m_plateau = 0.0;

    this->m_bandActive = false;
    this->m_bandwidth = 0.0;
    this->m_bandrelax = 1.0;

// #if MIMMO_ENABLE_MPI
//     this->m_ghostTag = 0;
//     this->m_pointGhostTag = 0;
// #endif
}

/*!
 * Destructor;
 */
template<std::size_t NCOMP>
PropagateField<NCOMP>::~PropagateField(){};

/*!
 * Copy constructor
 */
template<std::size_t NCOMP>
PropagateField<NCOMP>::PropagateField(const PropagateField<NCOMP> & other):BaseManipulation(other){
    setDefaults();
    this->m_thres   = other.m_thres;
    this->m_tol     = other.m_tol;
    this->m_print   = other.m_print;
    this->m_field   = other.m_field;

    this->m_dirichletPatches = other.m_dirichletPatches;
    this->m_dirichletBcs      = other.m_dirichletBcs;

    this->m_dampingActive= other.m_dampingActive;
    this->m_dampingType  = other.m_dampingType;
    this->m_decayFactor  = other.m_decayFactor;
    this->m_radius       = other.m_radius;
    this->m_plateau      = other.m_plateau;
    this->m_dampingSurfaces = other.m_dampingSurfaces;

    this->m_bandActive   = other.m_bandActive;
    this->m_bandwidth    = other.m_bandwidth;
    this->m_bandrelax    = other.m_bandrelax;
    this->m_bandSurfaces = other.m_bandSurfaces;

};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
template<std::size_t NCOMP>
void PropagateField<NCOMP>::swap(PropagateField<NCOMP> & x) noexcept {

    std::swap(this->m_thres, x.m_thres);
    std::swap(this->m_tol, x.m_tol);
    std::swap(this->m_print, x.m_print);
    this->m_field.swap(x.m_field);

    std::swap(this->m_dirichletPatches, x.m_dirichletPatches);
    std::swap(this->m_dirichletBcs, x.m_dirichletBcs);

    std::swap(this->m_dampingActive, x.m_dampingActive);
    std::swap(this->m_dampingType, x.m_dampingType);
    std::swap(this->m_decayFactor, x.m_decayFactor);
    std::swap(this->m_radius, x.m_radius);
    std::swap(this->m_plateau, x.m_plateau);
    std::swap(this->m_dampingSurfaces, x.m_dampingSurfaces);
    this->m_damping.swap(x.m_damping);

    std::swap(this->m_bandActive, x.m_bandActive);
    std::swap(this->m_bandwidth, x.m_bandwidth);
    std::swap(this->m_bandrelax, x.m_bandrelax);
    std::swap(this->m_bandSurfaces, x.m_bandSurfaces);
    this->m_banddistances.swap(x.m_banddistances);

    this->BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, PropagateField<NCOMP> >(this, &PropagateField<NCOMP>::setGeometry, M_GEOM, true));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, PropagateField<NCOMP> >(this, &PropagateField<NCOMP>::addDirichletBoundaryPatch, M_GEOM2, true));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, PropagateField<NCOMP> >(this, &PropagateField<NCOMP>::addDampingBoundarySurface, M_GEOM3));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, PropagateField<NCOMP> >(this, &PropagateField<NCOMP>::addNarrowBandBoundarySurface, M_GEOM7));
    m_arePortsBuilt = built;
};

/*!
 * It sets the tolerance on residuals for solver convergence for laplacian solver.
 * \param[in] tol Convergence tolerance.
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::setTolerance(double tol){
    m_tol = std::max(std::numeric_limits<double>::min(), tol);;
}

/*!
 * It sets a lower threshold such that all the cells who have a norm of solution field
 * above it will be marked for update purposes.
 * \param[in] thres update lower threshold
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::setUpdateThreshold(double thres){
    m_thres = std::max(std::numeric_limits<double>::min(), thres);
}

/*!
 * If true print residuals and info during system solving.
   The option is available only for MPI compilation.
 * \param[in] print Print flag
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::setPrint(bool print){
    bool check = false;
#if MIMMO_ENABLE_MPI
    check = print;
#else
    BITPIT_UNUSED(print);
#endif
    m_print = check;
}

/*!
 * Set pointer to your target bulk geometry. Reimplemented from mimmo::BaseManipulation::setGeometry().
 * Geometry must be a of volume or surface type (MimmoObject type = 2 and type = 1);
 * \param[in] geometry_ pointer to target geometry
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setGeometry(MimmoSharedPointer<MimmoObject> geometry_){

    if (geometry_ == nullptr) return;
    if (geometry_->getType() != 2 && geometry_->getType() != 1){
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Warning: "<<m_name<<" allows only surface or volume bulk mesh. Skip set geometry."<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        return;
    }

    m_geometry = geometry_;

    if(!m_geometry->areAdjacenciesBuilt()){
        m_geometry->buildAdjacencies();
    }

}

/*!
 * Add a portion of boundary mesh relative to geometry target
 * that must be constrained with Dirichlet conditions. See method
   addDirichletConditions to properly provide a field of bc.
 * \param[in] bsurface Dirichlet boundary patch.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::addDirichletBoundarySurface(MimmoSharedPointer<MimmoObject> bsurface){
    if (bsurface == nullptr)       return;

    m_dirichletPatches.insert(bsurface);
}

/*!
 * Add a portion of boundary mesh relative to geometry target
 * that must be constrained with Dirichlet conditions. See method
   addDirichletConditions to properly provide a field of bc.
 * \param[in] bsurface Dirichlet boundary patch.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::addDirichletBoundaryPatch(MimmoSharedPointer<MimmoObject> bsurface){
    if (bsurface == nullptr)       return;

    m_dirichletPatches.insert(bsurface);
}

/*!
 * Activate Narrow Band Control(see class doc).
 * \param[in] flag boolean true activate, false deactivate.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setNarrowBand(bool flag){
	m_bandActive = flag;
}

/*!
   Add a portion of boundary mesh relative to geometry target
   to activate Narrow Band Control in its neighborhood (relaxation of Laplace solution).
   These patches are optional. If nothing is linked, Narrow Band Control (NBC) remains unactive.
   Please rember to specify the width of the narrow band if the NBC is active.
   \param[in] surface Boundary patch.
 */
 template <std::size_t NCOMP>
 void
 PropagateField<NCOMP>::addNarrowBandBoundarySurface(MimmoSharedPointer<MimmoObject> surface){
	if (surface == nullptr)      return;
	if (surface->getType()!= 1 ){
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Warning: "<<m_name<<" allows only narrouw band boundary surfaces. Skipping input slip patch."<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        return;
	}

    m_bandSurfaces.insert(surface);
}

/*!
    Set the width of narrow band for Narrow Band Control method near
    boundary surfaces (those set with addNarrowBandBoundarySurface method).
    \param[in] width of the narrow band (>0.0)
*/
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setNarrowBandWidth(double width){
    m_bandwidth = std::max(width, 0.0);
}

/*!
    Set the value of relaxation parameter for Narrow Band Control near
    boundary surfaces (those set with addNarrowBandBoundarySurface method).
    \param[in] relax  relaxation parameter within [0,1], where 1 means no relaxation, 0 max relaxation.
*/
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setNarrowBandRelaxation(double relax){
    m_bandrelax = std::max(0.0,std::min(relax, 1.0));
}

/*!
 * Activate Damping control by means of artificial diffusivity(see class doc).
 * \param[in] flag boolean true activate, false deactivate.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDamping(bool flag){
	m_dampingActive = flag;
}

/*!
 * If Damping is active, set its type, if distance control based (0) or volume cell control based (1).
 * \param[in] type of control for dumping [0,1].
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDampingType(int type){
	m_dampingType = std::max(0, std::min(1, type));
}


/*!
 * Set the Damping decay factor.
 * \param[in] decay exponential decaying factor for damping function.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDampingDecayFactor(double decay){
	m_decayFactor = decay;
}


/*!
 * Add a portion of boundary mesh to be used for damping function/artificial
   diffusivity calculation.
 * \param[in] bdumping Boundary patch.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::addDampingBoundarySurface(MimmoSharedPointer<MimmoObject> bdumping){
	if (bdumping == nullptr)       return;
	if (bdumping->getType()!= 1 ){
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Warning: "<<m_name<<" allows only slip boundary surfaces. Skipping input slip patch."<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        return;
	}

	m_dampingSurfaces.insert(bdumping);
}


/*!
 * Set the Damping Inner distance p (see class doc).
 * \param[in] plateau inner distance.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDampingInnerDistance(double plateau){
	m_plateau = plateau;
}

/*!
 * Set the Damping Outer distance r (see class doc).
 * \param[in] radius outer distance.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::setDampingOuterDistance(double radius){
	m_radius = radius;
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){


    BITPIT_UNUSED(name);

    //start absorbing
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Tolerance")){
        std::string input = slotXML.get("Tolerance");
        input = bitpit::utils::string::trim(input);
        double value = 1.0E-12;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::fmax(0.0, value);
        }
        setTolerance(value);
    }

    if(slotXML.hasOption("UpdateThres")){
        std::string input = slotXML.get("UpdateThres");
        input = bitpit::utils::string::trim(input);
        double value = 1.0E-8;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setUpdateThreshold(value);
    }

    if(slotXML.hasOption("Print")){
        std::string input = slotXML.get("Print");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setPrint(value);
    }


    if(slotXML.hasOption("NarrowBand")){
        std::string input = slotXML.get("NarrowBand");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setNarrowBand(value);
    }

    if(m_bandActive){
        if(slotXML.hasOption("NarrowBandWidth")){
            std::string input = slotXML.get("NarrowBandWidth");
            input = bitpit::utils::string::trim(input);
            double value = 0.0;
            if(!input.empty()){
                std::stringstream ss(input);
                ss >> value;
                value = std::fmax(0.0, value);
            }
            setNarrowBandWidth(value);
        }
        if(slotXML.hasOption("NarrowBandRelax")){
            std::string input = slotXML.get("NarrowBandRelax");
            input = bitpit::utils::string::trim(input);
            double value = 1.0;
            if(!input.empty()){
                std::stringstream ss(input);
                ss >> value;
                value = std::fmax(0.0, value);
            }
            setNarrowBandRelaxation(value);
        }
    }

    std::string inputDamping;
    if(slotXML.hasOption("Damping")){
        inputDamping = slotXML.get("Damping");
    }else if(slotXML.hasOption("Dumping")){  //-->legacy check
        inputDamping = slotXML.get("Dumping");
    }
    if(!inputDamping.empty()){
        std::string input = bitpit::utils::string::trim(inputDamping);
        bool value = false;
        if(!input.empty()){
           std::stringstream ss(input);
           ss >> value;
        }
       setDamping(value);
    }

    if(m_dampingActive){
        if(slotXML.hasOption("DecayFactor")){
            std::string input = slotXML.get("DecayFactor");
            input = bitpit::utils::string::trim(input);
            double value = 1.0;
            if(!input.empty()){
                std::stringstream ss(input);
                ss >> value;
                value = std::fmax(0.0, value);
            }
            setDampingDecayFactor(value);
        }

        std::string inputDampingType;
        if(slotXML.hasOption("DampingType")){
            inputDampingType = slotXML.get("DampingType");
        }else if(slotXML.hasOption("DumpingType")){ // --> legacy check
            inputDampingType = slotXML.get("DumpingType");
        }
        if(!inputDampingType.empty()){
            std::string input = bitpit::utils::string::trim(inputDampingType);
            int value = 0;
            if(!input.empty()){
                std::stringstream ss(input);
                ss >> value;
            }
            setDampingType(value);
        }
        std::string inputDampingInnerDistance;
        if(slotXML.hasOption("DampingInnerDistance")){
            inputDampingInnerDistance = slotXML.get("DampingInnerDistance");
        }else if(slotXML.hasOption("DumpingInnerDistance")){ // --> legacy check
            inputDampingInnerDistance = slotXML.get("DumpingInnerDistance");
        }
        if(!inputDampingInnerDistance.empty()){
            std::string input = bitpit::utils::string::trim(inputDampingInnerDistance);
            double value = 0.0;
            if(!input.empty()){
                std::stringstream ss(input);
                ss >> value;
                value = std::fmax(0.0, value);
            }
            setDampingInnerDistance(value);
        }

        std::string inputDampingOuterDistance;
        if(slotXML.hasOption("DampingOuterDistance")){
            inputDampingOuterDistance = slotXML.get("DampingOuterDistance");
        }else if(slotXML.hasOption("DumpingOuterDistance")){ // --> legacy check
            inputDampingOuterDistance = slotXML.get("DumpingOuterDistance");
        }
        if(!inputDampingOuterDistance.empty()){
            std::string input = bitpit::utils::string::trim(inputDampingOuterDistance);
            double value = 0.0;
            if(!input.empty()){
                std::stringstream ss(input);
                ss >> value;
                value = std::fmax(0.0, value);
            }
            setDampingOuterDistance(value);
        }
    }

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
template <std::size_t NCOMP>
void PropagateField<NCOMP>::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("Tolerance",std::to_string(m_tol));
    slotXML.set("UpdateThres",std::to_string(m_thres));
    slotXML.set("Print",std::to_string(int(m_print)));

    slotXML.set("NarrowBand", std::to_string(int(m_bandActive)));
    if(m_bandActive){
        slotXML.set("NarrowBandWidth",std::to_string(m_bandwidth));
        slotXML.set("NarrowBandRelax",std::to_string(m_bandrelax));
    }
    slotXML.set("Damping", std::to_string(int(m_dampingActive)));
    if(m_dampingActive){
        slotXML.set("DampingInnerDistance",std::to_string(m_plateau));
        slotXML.set("DampingOuterDistance",std::to_string(m_radius));
        slotXML.set("DampingType",std::to_string(m_dampingType));
        slotXML.set("DecayFactor",std::to_string(m_decayFactor));
    }
}

/*!
 * Restore data as in class default construction.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::clear(){
    BaseManipulation::clear();
    m_field.clear();
    m_dirichletPatches.clear();
    m_bc_dir.clear();
    m_dirichletBcs.clear();
    m_damping.clear();
    m_dampingSurfaces.clear();
    m_dampingUniSurface = nullptr;
    m_banddistances.clear();
    m_bandSurfaces.clear();
    m_bandUniSurface = nullptr;

    setDefaults();
}

/*!
 * Check coherence of the input data of the class, in particular:
   - check if points of Dirichlet patches belongs to the bulk mesh.
   - check if points of NarrowBand surfaces belongs to the bulk mesh.
   - check if points of Damping surfaces belongs to the bulk mesh.
 * Note. In parallel the bulk and boundaries meshes must have the same spatial partitioning.
 * \return true if coherence is satisfied, false otherwise.
 */
template <std::size_t NCOMP>
bool
PropagateField<NCOMP>::checkBoundariesCoherence(){

    if(!m_geometry) return false;

    std::unordered_set<long> nodeList;
    for(MimmoSharedPointer<MimmoObject> obj : m_dirichletPatches){
        std::vector<long> temp = obj->getVerticesIds(true); //only on internals
        nodeList.insert(temp.begin(), temp.end());
    }

    if(m_bandActive){
        for(MimmoSharedPointer<MimmoObject> obj : m_bandSurfaces){
            std::vector<long> temp = obj->getVerticesIds(true); //only on internals
            nodeList.insert(temp.begin(), temp.end());
        }
    }

    if(m_dampingActive){
        for(MimmoSharedPointer<MimmoObject> obj : m_dampingSurfaces){
            std::vector<long> temp = obj->getVerticesIds(true); //only on internals
            nodeList.insert(temp.begin(), temp.end());
        }
    }

    bitpit::PiercedVector<bitpit::Vertex> meshVertices = m_geometry->getVertices();
    for(long id : nodeList){
        if(!meshVertices.exists(id)){
            nodeList.clear();
            return false;
        }
    }
    return true;
}

/*!
 * Distribute the dirichlet bc values on the points of the bulk mesh.
   i.e. fill the internal member m_bc_dir.
 */
template <std::size_t NCOMP>
void
PropagateField<NCOMP>::distributeBCOnBoundaryPoints(){
    //estimate the total number of internal nodes carrying dirichlet bc
    int ncount = 0;
    for(MimmoSharedPointer<MimmoObject> obj : m_dirichletPatches){
        ncount += obj->getNInternalVertices();
    }
    m_bc_dir.clear();
    m_bc_dir.reserve(ncount);
    m_bc_dir.setDataLocation(MPVLocation::POINT);
    m_bc_dir.setGeometry(m_geometry);

    //starting from fields, fill in m_bc_dir and track those patches without a proper field;
    std::array<double,NCOMP> zeroval;
    zeroval.fill(0.0);
    std::unordered_set<MimmoSharedPointer<MimmoObject> > withoutField = m_dirichletPatches;
    for(MimmoPiercedVector<std::array<double,NCOMP>> * field : m_dirichletBcs){
        MimmoSharedPointer<MimmoObject> ref = field->getGeometry();
        //remove it from unvisited patch list
        withoutField.erase(ref);
        std::vector<long> ids = ref->getVerticesIds(true);
        for(long id : ids){
            if(m_bc_dir.exists(id)) continue;
            if(!field->exists(id)){
                m_bc_dir.insert(id, zeroval);
            }else{
                m_bc_dir.insert(id, field->at(id));
            }
        }
    }

    for(MimmoSharedPointer<MimmoObject> ref : withoutField){
        std::vector<long> ids = ref->getVerticesIds(true); //only internals
        for(long id : ids){
            if(m_bc_dir.exists(id)) continue;
            m_bc_dir.insert(id, zeroval);
        }
    }
    //shrink to fit to be sure.
    m_bc_dir.shrinkToFit();
}

/*!
 * Utility to instantiate a unique surface stitching a list of surface patches.
 * In both Serial and MPI version the final surface is provided as list
   of boundary patches referring to the target bulk mesh .
   The aim of this method is to stitch them together in a unique common surface.
   Moreover, in MPI versions, serialization is needed in order to provide all ranks
   to work with exactly the same surface (so no ghosts cells/vertices are retained).
   \param[in] listSurf list of boundary patches
   \param[out] uniSurf unique_ptr to reconstructed mesh.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::initializeUniqueSurface(const std::unordered_set<MimmoSharedPointer<MimmoObject> > & listSurf, MimmoSharedPointer<MimmoObject> & uniSurf){

    uniSurf = nullptr;
    // check if list is empty
    if (listSurf.empty()) {
        return;
    };

    //put all patches in the list in a unique surface (NOTE TYPE SURFACE). For MPI version, surface
    //are most likely partitioned, so you need to have information on internals/ghosts updated.
    //reconstructed surface is made by internals element/nodes only.
    MimmoSharedPointer<MimmoObject> tempSurface(new MimmoObject(1));
    for(MimmoSharedPointer<MimmoObject> obj : listSurf){
#if MIMMO_ENABLE_MPI
        if(!obj->isInfoSync()) obj->buildPatchInfo();
        if(!obj->arePointGhostExchangeInfoSync()) obj->updatePointGhostExchangeInfo();
#endif
        tempSurface->getPatch()->reserveVertices(tempSurface->getPatch()->getVertexCount() + obj->getPatch()->getVertexCount());
        for(bitpit::Vertex & vertex : obj->getVertices()){
            long vertexId = vertex.getId();
            tempSurface->addVertex(vertex, vertexId);
        }
        tempSurface->getPatch()->reserveCells(tempSurface->getPatch()->getCellCount() + obj->getPatch()->getCellCount());
        for (bitpit::Cell &cell : obj->getCells()){
            long idCell = cell.getId();
#if MIMMO_ENABLE_MPI
            int rank = m_rank;
            if (!cell.isInterior()){
                rank = obj->getPatch()->getCellRank(idCell);
            }
            tempSurface->addCell(cell, idCell, rank);
#else
            tempSurface->addCell(cell, idCell);
#endif
        }
    }

    tempSurface->getPatch()->squeezeCells();
    tempSurface->getPatch()->squeezeVertices();

#if MIMMO_ENABLE_MPI
    // Set partitioned patch to build parallel information
    tempSurface->setPartitioned();
#endif

    uniSurf = std::move(tempSurface);

    // Now uniSurf is the union of the surface patches in the list.
    // The original connections are untouched and no new connections are imposed
    // Duplicated vertices from different surface patches are maintained

}

/*!
   It computes for the first time the damping function (artificial diffusivity)
   used to modulate Laplacian solution over the target mesh.
   The damping function is a scalar field data referred on CELLS of the bulk
  (ghost included) and it will be stored in internal member m_damping.
   Please be sure to create the unique Surface m_dampingUniSurface first.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::initializeDampingFunction(){

    bitpit::PatchKernel * patch_ = getGeometry()->getPatch();
    m_damping.clear();
    m_damping.reserve(getGeometry()->getNCells());
    m_damping.setGeometry(getGeometry());
    m_damping.setDataLocation(MPVLocation::CELL);

    for (auto it = patch_->cellBegin(); it!=patch_->cellEnd(); ++it){
        m_damping.insert(it.getId(), 1.0);
    }
    updateDampingFunction();
}

/*!
 * It updates an existent damping function (artificial diffusivity) stored in m_damping member.
 * It will extract the pool of cells with diffusivity > 1.0. Thus, it will update
 * the narrow band inside the radius of influence, and recalculate distance and volumes eventually.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::updateDampingFunction(){

    if(!m_dampingUniSurface.get()){
        //this unique_ptr structure is initialized with initializeUniqueSurface applied to damping surfaces
        // if damping is not active is set to nullptr
        // leaving damping as it is.
        return;
    }

    const double maxd(m_radius);
    //get the list of elements in m_damping with diffusivity > 1.0;
    // at the same time reset the dumping function values to 1.0;
    livector1D seedlist;
    seedlist.reserve(m_damping.size());
    for(auto it= m_damping.begin(); it!=m_damping.end(); ++it){
        if(*it > 1.0){
            seedlist.push_back(it.getId());
            //reset local value of m_damping to 1.0.
            *it = 1.0;
        }
    }

    // if seedlist is still empty have a try finding some seeds along
    // bulk volume mesh borders (ghost included)
    if (seedlist.empty()){
        // Initialize seeds to avoid use of skdtree in narrowband computing
        seedlist = getGeometry()->extractBoundaryCellID(true); // with ghost included
    }
    seedlist.shrink_to_fit();

    //reevaluate narrow band cells at distance d < maxd
    bitpit::PiercedVector<double> distFactor = getGeometry()->getCellsNarrowBandToExtSurfaceWDist(*(m_dampingUniSurface.get()), maxd, &seedlist);
    seedlist.clear();

    double distanceMax = std::pow((maxd/m_plateau), m_decayFactor);
    for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
        if(*it < m_plateau){
            (*it) = 1.0;
        }else{
            (*it) = (std::pow(maxd/(*it), m_decayFactor) -1.0) / (distanceMax -1.0);
        }
    }

    if(m_dampingType != 1) { // no volume stuff, fill dumping with distance part only
        for(auto it=distFactor.begin(); it!=distFactor.end(); ++it){
            m_damping.at(it.getId()) = (distanceMax - 1.0)*(*it) + 1.0;;
        }
    }else{
        // Evaluating the volume part.
        bitpit::PiercedVectorStorage<double> volFactor;
        volFactor.setStaticKernel(&distFactor); //, bitpit::PiercedSyncMaster::SyncMode::SYNC_MODE_DISABLED);
        volFactor.fill(1.0);

        //evaluating cell volumes
        double locvol;
        std::size_t countNegativeVolumes(0);

        double volmax = 0.0, volmin=std::numeric_limits<double>::max();
        for(auto it=distFactor.begin(); it!=distFactor.end(); ++it){
            locvol = getGeometry()->evalCellVolume(it.getId());
            if(locvol < std::numeric_limits<double>::min()){
                ++countNegativeVolumes;
                locvol = std::numeric_limits<double>::min(); //to assess myself around a 1.E-38 as minimum.
            }
            volFactor.rawAt(it.getRawIndex()) = locvol;
            volmin = std::min(volmin,locvol);
            volmax = std::max(volmax,locvol);
        }
        if(countNegativeVolumes > 0){
            m_log->setPriority(bitpit::log::Verbosity::DEBUG);
            (*m_log)<<"Warning in "<<m_name<<". Detected " << countNegativeVolumes<<" cells with almost zero or negative volume"<<std::endl;
            m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        }

#if MIMMO_ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &volmin, 1, MPI_DOUBLE, MPI_MIN, m_communicator);
        MPI_Allreduce(MPI_IN_PLACE, &volmax, 1, MPI_DOUBLE, MPI_MAX, m_communicator);
        MPI_Allreduce(MPI_IN_PLACE, &countNegativeVolumes, 1, MPI_INT, MPI_SUM, m_communicator);
#endif
        //evaluate the volume normalized function and store it in dumping.
        for(auto it = distFactor.begin(); it !=distFactor.end(); ++it){
            m_damping.at(it.getId()) = std::pow(1.0 + (volmax -volmin)/volFactor.rawAt(it.getRawIndex()), *it);
        }
    }
	//all done if you are here, you computed also the bulk part and put together all the stuffs.
    // you don't need to communicate ghost data for dumping. Everyone, ghost included had the correct info.
}

/*!
    Compute/Update the list of vertices in the narrow band and store it with their distance
    in m_banddistances member.
    If m_banddistances is not empty, use its information as starting point to update
    the narrow band list.
    It requires the initialization of a Unique Surface (m_bandUniSurface) for m_bandSurfaces list first.
*/
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::updateNarrowBand(){

    if(!m_bandUniSurface.get()){
        //this unique_ptr structure is initialized with initializeUniqueSurface applied to narrow band surfaces
        // if narrow band control is not active is set to nullptr
        // and narrow band is not computed.
        m_banddistances.clear();
        return;
    }

    //get the list of vertex elements in m_banddistances;
    livector1D seedlist = m_banddistances.getIds();

    // is seedlist is still empty have a try finding some vertex seeds along
    // bulk volume mesh borders (ghost included)
    if (seedlist.empty()){
        // Initialize seeds to avoid use of skdtree in narrowband computing
        seedlist = getGeometry()->extractBoundaryVertexID(true); // ghost included
    }
    seedlist.shrink_to_fit();

    //re-evaluate narrow band vertices at distance d < m_bandwidth
    m_banddistances = getGeometry()->getVerticesNarrowBandToExtSurfaceWDist(*(m_bandUniSurface.get()), m_bandwidth, &seedlist);
}


/*!
 * The method modifies raw Laplacian stencils (without bc already assigned) on
   points in Graph Laplacian approximation to take in account a relaxed
   solution in the narrow band near some target boundary surfaces.
   The idea is to manipulate stencils of a mesh point inside the narrow band so that
   nodes nearer to reference boundaries "weights" more than farther ones.
   This upwinding effect let the solution to diffuse more inside the bulk core.
   Outside the narrow band stencils are the usual ones.
 * \param[in, out] laplaceStencils unique pointer to original set of laplacian stencils in input, narrowband modified in output.
 */
 template<std::size_t NCOMP>
 void
 PropagateField<NCOMP>::modifyStencilsForNarrowBand(GraphLaplStencil::MPVStencilUPtr &laplaceStencils ){

    livector1D nbv_list = m_banddistances.getIds();
    //loop on stencils associated to nbv_list
   for (long idT : nbv_list){
       //if belongs to laplace stencils the point is not a ghost for sure.
       // So discard point not included in the laplaceStencils
       if(!laplaceStencils->exists(idT)) continue;

        double localdist = m_banddistances.at(idT);
        bitpit::StencilScalar &loc_stencil = laplaceStencils->at(idT);
        std::size_t size = loc_stencil.size();
        //check content of the stencil.
        double correction_old = 0.0;
        double correction_new = 0.0;
        //the last one is the point idT itself carrying the complement weight.
        for (std::size_t i=0; i<size-1; ++i){
            long &idL = loc_stencil.patternData()[i];
            if(!m_banddistances.exists(idL))         continue;
            if(m_banddistances.at(idL) > localdist ){
                correction_old += loc_stencil.weightData()[i];
                loc_stencil.weightData()[i] *= m_bandrelax; // use here the band relaxation factor.
                correction_new += loc_stencil.weightData()[i];
            }

        }
        loc_stencil.weightData()[size-1] += (correction_old - correction_new); //(-sum(w_i) + sum(w_excluded) - sum_recalc );
    }
}

/*!
 * Prepare your system solver, feeding the laplacian stencils you previosly calculated
 * with GraphLaplStencil::computeLaplacianStencil method of StencilFunctions.
 * Provide the map that get consecutive Index from Global Pierced vector Index system for POINTS
 * The stencil will be renumerated with the consecutiveIdIndexing provided.
 * The method requires the m_solver to be instantiated already
 *
 * param[in] laplacianStencils pointer to MPV structure of laplacian stencils.
 * param[in] map of consecutive points ID from Global PV indexing (typically get from MimmoObject::getMapDataInv)
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::initializeLaplaceSolver(GraphLaplStencil::MPVStencil * laplacianStencils, const lilimap & maplocals){

	bitpit::KSPOptions &solverOptions = m_solver->getKSPOptions();

	solverOptions.rtol      = m_tol;
	solverOptions.subrtol   = m_tol;
	// total number of local DOFS, determines size of matrix
	long nDOFs = laplacianStencils->size();

	// total number of non-zero elements in the stencils.
	long nNZ(0);
	for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
		nNZ += it->size();
	}

#if MIMMO_ENABLE_MPI==1
	//instantiate the SparseMatrix
	bitpit::SparseMatrix matrix(m_communicator, getGeometry()->getPatch()->isPartitioned(), nDOFs, nDOFs, nNZ);
#else
	//instantiate the SparseMatrix
	bitpit::SparseMatrix matrix(nDOFs, nDOFs, nNZ);
#endif

	std::vector<long> mapsort(laplacianStencils->size());
	long id, ind;
	for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
		id = it.getId();
		ind = maplocals.at(id);
#if MIMMO_ENABLE_MPI
		ind -= getGeometry()->getPointGlobalCountOffset();
#endif
		mapsort[ind] = id;
	}

	//Add ordered rows
	for(auto id : mapsort){
		bitpit::StencilScalar item = laplacianStencils->at(id);
		//renumber values on the fly.
		item.renumber(maplocals);
		matrix.addRow(item.size(), item.patternData(), item.weightData());
	}

	//assembly the matrix;
	matrix.assembly();

	//clean up the previous stuff in the solver.
	m_solver->clear();
	// now you can initialize the m_solver with this matrix.
	m_solver->getKSPOptions().restart = 30;
	m_solver->getKSPOptions().overlap = 1;
	m_solver->getKSPOptions().sublevels = 1;
	m_solver->assembly(matrix);
}

/*!
 * Update your system solver, feeding the cell based laplacian stencils you want to update in the matrix.
 * This method works with any valid subset of stencils in the mesh, but require the solver matrix to be initialized
 * and to have the new stencils with the same id pattern as they had at the time of the matrix initialization.
 * Provide the map that get consecutive Index from Global Pierced vector Index system for CELLS
 * The stencil will be renumerated with the consecutiveIdIndexing provided.
 *
 * param[in] laplacianStencils pointer to MPV structure of laplacian stencils subset to feed as update.
 * param[in] map of consecutive cell ID from Global PV indexing (typically get from MimmoObject::getMapcellInv)
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::updateLaplaceSolver(FVolStencil::MPVDivergence * laplacianStencils, const lilimap & maplocals){

	// total number of local DOFS, determines size of matrix
	long nDOFs = m_solver->getColCount();
	long nupdate = laplacianStencils->size();

	// total number of non-zero elements in the stencils.
	long nNZ(0);
	for(auto it=laplacianStencils->begin(); it!=laplacianStencils->end(); ++it){
		nNZ += it->size();
	}

#if MIMMO_ENABLE_MPI==1
	//instantiate the SparseMatrix
	bitpit::SparseMatrix upelements(m_communicator, getGeometry()->getPatch()->isPartitioned(), nupdate, nDOFs, nNZ);
#else
	//instantiate the SparseMatrix
	bitpit::SparseMatrix upelements(nupdate, nDOFs, nNZ);
#endif


	// store the local ind of the rows involved, while filling the matrix of update values.
	std::vector<long> rows_involved;
	rows_involved.reserve(laplacianStencils->size());

	long id, ind;
	for(auto it=laplacianStencils->begin(); it != laplacianStencils->end(); ++it){
		id = it.getId();
		ind = maplocals.at(id);
#if MIMMO_ENABLE_MPI
		ind -= getGeometry()->getPointGlobalCountOffset();
#endif
		rows_involved.push_back(ind);
		bitpit::StencilScalar item(*it);
		item.renumber(maplocals);
		upelements.addRow(item.size(), item.patternData(), item.weightData());
	}
	//assembly the update matrix;
	upelements.assembly();

	//call the solver update;
	m_solver->update(rows_involved.size(), rows_involved.data(), upelements);

}

/*!
 * This method evaluate the bc corrections for a singular run of the system solver,
 * update the system matrix in m_solver and evaluate the rhs part due to bc.
 * After you call this method, you are typically ready to solve the laplacian system.
 * bc @ nodes are directly taken from class member m_bc_dir.
 * The method requires the Laplacian m_solver to be initialized. No ghost are taken into account.
 *
 * Moreover it needs as input :
 * \param[in] comp target component of bc conditions.
 * \param[in] unused boolean
 * \param[in] borderLaplacianStencil list of laplacian Stencil on border nodes, where the bc is temporarily imposed as homogeneous Neumann
 * \param[in] maplocals map from global id numbering to local system solver numbering.
 * \param[in,out] rhs vector of right-hand-sides to append constant data from bc corrections.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::assignBCAndEvaluateRHS(std::size_t comp, bool unused,
                                              GraphLaplStencil::MPVStencil * borderLaplacianStencil,
                                              const lilimap & maplocals, dvector1D & rhs)
{
    BITPIT_UNUSED(unused);
    //resize rhs to the number of internal cells
    MimmoSharedPointer<MimmoObject> geo = getGeometry();
    rhs.resize(geo->getNInternalVertices(), 0.0);

    if (!m_solver->isAssembled()) {
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Warning in "<<m_name<<". Unable to assign BC to the system. The solver is not yet initialized."<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        return;
    }

    if (!borderLaplacianStencil) {
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Warning in "<<m_name<<". Unable to reach border nodes stencils data. Nothing to do."<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        return;
    }

    //copy laplacian stencils in a work mpv .
    GraphLaplStencil::MPVStencilUPtr lapwork(new GraphLaplStencil::MPVStencil(*borderLaplacianStencil));

    //correct the original border laplacian stencils applying the Dirichlet conditions.
    //Neumann are implicitely imposed by graph-laplacian scheme.
    //renumber it and update the laplacian matrix and fill the rhs.
    bitpit::StencilScalar correction;

    //loop on all dirichlet boundary nodes.
    for(long id : m_bc_dir.getIds()){
        if (geo->isPointInterior(id)){
            //apply the correction relative to bc @ dirichlet node.
            correction.clear(true);
            correction.appendItem(id, 1.);
            correction.sumConstant(-m_bc_dir[id][comp]);
            //Fix to zero the old stencil (the update of system solver doesn't substitute but modify or append new pattern item and weights)
            lapwork->at(id) *= 0.;
            lapwork->at(id) += correction;
        }
    }

    // now its time to update the solver matrix and to extract the rhs contributes.
    updateLaplaceSolver(lapwork.get(), maplocals);

    // now get the rhs
    for(auto it = lapwork->begin(); it != lapwork->end();++it){
        auto index = maplocals.at(it.getId());
#if MIMMO_ENABLE_MPI
        //correct index if in parallel
        index -= geo->getPointGlobalCountOffset();
#endif
        rhs[index] -= it->getConstant();
    }
}


/*!
 * It solves the laplacian problem, one field at a time.
 * Basically it is a wrapper to m_solver method solve();
 * Before calling this method be sure to have initialized the Laplacian linear system
 *
 * \param[in] rhs on internal cells;
 * \param[in,out] result on internal cells. In input is an initial solution, in output is the reusl of computation.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::solveLaplace(const dvector1D &rhs, dvector1D &result){

    result.resize(getGeometry()->getNInternalVertices(), 0.);

    // Check if the internal solver is initialized
    if (!m_solver->isAssembled()) {
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Warning in "<<m_name<<". Unable to solve the system. The solver is not yet initialized."<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        return;
    }

    // Solve the system
    m_solver->solve(rhs, &result);

    //think I've done my job.
}

/*!
 * Utility to put laplacian solution directly into m_field (cleared and refreshed).
 * Ghost communication is already taken into account in case of MPI version.
 * \param[in] results data of laplacian solutions collected in raw vectors.
 * \param[in] mapglobals map to retrieve cell/node global id from locals raw laplacian indexing
 * \param[out] marked (OPTIONAL) list of cells/nodes whose field norm is greater then m_thres value.
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::reconstructResults(const dvector2D & results, const lilimap & mapglobals, livector1D * marked)
{
    if(results.size() != NCOMP){
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log) << "WARNING in "<<m_name<<" . A field with dimension different from" <<NCOMP<<" is feeded to reconstructResults. m_field is not touched"<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);
        return;
    }
    // push result in a mpv linked to target mesh and on cell location.
    MimmoSharedPointer<MimmoObject> geo = getGeometry();
    std::unique_ptr<MimmoPiercedVector<std::array<double,NCOMP> > > mpvres(new MimmoPiercedVector<std::array<double,NCOMP> >(geo, MPVLocation::POINT));
    mpvres->reserve(geo->getNVertices());

    long id,counter;
    std::array<double, NCOMP> temp;
    for(auto & pair : mapglobals){
        id = pair.second;
        if (geo->isPointInterior(id)){
            counter = pair.first;
#if MIMMO_ENABLE_MPI
            counter -= geo->getPointGlobalCountOffset();
#endif
            for(int i=0; i<int(NCOMP); ++i){
                temp[i] = results[i][counter];
            }
            mpvres->insert(id, temp);
        }else{
            for(int i=0; i<int(NCOMP); ++i){
                temp[i] = 0.0;
            }
            mpvres->insert(id, temp);
        }
    }

#if MIMMO_ENABLE_MPI
    //communicate ghosts
    communicatePointGhostData(mpvres.get());
#endif

    //mark point with solution norm above m_thres
    if(marked){
        marked->clear();
        marked->reserve(mpvres->size());
        for(auto it=mpvres->begin(); it != mpvres->end(); ++it){
            if(norm2(*it) > m_thres){
                marked->push_back(it.getId());
            }
        }
        marked->shrink_to_fit();
    }

    // store all in m_field.
    m_field.swap(*(mpvres.get()));
}

//******************
//EXCLUSIVE MPI METHODS
//******************
#if MIMMO_ENABLE_MPI
/*!
    Creates a new ghost communicator and return its tag. if already exists, do nothing
    and return its current tag.

    \param[in] refGeo pointer to reference partitioned MimmoObject
    \param[in] continuous defines if the communicator will be set in continuous mode
    \return The tag associated to the newly created/or already existent communicator.
 */
template<std::size_t NCOMP>
int
PropagateField<NCOMP>::createGhostCommunicator(MimmoObject* refGeo, bool continuous){

    if(m_ghostCommunicators.count(refGeo) == 0){
        // Create communicator
        m_ghostCommunicators[refGeo] = std::unique_ptr<GhostCommunicator>(new GhostCommunicator(refGeo->getPatch()));
        m_ghostCommunicators[refGeo]->resetExchangeLists();
        m_ghostCommunicators[refGeo]->setRecvsContinuous(continuous);
    }
    // Return Communicator tag
    return int(m_ghostCommunicators[refGeo]->getTag());
}

/*!
    Creates a new point ghost communicator.

    \param[in] refGeo pointer to reference partitioned MimmoObject
    \param[in] continuous defines if the communicator will be set in continuous mode
    \return The tag associated to the newly created communicator.
 */
template<std::size_t NCOMP>
int
PropagateField<NCOMP>::createPointGhostCommunicator(MimmoObject* refGeo, bool continuous){
    if(m_pointGhostCommunicators.count(refGeo) == 0){
        // Create communicator
        m_pointGhostCommunicators[refGeo] = std::unique_ptr<PointGhostCommunicator>(new PointGhostCommunicator(refGeo));
        m_pointGhostCommunicators[refGeo]->resetExchangeLists();
        m_pointGhostCommunicators[refGeo]->setRecvsContinuous(continuous);
    }
    // Return Communicator tag
    return int(m_pointGhostCommunicators[refGeo]->getTag());
}

/*!
    Communicate MPV data on ghost cells on the reference geometry linked by data itself.
    The method creates a new communicator and streamer if not already allocated.
    Otherwise the pointer of the communicator to the data is updated with the input argument.
    \param[in] data Pointer to field with data to communicate
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::communicateGhostData(MimmoPiercedVector<std::array<double, NCOMP> > *data){
    // Creating cell ghost communications for exchanging interpolated values
    MimmoObject * geo = data->getGeometry().get();
    if(!geo){
        throw std::runtime_error("Propagate Class ::communicateGhostData no ref Geometry in mpv data!");
    }
    //if geo is not partitioned you have nothing to communicate.
    if (!geo->isPartitioned()) return;

    //check for communicator on geometry, if not exists create it
    m_ghostTags[geo] = createGhostCommunicator(geo, true);
    //after this call a communicator dedicated to geo surely exists
    //check if streamer for data type exists
    if(m_ghostStreamers.count(geo) > 0){
        //set data to the streamer. If you created the streamer
        //you have already add it to its comunicator.
        m_ghostStreamers[geo]->setData(data);
    }else{
        //you need to create it brand new. Attach data directly.
        m_ghostStreamers[geo] = std::unique_ptr<MimmoDataBufferStreamer<NCOMP>>(new MimmoDataBufferStreamer<NCOMP>(data));
        m_ghostCommunicators[geo]->addData(m_ghostStreamers[geo].get());
    }

    // Send data
    m_ghostCommunicators[geo]->startAllExchanges();
    // Receive data
    m_ghostCommunicators[geo]->completeAllExchanges();
}

/*!
    Communicate MPV data on ghost nodes on the reference geometry linked by data itself.
    The method creates a new communicator and streamer if not already allocated.
    Otherwise the pointer of the communicator to the data is updated with the input argument.
    \param[in] data Pointer to field with data to communicate
 */
template<std::size_t NCOMP>
void
PropagateField<NCOMP>::communicatePointGhostData(MimmoPiercedVector<std::array<double, NCOMP> > *data){
    // Creating cell ghost communications for exchanging interpolated values
    MimmoObject * geo = data->getGeometry().get();
    if(!geo){
        throw std::runtime_error("Propagate Class ::communicatePointGhostData no ref Geometry in mpv data!");
    }
    //if geo is not partitioned you have nothing to communicate.
    if (!geo->isPartitioned()) return;

    //check for communicator on geometry, if not exists create it
    m_pointGhostTags[geo] = createPointGhostCommunicator(geo, true);
    //after this call a communicator dedicated to geo surely exists
    //check if streamer for data type exists
    if(m_pointGhostStreamers.count(geo) > 0){
        //set data to the streamer. If you created the streamer
        //you have already add it to its comunicator.
        m_pointGhostStreamers[geo]->setData(data);
    }else{
        //you need to create it brand new. Attach data directly.
        m_pointGhostStreamers[geo] = std::unique_ptr<MimmoPointDataBufferStreamer<NCOMP>>(new MimmoPointDataBufferStreamer<NCOMP>(data));
        m_pointGhostCommunicators[geo]->addData(m_pointGhostStreamers[geo].get());
    }

    // Send data
    m_pointGhostCommunicators[geo]->startAllExchanges();
    // Receive data
    m_pointGhostCommunicators[geo]->completeAllExchanges();
}

#endif


} //end of mimmo namespace
