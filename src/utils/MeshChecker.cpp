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

#include <MeshChecker.hpp>
#include <bitpit_common.hpp>

namespace mimmo{

/*!
 * Constructor
 */
MeshChecker::MeshChecker()
{
	m_name = "mimmo.MeshChecker";
	setDefault();
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
MeshChecker::MeshChecker(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.MeshChecker";
	setDefault();

	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);
	if(input == "mimmo.MeshChecker"){
		absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*! Destructor
 *
 */
MeshChecker::~MeshChecker(){};

/*! Copy Constructor. Result displacements are never copied.
 *\param[in] other MeshChecker where copy from
 */
MeshChecker::MeshChecker(const MeshChecker & other):BaseManipulation(other){
	m_minVolume = other.m_minVolume;
	m_maxVolume = other.m_maxVolume;
	m_maxSkewness = other.m_maxSkewness;
	m_maxSkewnessBoundary = other.m_maxSkewnessBoundary;
	m_minFaceValidity = other.m_minFaceValidity;
	m_maxFaceValidity = other.m_maxFaceValidity;
	m_minVolumeChange = other.m_minVolumeChange;
	m_minVolumeTol = other.m_minVolumeTol;
	m_maxVolumeTol = other.m_maxVolumeTol;
	m_maxSkewnessTol = other.m_maxSkewnessTol;
	m_maxSkewnessBoundaryTol = other.m_maxSkewnessBoundaryTol;
	m_minFaceValidityTol = other.m_minFaceValidityTol;
	m_minVolumeChangeTol = other.m_minVolumeChangeTol;
	m_isGood = other.m_isGood;
    m_qualityStatus = other.m_qualityStatus;
    m_printResumeFile = other.m_printResumeFile;
};

/*!
 * Assignement operator.  Result displacements are never copied.
 */
MeshChecker & MeshChecker::operator=(MeshChecker other){
	swap(other);
	return *this;
}

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void
MeshChecker::swap(MeshChecker & x) noexcept
{
	std::swap(m_minVolume, x.m_minVolume);
	std::swap(m_maxVolume, x.m_maxVolume);
	std::swap(m_maxSkewness, x.m_maxSkewness);
	std::swap(m_maxSkewnessBoundary, x.m_maxSkewnessBoundary);
	std::swap(m_minFaceValidity, x.m_minFaceValidity);
	std::swap(m_maxFaceValidity, x.m_maxFaceValidity);
	std::swap(m_minVolumeChange, x.m_minVolumeChange);
	std::swap(m_minVolumeTol, x.m_minVolumeTol);
	std::swap(m_maxVolumeTol, x.m_maxVolumeTol);
	std::swap(m_maxSkewnessTol, x.m_maxSkewnessTol);
	std::swap(m_maxSkewnessBoundaryTol, x.m_maxSkewnessBoundaryTol);
	std::swap(m_minFaceValidityTol, x.m_minFaceValidityTol);
	std::swap(m_minVolumeChangeTol, x.m_minVolumeChangeTol);
	std::swap(m_isGood, x.m_isGood);
    std::swap(m_qualityStatus, x.m_qualityStatus);
    std::swap(m_printResumeFile, x.m_printResumeFile);
}

/*!
 * Set Default Value
 */
void
MeshChecker::setDefault()
{
	m_minVolume = 1.e+18;
	m_maxVolume = 0.;
	m_maxSkewness = 0.;
	m_maxSkewnessBoundary = 0.;
	m_minFaceValidity = 1.;
	m_maxFaceValidity = 0.;
	m_minVolumeChange = 1.;
	m_minVolumeTol = 1.0e-14;
	m_maxVolumeTol = 1.0e+10;
	m_maxSkewnessTol = 85. * BITPIT_PI / 180.;
	m_maxSkewnessBoundaryTol = m_maxSkewnessTol;
	m_minFaceValidityTol = 0.5;
	m_minVolumeChangeTol = 1.0e-05;
	m_isGood = false;
    m_qualityStatus = CMeshOutput::NOTRUN;
    m_printResumeFile = false;
	m_volumes.clear();
}

/*!
 * Set Default Value
 */
void
MeshChecker::execute()
{
	if (getGeometry() == nullptr){
		throw std::runtime_error (m_name + " : nullptr pointer to linked geometry");
	}
	if (getGeometry()->isEmpty()){
        (*m_log)<<"WARNING: "<<m_name<< " : empty geometry found in execution"<<std::endl;
	}

    m_volume = std::unique_ptr<MimmoObject>(new MimmoObject(getGeometry()->getType()));
    m_volumechange = std::unique_ptr<MimmoObject>(new MimmoObject(getGeometry()->getType()));
    m_skewness = std::unique_ptr<MimmoObject>(new MimmoObject(getGeometry()->getType()));
    m_facevalidity =std::unique_ptr<MimmoObject>(new MimmoObject(getGeometry()->getType()));

    m_qualityStatus = checkMeshQuality();
	m_isGood = ( m_qualityStatus == CMeshOutput::GOOD );

    if(m_printResumeFile){
        printResumeFile();
    }
	//Clear auxiliary variables
	clear();

}

/*! It builds the input/output ports of the object
 */
void
MeshChecker::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, MeshChecker>(this, &MeshChecker::setGeometry, M_GEOM, true));
    built = (built && createPortOut<bool, MeshChecker>(this, &mimmo::MeshChecker::isGood, M_VALUEB));
    built = (built && createPortOut<int, MeshChecker>(this, &mimmo::MeshChecker::getQualityStatusInt, M_VALUEI));
	m_arePortsBuilt = built;
};

/*!
    Overload BaseManipulation setGeometry
    \param[in] geo target geometry
*/
void MeshChecker::setGeometry(MimmoSharedPointer<MimmoObject> geo){
    if(geo){
        BaseManipulation::setGeometry(geo);
        //Clear auxiliary variables
        clear();
    }
}

/*!
 * Set the tolerance
 */
void
MeshChecker::setMinimumVolumeTolerance(double tol)
{
	m_minVolumeTol = tol;
}

/*!
 * Set the tolerance
 */
void
MeshChecker::setMaximumVolumeTolerance(double tol)
{
	m_maxVolumeTol = tol;
}

/*!
 * Set the tolerance in degrees
 */
void
MeshChecker::setMaximumSkewnessTolerance(double tol)
{
	m_maxSkewnessTol = tol / 180. * BITPIT_PI;
}

/*!
 * Set the tolerance in degrees
 */
void
MeshChecker::setMaximumBoundarySkewnessTolerance(double tol)
{
	m_maxSkewnessBoundaryTol = tol / 180. * BITPIT_PI;
}

/*!
 * Set the tolerance
 */
void
MeshChecker::setMinimumFaceValidityTolerance(double tol)
{
	m_minFaceValidityTol = tol;
}

/*!
 * Set the tolerance
 */
void
MeshChecker::setMinimumVolumeChangeTolerance(double tol)
{
	m_minVolumeChangeTol = tol;
}

/*!
 * set true to print resume file
 */
void
MeshChecker::setPrintResumeFile(bool flag)
{
	m_printResumeFile = flag;
}


/*!
 * \return true is quality is good
 */
bool
MeshChecker::isGood()
{
	return int(m_isGood);
}

/*!
 * \return quality flag (see class doc)
 */
MeshChecker::CMeshOutput
MeshChecker::getQualityStatus()
{
	return m_qualityStatus;
}

/*!
 * \return quality flag (see class doc)
 */
int
MeshChecker::getQualityStatusInt()
{
	return static_cast<int>(m_qualityStatus);
}


/*! It check the quality of the mesh
 * \return error flag. If mesh is good return 0. Otherwise the error with highest priority. See enum CMeshOutput for errors and their priority.
 */
MeshChecker::CMeshOutput
MeshChecker::checkMeshQuality()
{
    CMeshOutput check = CMeshOutput::NOTRUN;

    (*m_log) << bitpit::log::priority(bitpit::log::NORMAL);
    (*m_log) << bitpit::log::context("mimmo");

    try{
        checkVolume();
        checkSkewness();
    	checkFaceValidity();
    }catch(std::exception & e){
        (*m_log)<<m_name<<" : FAILED to calculate quality indicators. Check NOT RUN"<< std::endl;
        return check;
    }

    check = CMeshOutput::GOOD;
#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif
	if (m_maxVolume > m_maxVolumeTol){
        check = CMeshOutput::MAXIMUMVOLUME;
		(*m_log)<<m_name<<" : FAILED maximum cell volume check, found value : "<< m_maxVolume << std::endl;
	}
	else{
		(*m_log)<<m_name<<" : maximum cell volume check PASSED with value : "<< m_maxVolume << std::endl;
	}

#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif
	if (m_minVolume < m_minVolumeTol){
        check = CMeshOutput::MINIMUMVOLUME;
		(*m_log)<<m_name<<" : FAILED minimum cell volume check, found value : "<< m_minVolume << std::endl;
	}
	else{
		(*m_log)<<m_name<<" : minimum cell volume check PASSED with value : "<< m_minVolume << std::endl;
	}

#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif
	if (m_maxSkewness > m_maxSkewnessTol){
        check = CMeshOutput::SKEWNESS;
		(*m_log)<<m_name<<" : FAILED cell skewness, found value [deg] : "<< m_maxSkewness / BITPIT_PI * 180. << std::endl;
	}
	else{
		(*m_log)<<m_name<<" : cell skewness PASSED with value [deg] : "<< m_maxSkewness / BITPIT_PI * 180. << std::endl;
	}

#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif
	if (m_maxSkewnessBoundary > m_maxSkewnessBoundaryTol){
        check = CMeshOutput::BOUNDARYSKEWNESS;
		(*m_log)<<m_name<<" : FAILED boundary skewness, found value [deg] : "<< m_maxSkewnessBoundary / BITPIT_PI * 180. << std::endl;
	}
	else{
		(*m_log)<<m_name<<" : boundary skewness PASSED with value [deg] : "<< m_maxSkewnessBoundary / BITPIT_PI * 180. << std::endl;
	}

#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif
	if (m_minVolumeChange < m_minVolumeChangeTol){
        check = CMeshOutput::VOLUMECHANGERATIO;
		(*m_log)<<m_name<<" : FAILED cell volume change check, found value : "<< m_minVolumeChange << std::endl;
	}
	else{
		(*m_log)<<m_name<<" : cell volume change check PASSED with value : "<< m_minVolumeChange << std::endl;
	}

#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif
	if (m_minFaceValidity < m_minFaceValidityTol){
		check = CMeshOutput::FACEVALIDITY;
		(*m_log)<<m_name<<" : FAILED cell face validity check, found value : "<< m_minFaceValidity << std::endl;
	}
	else{
		(*m_log)<<m_name<<" : cell face validity check PASSED with value : "<< m_minFaceValidity << std::endl;
	}

#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif
	if (check == CMeshOutput::GOOD){
		(*m_log)<< m_name << " : ALL checks PASSED " << std::endl;
	}

    (*m_log) << bitpit::log::priority(bitpit::log::DEBUG);


	return check;
}

/*! It check the quality of the volume cells of the mesh
 * \return false if one of the volumes tolerance is not satisfied
 */
bool
MeshChecker::checkVolume()
{

	m_volumes.clear();
    initializeVolumes();
    getGeometry()->updateAdjacencies();

    livector1D mvolcells;
    livector1D mvolchangecells;

	for (auto it = getGeometry()->getPatch()->internalBegin(); it!= getGeometry()->getPatch()->internalEnd(); ++it){
		long id = it.getId();
		double vol = m_volumes[id];

		m_minVolume = std::min(m_minVolume, vol);
		if (vol < m_minVolumeTol){
			mvolcells.push_back(id);
		}

		m_maxVolume = std::max(m_maxVolume, vol);
		if (vol > m_maxVolumeTol){
			mvolcells.push_back(id);
		}

		double ratio = 1.;
		int nneigh = it->getAdjacencyCount();
		long* neighs = it->getAdjacencies();
		for (int i=0; i<nneigh; i++){
			long neigh = neighs[i];
			if (neigh > -1)
				ratio = std::min(ratio, vol/m_volumes[neigh]);
		}
		m_minVolumeChange = std::min(m_minVolumeChange, ratio);

        if (ratio < m_minVolumeChangeTol){
			mvolchangecells.push_back(id);
		}
	}

    bool passed = ( mvolcells.empty() && mvolchangecells.empty() );
	if (!passed){

        bitpit::PiercedVector<bitpit::Vertex> &vertices = getGeometry()->getVertices();
        bitpit::PiercedVector<bitpit::Cell> &cells    = getGeometry()->getCells();
        bitpit::PiercedVector<bitpit::Vertex> &localVolumeVerts = m_volume->getVertices();
        bitpit::PiercedVector<bitpit::Vertex> &localVolumeChangeVerts = m_volumechange->getVertices();

        for(long idCell : mvolcells){
            bitpit::Cell & cell = cells.at(idCell);
            bitpit::ConstProxyVector<long> connIds = cell.getVertexIds();

            for(long idV : connIds){
                if(!localVolumeVerts.exists(idV)){
                    m_volume->addVertex(vertices.at(idV), idV);
                }
            }
            m_volume->addCell(cell, idCell);
        }

        for(long idCell : mvolchangecells){
            bitpit::Cell & cell = cells.at(idCell);
            bitpit::ConstProxyVector<long> connIds = cell.getVertexIds();

            for(long idV : connIds){
                if(!localVolumeChangeVerts.exists(idV)){
                    m_volumechange->addVertex(vertices.at(idV), idV);
                }
            }
            m_volumechange->addCell(cell, idCell);
        }
	}

#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &passed, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &m_minVolumeChange, 1, MPI_DOUBLE, MPI_MIN, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &m_maxVolume, 1, MPI_DOUBLE, MPI_MAX, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &m_minVolume, 1, MPI_DOUBLE, MPI_MIN, m_communicator);

#endif

	return passed;
}

/*! It check the quality of the skewness of the cells of the mesh
 * \return false if skewness or boundary skewness tolerance is not satisfied
 */
bool
MeshChecker::checkSkewness()
{
	if (getGeometry()->getType() != 2){
		(*m_log)<<m_name<<" : skewness check active only for volume mesh -> check disabled" << std::endl;
		return true;
	}

    getGeometry()->updateInterfaces();

    livector1D listInterfaces;
    livector1D listBoundInterfaces;
	for (bitpit::Interface & interface : getGeometry()->getInterfaces()){
		if (!interface.isBorder()){
			std::array<double,3> centroidsVector = getGeometry()->evalCellCentroid(interface.getNeigh()) - getGeometry()->evalCellCentroid(interface.getOwner());
			std::array<double,3> normal = static_cast<bitpit::VolUnstructured*>(getGeometry()->getPatch())->evalInterfaceNormal(interface.getId());
			double angle = std::acos(dotProduct(centroidsVector,normal) / (norm2(centroidsVector)*norm2(normal)));
			m_maxSkewness = std::max(m_maxSkewness, angle);
			if (angle > m_maxSkewnessTol){
                listInterfaces.push_back(interface.getId());
			}
		}else{
			std::array<double,3> centroidsVector = getGeometry()->evalInterfaceCentroid(interface.getId()) -getGeometry()->evalCellCentroid(interface.getOwner());
			std::array<double,3> normal = static_cast<bitpit::VolUnstructured*>(getGeometry()->getPatch())->evalInterfaceNormal(interface.getId());
			double angle = std::acos(dotProduct(centroidsVector,normal) / (norm2(centroidsVector)*norm2(normal)));
			m_maxSkewnessBoundary = std::max(m_maxSkewnessBoundary, angle);
			if (angle > m_maxSkewnessBoundaryTol){
                listBoundInterfaces.push_back(interface.getId());
			}
		}
	}

    bool passed = ( listInterfaces.empty() && listBoundInterfaces.empty() );

    if (!passed){

        bitpit::PiercedVector<bitpit::Interface> & interfaces = getGeometry()->getInterfaces();
        bitpit::PiercedVector<bitpit::Cell> & cells = getGeometry()->getCells();
        bitpit::PiercedVector<bitpit::Vertex> & vertices = getGeometry()->getVertices();

        bitpit::PiercedVector<bitpit::Vertex> & skewVerts = m_skewness->getVertices();
        bitpit::PiercedVector<bitpit::Cell> & skewCells = m_skewness->getCells();

        for( long idI : listInterfaces){
                long idOwner = interfaces.at(idI).getOwner();
                long idNeigh = interfaces.at(idI).getNeigh();
                bitpit::Cell & cellOwner = cells.at(idOwner);
                bitpit::Cell & cellNeigh = cells.at(idNeigh);

                //push the cell if it is not ghost and if it's not already inside m_skewness
                if(cellOwner.isInterior() && !skewCells.exists(idOwner)){
                    bitpit::ConstProxyVector<long> vIds = cellOwner.getVertexIds();
                    for(long idV: vIds){
                        if(!skewVerts.exists(idV)){
                            m_skewness->addVertex(vertices.at(idV), idV);
                        }
                    }
                    m_skewness->addCell(cellOwner, idOwner);
                }

                //push the cell if it is not ghost and if it's not already inside m_skewness
                if(cellNeigh.isInterior() && !skewCells.exists(idNeigh)){
                    bitpit::ConstProxyVector<long> vIds = cellNeigh.getVertexIds();
                    for(long idV: vIds){
                        if(!skewVerts.exists(idV)){
                            m_skewness->addVertex(vertices.at(idV), idV);
                        }
                    }
                    m_skewness->addCell(cellNeigh, idNeigh);
                }
        }
	} // end of if passed.

#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &passed, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &m_maxSkewness, 1, MPI_DOUBLE, MPI_MAX, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &m_maxSkewnessBoundary, 1, MPI_DOUBLE, MPI_MAX, m_communicator);
#endif

	return passed;
}

/*! It check the face validity of the cells of the mesh
 * \return false if face validity tolerance is not satisfied
 */
bool
MeshChecker::checkFaceValidity()
{

	if (getGeometry()->getType() != 2){
		(*m_log)<<m_name<<" : face validity check active only for volume mesh -> check disabled" << std::endl;
		return true;
	}


    getGeometry()->updateInterfaces();

	//Prepare two structure areaGood and sumArea
	std::unordered_map<long, double> sumArea;
	std::unordered_map<long, double> areaGood;
	for (bitpit::Cell & cell : getGeometry()->getCells()){
		sumArea[cell.getId()] =  0.0;
		areaGood[cell.getId()] =  0.0;
	}
	for (bitpit::Interface & interface : getGeometry()->getInterfaces()){
		int n=1;
		if (!interface.isBorder())
			n=2;
		for (int i=0; i<n; i++){
			long idcell = interface.getOwnerNeigh()[i];
			long idinter = interface.getId();
			std::array<double,3> centroid = getGeometry()->getPatch()->evalCellCentroid(idcell);
			std::array<double,3> centerinterface = getGeometry()->getPatch()->evalInterfaceCentroid(idinter);
			std::array<double,3> normal = static_cast<bitpit::VolUnstructured*>(getGeometry()->getPatch())->evalInterfaceNormal(idinter);
			if (i==1){
				normal = -1.*normal;
			}
			double product = dotProduct((centerinterface-centroid),normal);
			double area = static_cast<bitpit::VolUnstructured*>(getGeometry()->getPatch())->evalInterfaceArea(idinter);
			if (product > 0.0){
				areaGood[idcell] = areaGood[idcell] + area;
			}
			sumArea[idcell] = sumArea[idcell] + area;
		}
	}

    // save the sick elements in a list.
    livector1D listSickCell;

	for (auto it = getGeometry()->getPatch()->internalBegin(); it!=getGeometry()->getPatch()->internalEnd(); ++it){
        long idcell = it.getId();
		double area = areaGood[idcell] / sumArea[idcell];
		m_minFaceValidity = std::min(m_minFaceValidity, area);
		if (area < m_minFaceValidityTol){
			listSickCell.push_back(idcell);
		}
	}
    //create the boolean
    bool passed = listSickCell.empty();

    // fill the structure.
    if (!passed){
        bitpit::PiercedVector<bitpit::Vertex> & vertices = getGeometry()->getVertices();
        bitpit::PiercedVector<bitpit::Cell> & cells = getGeometry()->getCells();

        bitpit::PiercedVector<bitpit::Vertex> & fvalVerts = m_facevalidity->getVertices();

        for(long idCell : listSickCell){
            bitpit::Cell & cell = cells.at(idCell);
            bitpit::ConstProxyVector<long> connIds = cell.getVertexIds();

            for(long idV : connIds){
                if(!fvalVerts.exists(idV)){
                    m_facevalidity->addVertex(vertices.at(idV), idV);
                }
            }
            m_facevalidity->addCell(cell, idCell);
        }
	}

#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &passed, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &m_minFaceValidity, 1, MPI_DOUBLE, MPI_MIN, m_communicator);
#endif

	return passed;
}


/*!
 * Initialize cells volumes
 */
void
MeshChecker::initializeVolumes()
{
	m_volumes.initialize(getGeometry(), MPVLocation::CELL, 0.);
	for (bitpit::Cell & cell : getGeometry()->getCells()){
		long id = cell.getId();
		m_volumes[id] = getGeometry()->evalCellVolume(id);
	}
}

/*! Clear auxiliary variables
 *
 */
void
MeshChecker::clear(){
	m_volumes.clear();
}

/*!
 * Plot mesh with failed check fields
 */
void
MeshChecker::plotOptionalResults(){
	if (getGeometry() == nullptr) return;
    if(m_qualityStatus == CMeshOutput::NOTRUN )  return;

    bool volumeEmpty = m_volume->isEmpty();
    bool volchangeEmpty = m_volumechange->isEmpty();
    bool skewEmpty = m_skewness->isEmpty();
    bool fvalEmpty = m_facevalidity->isEmpty();

#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &volumeEmpty, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &volchangeEmpty, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &skewEmpty, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &fvalEmpty, 1, MPI_C_BOOL, MPI_LAND, m_communicator);
#endif

	if (!volumeEmpty){
//#if MIMMO_ENABLE_MPI
//        m_volume->setPartitioned();
//#endif
        m_volume->getPatch()->getVTK().setDirectory(m_outputPlot+"/");
        m_volume->getPatch()->getVTK().setName(m_name+std::to_string(getId())+".volume");
        m_volume->getPatch()->write();
    }
	if (!volchangeEmpty){
//#if MIMMO_ENABLE_MPI
//		m_volumechange->setPartitioned();
//#endif
        m_volumechange->getPatch()->getVTK().setDirectory(m_outputPlot+"/");
        m_volumechange->getPatch()->getVTK().setName(m_name+std::to_string(getId())+".volume.ratio");
		m_volumechange->getPatch()->write();
    }
	if (!skewEmpty){
//#if MIMMO_ENABLE_MPI
//		m_skewness->setPartitioned();
//#endif
        m_skewness->getPatch()->getVTK().setDirectory(m_outputPlot+"/");
        m_skewness->getPatch()->getVTK().setName(m_name+std::to_string(getId())+".skewness");
		m_skewness->getPatch()->write();
    }
	if (!fvalEmpty){
//#if MIMMO_ENABLE_MPI
//		m_facevalidity->setPartitioned();
//#endif
        m_facevalidity->getPatch()->getVTK().setDirectory(m_outputPlot+"/");
        m_facevalidity->getPatch()->getVTK().setName(m_name+std::to_string(getId())+".face.validity");
		m_facevalidity->getPatch()->write();
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
MeshChecker::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	BaseManipulation::absorbSectionXML(slotXML, name);

	//TODO XML SECTIONS

   if(slotXML.hasOption("MinVolTOL")){
       std::string input = slotXML.get("MinVolTOL");
       input = bitpit::utils::string::trim(input);
       double val = 1.0E-14;
       if(!input.empty()){
           std::stringstream ss(input);
           ss>>val;
       }
       setMinimumVolumeTolerance(val);
   }
   if(slotXML.hasOption("MaxVolTOL")){
       std::string input = slotXML.get("MaxVolTOL");
       input = bitpit::utils::string::trim(input);
       double val = 1.0E10;
       if(!input.empty()){
           std::stringstream ss(input);
           ss>>val;
       }
       setMaximumVolumeTolerance(val);
   }

   if(slotXML.hasOption("MaxSkewTOL")){
       std::string input = slotXML.get("MaxSkewTOL");
       input = bitpit::utils::string::trim(input);
       double val = 85.0*BITPIT_PI/180;
       if(!input.empty()){
           std::stringstream ss(input);
           ss>>val;
       }
       setMaximumSkewnessTolerance(val);
   }

   if(slotXML.hasOption("MaxBoundarySkewTOL")){
       std::string input = slotXML.get("MaxBoundarySkewTOL");
       input = bitpit::utils::string::trim(input);
       double val = 85.0*BITPIT_PI/180;
       if(!input.empty()){
           std::stringstream ss(input);
           ss>>val;
       }
       setMaximumBoundarySkewnessTolerance(val);
   }

   if(slotXML.hasOption("MinFaceValidityTOL")){
       std::string input = slotXML.get("MinFaceValidityTOL");
       input = bitpit::utils::string::trim(input);
       double val = 0.5;
       if(!input.empty()){
           std::stringstream ss(input);
           ss>>val;
       }
       setMinimumFaceValidityTolerance(val);
   }
   if(slotXML.hasOption("MinVolChangeTOL")){
       std::string input = slotXML.get("MinVolChangeTOL");
       input = bitpit::utils::string::trim(input);
       double val = 1.0E-5;
       if(!input.empty()){
           std::stringstream ss(input);
           ss>>val;
       }
       setMinimumVolumeChangeTolerance(val);
   }
   if(slotXML.hasOption("ResumeFile")){
       std::string input = slotXML.get("ResumeFile");
       input = bitpit::utils::string::trim(input);
       bool val = false;
       if(!input.empty()){
           std::stringstream ss(input);
           ss>>val;
       }
       setPrintResumeFile(val);
   }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
MeshChecker::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	BaseManipulation::flushSectionXML(slotXML, name);

    {
        std::stringstream ss;
        ss << std::scientific<<std::setprecision(12)<<m_minVolumeTol;
        slotXML.set("MinVolTOL", ss.str());
    }

    {
        std::stringstream ss;
        ss << std::scientific<<std::setprecision(12)<<m_maxVolumeTol;
        slotXML.set("MaxVolTOL", ss.str());
    }
    {
        std::stringstream ss;
        ss << std::scientific<<std::setprecision(12)<<m_maxSkewnessTol;
        slotXML.set("MaxSkewTOL", ss.str());
    }

    {
        std::stringstream ss;
        ss << std::scientific<<std::setprecision(12)<<m_maxSkewnessBoundaryTol;
        slotXML.set("MaxBoundarySkewTOL", ss.str());
    }
    {
        std::stringstream ss;
        ss << std::scientific<<std::setprecision(12)<<m_minFaceValidityTol;
        slotXML.set("MinFaceValidityTOL", ss.str());
    }

    {
        std::stringstream ss;
        ss << std::scientific<<std::setprecision(12)<<m_minVolumeChangeTol;
        slotXML.set("MinVolChangeTOL", ss.str());
    }

    slotXML.set("ResumeFile", std::to_string(int(m_printResumeFile)));
};

/*!
    Print resume file, of performing mesh check
*/
void MeshChecker::printResumeFile(){

    long countVolumeCells = m_volume->getNCells();
    long countVolumeChangeCells = m_volumechange->getNCells();
    long countSkewCells = m_skewness->getNCells();
    long countFValCells = m_facevalidity->getNCells();

#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &countVolumeCells,1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &countVolumeChangeCells,1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &countSkewCells,1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, &countFValCells,1, MPI_LONG, MPI_SUM, m_communicator);

    if(m_rank == 0)
#endif
    {
        std::string dir = m_outputPlot + "/";
        std::string name = m_name +std::to_string(getId())+"_ResumeFile.dat";
        std::string file = dir + name;
        std::ofstream out;
        out.open(file);
        if(out.is_open()){

            out<<"#mimmo MeshChecker Resume File"<<std::endl;
            out<<std::endl;
            out<<std::endl;
            out<<"#Sick Cells Count: "<<std::endl;
            out<<std::endl;
            out<<"VOLUME       : "<<countVolumeCells<<std::endl;
            out<<"VOLUMECHANGE : "<<countVolumeChangeCells<<std::endl;
            out<<"SKEWNESS     : "<<countSkewCells<<std::endl;
            out<<"FACEVALIDITY : "<<countFValCells<<std::endl;

            out<<std::endl;
            out<<std::endl;
            out<<"#Found indicators values vs Thresholds: "<<std::endl;
            out<<std::endl;
            out<<std::scientific<<"MinVol Threshold          : "<< m_minVolumeTol           <<" , MinVol found          : "<<m_minVolume<<std::endl;
            out<<std::scientific<<"MaxVol Threshold          : "<< m_maxVolumeTol           <<" , MaxVol found          : "<<m_maxVolume<<std::endl;
            out<<std::scientific<<"MinVolChange Threshold    : "<< m_minVolumeChangeTol     <<" , MinVolChange found    : "<<m_minVolumeChange<<std::endl;
            out<<std::scientific<<"MaxSkewness Threshold     : "<< m_maxSkewnessTol         <<" , MaxSkewness found     : "<<m_maxSkewness<<std::endl;
            out<<std::scientific<<"MaxBndSkewness Threshold  : "<< m_maxSkewnessBoundaryTol <<" , MaxBndSkewness found  : "<<m_maxSkewnessBoundary <<std::endl;
            out<<std::scientific<<"MinFaceVal Threshold      : "<< m_minFaceValidityTol     <<" , MinFaceVal found      : "<<m_minFaceValidity <<std::endl;

            out.close();
        }else{
            throw std::runtime_error("MeshChecker cannot open ResumeFile to write.");
        }

    }
#if MIMMO_ENABLE_MPI
    MPI_Barrier(m_communicator);
#endif

}



}
