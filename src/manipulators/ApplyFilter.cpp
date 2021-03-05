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
\*---------------------------------------------------------------------------*/
#include "ApplyFilter.hpp"

namespace mimmo{

/*!
 * Default constructor of ApplyFilter
 */
ApplyFilter::ApplyFilter():Apply(){
    m_name = "mimmo.ApplyFilter";
    m_factor = std::numeric_limits<double>::max();
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ApplyFilter::ApplyFilter(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.ApplyFilter";
	m_factor = 1.;

	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);
	if(input == "mimmo.ApplyFilter"){
		absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!Default destructor of ApplyFilter
 */
ApplyFilter::~ApplyFilter(){};

/*!
    Copy constructor of ApplyFilter.
 */
ApplyFilter::ApplyFilter(const ApplyFilter & other):Apply(other){
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void ApplyFilter::swap(ApplyFilter & x) noexcept
{
	Apply::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
ApplyFilter::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dmpvecarr3E*, ApplyFilter>(this, &ApplyFilter::setInput, M_GDISPLS, true, 1));
    built = (built && createPortIn<dmpvector1D*, ApplyFilter>(this, &ApplyFilter::setScalarInput, M_SCALARFIELD, true, 1));
    built = (built && createPortIn<dmpvector1D*, ApplyFilter>(this, &ApplyFilter::setFilter, M_FILTER));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, Apply>(this, &BaseManipulation::setGeometry, M_GEOM));

	built = (built && createPortOut<dmpvecarr3E*, ApplyFilter>(this, &ApplyFilter::getOutput, M_GDISPLS));

	m_arePortsBuilt = built;
};


/*!Execution command.
 * It applies the filter field to the deformation field stored in the input of base class (casting the input
 * for apply object to dvecarr3E). Scaling factor is applied too if passed by the user.
 * After exec() the linked geometry is not modified.
 */
void
ApplyFilter::execute(){

	checkInput();

	m_output = m_input;

	if (m_factor != std::numeric_limits<double>::max()){
	    for (auto & val : m_output){
	        val *= m_factor;
	    }
	}

};


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ApplyFilter::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Scaling")){
        std::string input = slotXML.get("Scaling");
        input = bitpit::utils::string::trim(input);
        double val = 1.;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>val;
        }
        setScaling(val);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
ApplyFilter::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("Scaling", std::to_string(m_factor));

};

/*!
 * Check if the input is related to the target geometry.
 * Check the coherence of filter field.
 * If not create a zero input field.
 */
void
ApplyFilter::checkInput(){

	bool check = false;
	MimmoSharedPointer<MimmoObject> linked_geometry = nullptr;
	if (m_geometry) linked_geometry = m_geometry;
	if (m_input.size()){
	    if (!linked_geometry){
	        linked_geometry = m_input.getGeometry();
	    }
	    else{
	        m_input.setGeometry(linked_geometry);
	    }
		check = m_input.getDataLocation() == mimmo::MPVLocation::POINT;
		check = check && m_input.completeMissingData({{0.0,0.0,0.0}});
	}
	if (m_scalarinput.size()){
        if (!linked_geometry){
            linked_geometry = m_scalarinput.getGeometry();
        }
        else{
            m_scalarinput.setGeometry(linked_geometry);
        }
		check = m_scalarinput.getDataLocation() == mimmo::MPVLocation::POINT;
		check = check && m_scalarinput.completeMissingData(0.0);
		if (check){
			bitpit::SurfaceKernel* skernel = static_cast<bitpit::SurfaceKernel*>(m_scalarinput.getGeometry()->getPatch());
			bitpit::PiercedVector<darray3E> vNormals;
			bitpit::ConstProxyVector<long> verts;
			std::size_t size;
			long idN;
			for(const auto & cell: m_scalarinput.getGeometry()->getCells()){
				verts = cell.getVertexIds();
				size = verts.size();
				for(std::size_t i=0; i<size; ++i){
					idN = verts[i];
					if(!vNormals.exists(idN)){
						vNormals.insert(idN, skernel->evalVertexNormal(cell.getId(), i));
					}
				}
			}
			m_input.clear();
			for(const auto & vertex: m_scalarinput.getGeometry()->getVertices()){
				m_input.setDataLocation(mimmo::MPVLocation::POINT);
				m_input.setGeometry(m_scalarinput.getGeometry());
				long id = vertex.getId();
				m_input.insert(id, darray3E(m_scalarinput[id]*vNormals[id]));
			}
		}
	}
	if (m_filter.size()){
	    check = m_filter.getDataLocation() == mimmo::MPVLocation::POINT;
        // Force filter to link the geometry linked by deformation field
	    bool geometry_check = m_filter.getGeometry() == linked_geometry;
	    if (!geometry_check){
	        m_filter.setGeometry(linked_geometry);
	    }
	    // Complete missing data force filter to 0.
	    check = check && m_filter.completeMissingData(0.0);
	    if (check){
	        // ApplyFilter filter to input deformation field
	        for (const auto & vertex : linked_geometry->getVertices()){
	            long ID = vertex.getId();
	            m_input[ID] *= m_filter[ID];
	        }
	    }
	}
	if (!check){
		m_log->setPriority(bitpit::log::Verbosity::NORMAL);
		(*m_log) << m_name + " : not valid inputs found" << std::endl;
		throw std::runtime_error (m_name + " : not valid inputs found");
	}
};

}
