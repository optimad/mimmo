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

#include "SelectField.hpp"
#include "SkdTreeUtils.hpp"
#include "ExtractFields.hpp"
#include <unordered_map>

namespace mimmo{

/*!
 * Default constructor.
 * \param[in] loc MPVLocation of fields data: POINT, CELL or INTERFACE.
 */
SelectVectorField::SelectVectorField(MPVLocation loc):SelectField(){
    m_name = "mimmo.SelectVectorField";
    m_loc = loc;
    if(m_loc == MPVLocation::UNDEFINED) m_loc= MPVLocation::POINT;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectVectorField::SelectVectorField(const bitpit::Config::Section & rootXML){

    std::string fallback_name = "ClassNONE";
    std::string fallback_loc  = "-1";

    std::string input_name = rootXML.get("ClassName", fallback_name);
    input_name = bitpit::utils::string::trim(input_name);

    std::string input_loc = rootXML.get("DataLocation", fallback_loc);
    input_loc = bitpit::utils::string::trim(input_loc);

    int loc = std::stoi(input_loc);
    if(loc > 0 && loc < 4){
        m_loc  =static_cast<MPVLocation>(loc);
    }else{
        m_loc = MPVLocation::POINT;
    }
    m_name = "mimmo.SelectVectorField";

    if(input_name == "mimmo.SelectVectorField" ){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor
 */
SelectVectorField::~SelectVectorField(){}

/*!
 * Copy constructor
 */
SelectVectorField::SelectVectorField(const SelectVectorField & other):SelectField(other){
    m_loc = other.m_loc;
    m_fields = other.m_fields;
    m_result = other.m_result;
}

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void SelectVectorField::swap(SelectVectorField & x) noexcept
{
    std::swap(m_loc, x.m_loc);
    std::swap(m_fields, x.m_fields);
    //std::swap(m_result, x.m_result);
    m_result.swap(x.m_result);
    SelectField::swap(x);
}

/*!
 * Build the ports of the class;
 */
void
SelectVectorField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<std::vector<dmpvecarr3E*>, SelectVectorField>(this, &mimmo::SelectVectorField::setFields, M_VECVECTORFIELDS, true, 1));
    built = (built && createPortIn<dmpvecarr3E*, SelectVectorField>(this, &mimmo::SelectVectorField::addField, M_VECTORFIELD, true, 1));
    built = (built && createPortOut<dmpvecarr3E*, SelectVectorField>(this, &mimmo::SelectVectorField::getSelectedField, M_VECTORFIELD));

    SelectField::buildPorts();
    m_arePortsBuilt = built;
}

/*!
 * Get Selected field.
 * \return Selected field
 */
dmpvecarr3E*
SelectVectorField::getSelectedField(){
    return &m_result;
}

/*!
 * Set list of input fields.
 * \param[in] fields scalar fields
 */
void
SelectVectorField::setFields(std::vector<dmpvecarr3E*> fields){
    m_fields.clear();
    m_fields.reserve(fields.size());
    for(dmpvecarr3E * ff : fields){
        if(!ff) continue;
        if(ff->getDataLocation() == m_loc && ff->getGeometry()!= nullptr){
            m_fields.push_back(*ff);
        }
    }
    m_fields.shrink_to_fit();
}

/*!
 * Add a field to the list of input fields.
 * \param[in] field scalar field
 */
void
SelectVectorField::addField(dmpvecarr3E *field){
    if(!field) return;
    if(field->getDataLocation() == m_loc && field->getGeometry() != nullptr){
        m_fields.push_back(*field);
    }
}

/*!
 * Clear content of the class
 */
void
SelectVectorField::clear(){
    m_fields.clear();
    m_result.clear();
    SelectField::clear();
}

/*!
 * Plot Selected field on its geometry.
 */
void
SelectVectorField::plotOptionalResults(){

	m_result.setName(m_fieldname);
	write(m_result.getGeometry(), m_result);

}

/*!
 * Select your original field along the input geometry provided
 * \return true if Select without errors
 */
bool
SelectVectorField::mSelect(){

	m_result.clear();

	switch(m_mode) {

	case SelectType::GEOMETRY:
	{
		//Check if geometry is linked
		if (getGeometry() == nullptr) return false;

		//Extract by link to geometry
		for (const auto & field : m_fields){
			if (field.getGeometry() == getGeometry()){
				m_result = field;
				//geometry, location and name are copied from field.
				return true;
			}
		}
		break;
	}

	case SelectType::NAME:
	{
		//Extract by link to geometry
		for (const auto & field : m_fields){
			if (field.getName() == getFieldName()){
				m_result = field;
				//geometry, location and name are copied from field.
				return true;
			}
		}
		break;
	}

	case SelectType::MAPPING:
	{
		//Check if geometry is linked
		if (getGeometry() == nullptr) return false;

		//Extract by geometric mapping if active and if no positive match is found.
		if (getGeometry()->getType() != 3){

			m_result.setGeometry(getGeometry());
			m_result.setDataLocation(m_loc);

			ExtractVectorField * ef = new ExtractVectorField();
			ef->setGeometry(getGeometry());
			ef->setMode(ExtractMode::MAPPING);
			ef->setTolerance(m_tol);

			//create map for overlapping ids purpose;
			std::unordered_map<long, int> idRepetition;

			for (dmpvecarr3E & field : m_fields){
				ef->setField(&field);
				bool check = ef->extract();
				if(!check) continue;

				dmpvecarr3E * temp = ef->getExtractedField();
				dmpvecarr3E::iterator itB;
				auto itE = temp->end();
				long id;
				for(itB = temp->begin(); itB != itE; ++itB){
					id = itB.getId();
					if(!m_result.exists(id)){
						m_result.insert(id, *itB);
					}else{
						m_result[id] += *itB;
						if(idRepetition.count(id)>0){
							++idRepetition[id];
						}else{
							idRepetition[id] = 2;
						}
					}
				}
			}

			//resolve overlapping ids by averaging correspondent value;

			for(auto &itval : idRepetition){
				m_result[itval.first] /= double(itval.second);
			}

			delete ef;
            if (m_result.size() > 0) return true;
        }
        break;
    }

    } // end switch mode

    //if you are here something was wrong.
    return false;
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectVectorField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    if(slotXML.hasOption("DataLocation")){
        std::string input = slotXML.get("DataLocation");
        input = bitpit::utils::string::trim(input);
        int temp = -1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        if(int(m_loc) != temp){
            (*m_log)<<"Error absorbing DataLocation in "<<m_name<<". Class and read locations mismatch"<<std::endl;
            throw std::runtime_error (m_name + " : xml absorbing failed");
        }
    }

    SelectField::absorbSectionXML(slotXML, name);


};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectVectorField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    slotXML.set("DataLocation", std::to_string(int(m_loc)));
    SelectField::flushSectionXML(slotXML, name);
};



}
