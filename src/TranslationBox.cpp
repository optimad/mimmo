/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#include "TranslationBox.hpp"

using namespace mimmo;

REGISTER_MANIPULATOR("MiMMO.TranslationBox", "translationbox");
/*!Default constructor of TranslationBox
 */
TranslationBox::TranslationBox(darray3E direction){
	m_direction = direction;
	m_name = "MiMMO.TranslationBox";
};

/*!Default destructor of TranslationBox
 */
TranslationBox::~TranslationBox(){};

/*!Copy constructor of TranslationBox.
 */
TranslationBox::TranslationBox(const TranslationBox & other):BaseManipulation(other){
	m_direction = other.m_direction;
};

/*!Assignement operator of TranslationBox.
 */
TranslationBox & TranslationBox::operator=(const TranslationBox & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_direction = other.m_direction;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void TranslationBox::buildPorts(){
	bool built = true;
	built = (built && createPortIn<darray3E, TranslationBox>(&m_origin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, TranslationBox>(&m_direction, PortType::M_AXIS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, TranslationBox>(&m_alpha, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<darray3E, TranslationBox>(this, &mimmo::TranslationBox::getOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<darray3E, TranslationBox>(this, &mimmo::TranslationBox::getDirection, PortType::M_AXIS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<double, TranslationBox>(this, &mimmo::TranslationBox::getTranslation, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!It gets the direction of the translation.
 * \return Direction of translation.
 */
darray3E
TranslationBox::getDirection(){
	return(m_direction);
}

/*!It gets the value of the translation.
 * \return Value of translation.
 */
double
TranslationBox::getTranslation(){
	return(m_alpha);
}

/*!It gets the original position of the point to be translated (before the execution)
 * or the position of the translated point (after the execution of the object).
 * \return Position of the point.
 */
darray3E
TranslationBox::getOrigin(){
	return(m_origin);
}

/*!It sets the direction of the translation.
 * \param[in] direction Direction of translation.
 */
void
TranslationBox::setDirection(darray3E direction){
	m_direction = direction;
}

/*!It sets the value of the translation.
 * \param[in] alpha Value of translation.
 */
void
TranslationBox::setTranslation(double alpha){
	m_alpha = alpha;
}

/*!It sets the original coordinates of the point to be translated.
 * \param[in] origin Position of the point.
 */
void
TranslationBox::setOrigin(darray3E origin){
	m_origin = origin;
}

/*!Execution command. It modifies the coordinates of the origin
 * with the translation conditions.
 * The result of the translation is stored in member result of base class
 * and in the member m_origin.
 * After exec() the original point coordinates will be permanently modified.
 */
void
TranslationBox::execute(){
	for (int i=0; i<3; i++){
			m_origin[i] += m_alpha * m_direction[i];
	}
	return;
};

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to perform translation of 3D point
 * 
 * --> Absorbing data:
 * 		Origin: point tha need to be translated
 * 		Direction: translation direction
 * 		Translation : entity of translation
 * 
 * \param[in] slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void TranslationBox::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	if(slotXML.hasOption("Origin")){
		std::string input = slotXML.get("Origin");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setOrigin(temp);
	} 
	
	if(slotXML.hasOption("Direction")){
		std::string input = slotXML.get("Direction");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setDirection(temp);
	} 
	
	if(slotXML.hasOption("Translation")){
		std::string input = slotXML.get("Translation");
		input = bitpit::utils::trim(input);
		double temp = 0.0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
		setTranslation(temp);
	} 
	
};	
/*!
 * Write settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to perform translation of a 3D point.
 * 
 * --> Flushing data// how to write it on XML:
 * 		ClassName : name of the class as "MiMMO.RotationBox"
 * 		ClassID	  : integer identifier of the class	
 * 		Origin: point tha need to be translated
 * 		Direction: translation direction
 * 		Translation : entity of translation
 * 
 * \param[in] slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void TranslationBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_origin[0]<<'\t'<<m_origin[1]<<'\t'<<m_origin[2];
		slotXML.set("Origin", ss.str());
	}
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_direction[0]<<'\t'<<m_direction[1]<<'\t'<<m_direction[2];
		slotXML.set("Direction", ss.str());
	}
	
	slotXML.set("Translation", std::to_string(m_alpha));	
	
};	




