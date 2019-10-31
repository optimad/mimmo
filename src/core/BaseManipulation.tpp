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
namespace mimmo{

/*!It adds an output port to the object.The input ports with compatible tags can receive data
 * form this output port. The compatibility between ports with same tag is automatic.
 *
 * \param[in] var_ Pointer to T member to communicate.
 * \param[in] portS marker of the port; the port will be created at portS-th slot of the output ports of the object.
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortOut(T* var_, PortID portS){
	bool check = false;

	//checking if portS containerTAG and dataTAG are registered in mimmo::PortManager

	if(!mimmo::PortManager::instance().containsPort(portS) || m_portOut.count(portS) != 0 ){
		(*m_log)<<"Unable to physically create the Output Port."<<std::endl;
		(*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
		return check;
	}

	auto info = mimmo::PortManager::instance().getPortData(portS);
	DataType datat(info.container, info.datatype);
	PortOutT<T, O>* portOut = new PortOutT<T, O>(var_, datat);
	m_portOut[portS] = portOut;
	check = true;
	return(check);
}

/*!It adds an output port to the object.The input ports with compatible tags can receive data
 * form this output port. The compatibility between ports with same tag is automatic.
 *
 * \param[in] obj_ Pointer to BaseManipulation object.
 * \param[in] getVar_ Get function of the object (copy return) linked by the sender port.
 * \param[in] portS marker of the port; the port will be created at portS-th slot of the output ports of the object.
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortOut(O* obj_, T (O::*getVar_)(), PortID portS){
	bool check = false;

	//checking if portS containerTAG and dataTAG are registered in mimmo::PortManager

	if(!mimmo::PortManager::instance().containsPort(portS) || m_portOut.count(portS) != 0 ){
		(*m_log)<<"Unable to physically create the Output Port."<<std::endl;
		(*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
		return check;
	}

	auto info = mimmo::PortManager::instance().getPortData(portS);
	DataType datat(info.container, info.datatype);
	PortOutT<T, O>* portOut = new PortOutT<T, O>(obj_, getVar_, datat);
	m_portOut[portS] = portOut;
	check = true;
	return(check);
}

/*!It adds an input port to the object. The output ports with compatible container/data tags can send data
 * to this input port. The compatibility between ports with same tag is automatic.
 *
 * \param[in] var_ Pointer to member to fill.
 * \param[in] portR marker of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] mandatory does the port have to be mandatorily linked?
 * \param[in] family tag of family of mandatory ports (family=0 independent ports with no alternative).
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(T* var_, PortID portR, bool mandatory, int family ){
	bool check = false;

	//checking if portR containerTAG and dataTAG are registered in mimmo::PortManager

	if(!mimmo::PortManager::instance().containsPort(portR) || m_portOut.count(portR) != 0 ){
		(*m_log)<<"Unable to physically create the Input Port."<<std::endl;
		(*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
		return check;
	}

	auto info = mimmo::PortManager::instance().getPortData(portR);
	DataType datat(info.container, info.datatype);
	PortInT<T, O>* portIn = new PortInT<T, O>(var_, datat, mandatory, family);
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}

/*!It adds an input port to the object.The output ports compatible tags can send data
 * to this input port. The compatibility between ports with same tag is automatic.
 *
 * \param[in] obj_ Pointer to BaseManipulation object.
 * \param[in] setVar_ Get function of the object (copy return) linked by the sender port.
 * \param[in] portR marker of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] mandatory does the port have to be mandatorily linked?
 * \param[in] family tag of family of mandatory ports (family=0 independent ports with no alternative).
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(O* obj_, void (O::*setVar_)(T), PortID portR, bool mandatory, int family ){
	bool check = false;

	//checking if portR containerTAG and dataTAG are registered in mimmo::PortManager
	if(!mimmo::PortManager::instance().containsPort(portR) || m_portOut.count(portR) != 0 ){
		(*m_log)<<"Unable to physically create the Input Port."<<std::endl;
		(*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
		return check;
	}

	auto info = mimmo::PortManager::instance().getPortData(portR);
	DataType datat(info.container, info.datatype);
	PortInT<T, O>* portIn = new PortInT<T, O>(obj_, setVar_, datat, mandatory, family);
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}

/*!
 * Write an input geometry given as MimmoObject pointer with an input data field. Input geometry and the geometry linked in the
 * MimmoPiercedVector of data have to be consistent, otherwise the data is skipped.
 * It writes into the m_output directory with m_name+m_counter. It adds the input data field by deducing the type.
 * Allowed types are : all scalars that verify is_floating_point or is_integral functions; mimmo vector types, i.e. std::array<double,3>.
 * \param[in] geometry MimmoObject pointer to target geometry
 * \param[in] data MimmoPiercedVector with data field to write
 */
template<typename mpv_t>
void
BaseManipulation::write(MimmoObject* geometry, MimmoPiercedVector<mpv_t> & data)
{

	//Check data and geometry linked
	if (data.size() == 0 || data.getGeometry() == nullptr){
		(*m_log) << " Warning: data to write not consistent with geometry in " << m_name << "; skip data" << std::endl;
		write(geometry);
		return;
	}
	if (geometry == nullptr){
		(*m_log) << " Warning: geometry null during writing " << m_name << std::endl;
		return;
	}

	//Recover location of data
	bitpit::VTKLocation loc = bitpit::VTKLocation::UNDEFINED;
	switch(data.getDataLocation()){
	case MPVLocation::POINT :
		loc = bitpit::VTKLocation::POINT;
		break;
	case MPVLocation::CELL :
		loc = bitpit::VTKLocation::CELL;
		break;
	default:
		(*m_log)<<" Warning: Undefined Reference Location in plotOptionalResults of "<<m_name<<std::endl;
		(*m_log)<<" Interface or Undefined locations are not supported in VTU writing." <<std::endl;
		return;
		break;
	}

	//check size of field and adjust missing values to zero for writing purposes only.
	if(!data.completeMissingData(mpv_t())){
		(*m_log) << " Warning: error during complete missing data to write in " << m_name << "; skip data" << std::endl;
		write(geometry);
		return;
	}

	//Deduce data type and add data to vtk object
	bitpit::VTKFieldType fieldtype;
	if (std::is_integral<mpv_t>::value == true || std::is_floating_point<mpv_t>::value == true){
		fieldtype = bitpit::VTKFieldType::SCALAR;
	}
	else if (std::is_same<mpv_t, std::array<double,3>>::value == true){
		fieldtype = bitpit::VTKFieldType::VECTOR;
	}
	else{
		(*m_log) << " Warning: data type to write not allowed in " << m_name << "; skip data" << std::endl;
		write(geometry);
		return;
	}

	std::vector<mpv_t> field = data.getDataAsVector();

	//Check if geometry is a point cloud or not
//	if(geometry->getType() != 3){

		geometry->getPatch()->getVTK().addData(data.getName(), fieldtype, loc, field);

//	}else{
//
//		if(loc == bitpit::VTKLocation::CELL){
//			(*m_log) << " Warning: attempt writing Cell data field on cloud of points in " << m_name << "; skip data" << std::endl;
//			write(geometry);
//			return;
//		}
		geometry->getPatch()->getVTK().addData(data.getName(), bitpit::VTKFieldType::SCALAR, loc, field);

//	}

	write(geometry);

}

/*!
 * Write an input geometry given as MimmoObject pointer with a series of input data field (template variadic function).
 * Input geometry and the geometry linked in the MimmoPiercedVectors of data have to be consistent, otherwise the data
 * that doesn't satisfy the coherence is skipped.
 * It writes into the m_output directory with m_name+m_counter. It adds the input data fields by deducing the type.
 * Allowed types are : all scalars that verify is_floating_point or is_integral functions; mimmo vector types, i.e. std::array<double,3>.
 * \param[in] geometry MimmoObject pointer to target geometry
 * \param[in] data MimmoPiercedVector with data fields to write
 * \param[in] ... MimmoPiercedVector series with data fields to write (template variadic function)
 */
template<typename mpv_t, typename... Args>
void
BaseManipulation::write(MimmoObject* geometry, MimmoPiercedVector<mpv_t> & data, Args ... args)
{
	//Check data and geometry linked
	if (data.size() == 0 || data.getGeometry() == nullptr){
		(*m_log) << " Warning: data to write not consistent with geometry in " << m_name << "; skip data" << std::endl;
		write(geometry, args...);
		return;
	}
	if (geometry == nullptr){
		(*m_log) << " Warning: geometry null during writing " << m_name << std::endl;
		return;
	}

	bitpit::VTKLocation loc = bitpit::VTKLocation::UNDEFINED;
	switch(data.getDataLocation()){
	case MPVLocation::POINT :
		loc = bitpit::VTKLocation::POINT;
		break;
	case MPVLocation::CELL :
		loc = bitpit::VTKLocation::CELL;
		break;
	default:
		(*m_log)<<" Warning: Undefined Reference Location in plotOptionalResults of "<<m_name<<std::endl;
		(*m_log)<<" Interface or Undefined locations are not supported in VTU writing." <<std::endl;
		return;
		break;
	}

	//check size of field and adjust missing values to zero for writing purposes only.
	if(!data.completeMissingData(mpv_t())){
		(*m_log) << " Warning: error during complete missing data to write in " << m_name << "; skip data" << std::endl;
		write(geometry, args...);
		return;
	}

	//Deduce data type and add data to vtk object
	bitpit::VTKFieldType fieldtype;
	if (std::is_integral<mpv_t>::value == true || std::is_floating_point<mpv_t>::value == true){
		fieldtype = bitpit::VTKFieldType::SCALAR;
	}
	else if (std::is_same<mpv_t, std::array<double,3>>::value == true){
		fieldtype = bitpit::VTKFieldType::VECTOR;
	}
	else{
		(*m_log) << " Warning: data type to write not allowed in " << m_name << "; skip data" << std::endl;
		write(geometry, args...);
		return;
	}

	std::vector<mpv_t> field = data.getDataAsVector();

	//Check if geometry is a point cloud or not
//	if(geometry->getType() != 3){

		geometry->getPatch()->getVTK().addData(data.getName(), bitpit::VTKFieldType::SCALAR, loc, field);

//	}else{
//
//		if(loc == bitpit::VTKLocation::CELL){
//			(*m_log) << " Warning: attempt writing Cell data field on cloud of points in " << m_name << "; skip data" << std::endl;
//			write(geometry, args...);
//			return;
//		}
		geometry->getPatch()->getVTK().addData(data.getName(), bitpit::VTKFieldType::SCALAR, loc, field);

//	}

	write(geometry, args...);

	geometry->getPatch()->getVTK().removeData(data.getName());

}


};
