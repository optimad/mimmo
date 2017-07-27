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

#include "ExtractFields.hpp"

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!
 * Default constructor.
 */
ExtractVectorField::ExtractVectorField():ExtractField(){
    m_name = "mimmo.ExtractVectorField";
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ExtractVectorField::ExtractVectorField(const bitpit::Config::Section & rootXML){

    std::string fallback_name = "ClassNONE";
    std::string input_name = rootXML.get("ClassName", fallback_name);
    input_name = bitpit::utils::string::trim(input_name);

    m_name = "mimmo.ExtractVectorField";

    if(input_name == "mimmo.ExtractVectorField"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor
 */
ExtractVectorField::~ExtractVectorField(){
    m_field.clear();
    m_result.clear();
}

/*!
 * Build the ports of the class;
 */
void
ExtractVectorField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dmpvecarr3E, ExtractVectorField>(this, &mimmo::ExtractVectorField::setField, PortType::M_VECTORFIELD, mimmo::pin::containerTAG::MPVECARR3, mimmo::pin::dataTAG::FLOAT, true, 1));
    built = (built && createPortOut<dmpvecarr3E, ExtractVectorField>(this, &mimmo::ExtractVectorField::getExtractedField, PortType::M_VECTORFIELD, mimmo::pin::containerTAG::MPVECARR3, mimmo::pin::dataTAG::FLOAT));

    ExtractField::buildPorts();
    m_arePortsBuilt = built;
}

/*!
 * Get extracted field.
 * \return extracted field
 */
dmpvecarr3E
ExtractVectorField::getExtractedField(){
    return m_result;
}

/*!
 * Set input field for extraction.
 * \param[in] field input field related to whole geometry
 */
void
ExtractVectorField::setField(dmpvecarr3E field){
    m_field = field;
}

/*!
 * Clear content of the class
 */
void
ExtractVectorField::clear(){
    m_field.clear();
    m_result.clear();
    ExtractField::clear();
}

/*!
 * Plot extracted field over the target geometry
 */
void 
ExtractVectorField::plotOptionalResults(){

    if (m_result.size() == 0 || getGeometry() == NULL) return;

    bitpit::VTKLocation loc = bitpit::VTKLocation::POINT;

    dvecarr3E field;
    for (const auto & v : getGeometry()->getVertices()){
        field.push_back(m_result[v.getId()]);
    }
    
    bitpit::VTKElementType cellType = getGeometry()->desumeElement();
    liimap mapDataInv;
    dvecarr3E points = getGeometry()->getVertexCoords(&mapDataInv);

    if (cellType == bitpit::VTKElementType::UNDEFINED) return;

    if(cellType != bitpit::VTKElementType::VERTEX){
        ivector2D connectivity = getGeometry()->getCompactConnectivity(mapDataInv);
        bitpit::VTKUnstructuredGrid output(".",m_name+std::to_string(getClassCounter()),cellType);
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
        output.setDimensions(connectivity.size(), points.size());
        output.addData("field", bitpit::VTKFieldType::VECTOR, loc, field);
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
        output.addData("field", bitpit::VTKFieldType::VECTOR, loc, field);
        output.setCodex(bitpit::VTKFormat::APPENDED);
        output.write();
    }

}


/*!
 * Extract field from your original field by using the provided target geometry
 * \return true if extract without errors
 */
bool
ExtractVectorField::extract(){

    if (getGeometry() == NULL || m_field.getGeometry() == NULL) return false;

    m_result.clear();

    switch(m_mode){
    case ExtractMode::ID :
    {
        //Extract by IDs
        for (const auto & ID : getGeometry()->getVertices().getIds()){
            if (m_field.exists(ID)){
                m_result.insert(ID, m_field[ID]);
            }
        }
    }
    break;

    case ExtractMode::PID :
    {
        //Extract by PIDs
        std::set<int> pids;
        for (const auto & cell : getGeometry()->getCells()){
            pids.insert(cell.getPID());
        }
        for (const auto & cell : m_field.getGeometry()->getCells()){
            if (pids.count(cell.getPID())){
                for (const auto & id : m_field.getGeometry()->getCellConnectivity(cell.getId())){
                    if (!m_result.exists(id))
                        m_result.insert(id, m_field[id]);
                }
            }
        }
    }
    break;

    case ExtractMode::MAPPING :
    {
        //Extract by geometric mapping
        double tol = 1.0e-08;
        livector1D result = mimmo::bvTreeUtils::selectByPatch(m_field.getGeometry()->getBvTree(), getGeometry()->getBvTree(), tol);
        for (const auto & idC : result){
            for (const auto & id : getGeometry()->getCellConnectivity(idC)){
                if (!m_result.exists(id))
                    m_result.insert(id, m_field[id]);
            }
        }
    }
    break;

    default : //never been reached
        break;
    }

    if (m_result.size() == 0) return false;

    m_result.setGeometry(getGeometry());
//     m_result.setName(m_field.getName());

    return true;

}

}
