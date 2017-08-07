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
#include "GenericOutput.hpp"

using namespace std;
namespace mimmo{

/*!
 * Default constructor of GenericOutput
 * \param[in] dir Directory of the output file (default value = "./").
 * \param[in] filename Name of the output file (default value = "output.txt").
 * \param[in] csv True if the output file is a csv format file (default value false).
 */
GenericOutput::GenericOutput(std::string dir, std::string filename, bool csv){
    m_dir       = dir;
    m_filename  = filename;
    m_csv       = csv;
    m_portsType    = ConnectionType::BACKWARD;
    m_name         = "mimmo.GenericOutput";
    m_binary        = false;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
GenericOutput::GenericOutput(const bitpit::Config::Section & rootXML){

    m_dir       = "./";
    m_filename  = "output.txt";
    m_portsType    = ConnectionType::BACKWARD;
    m_name         = "mimmo.GenericOutput";
    m_csv       = false;
    m_binary        = false;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.GenericOutput"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor of GenericOutput.
 */
GenericOutput::~GenericOutput(){
    m_portsType = ConnectionType::BACKWARD;
};

/*!
 * Copy constructor of GenericOutput. m_input and m_result member are not copied.
 */
GenericOutput::GenericOutput(const GenericOutput & other):BaseManipulation(other){
    m_dir           = other.m_dir;
    m_filename      = other.m_filename;
    m_csv           = other.m_csv;
    m_binary        = other.m_binary;
};

/*!
 * Assignment operator
 */
GenericOutput & GenericOutput::operator=(GenericOutput other){
    swap(other);
    return *this;
}

/*!
 * Swap method
 * \param[in] x object to be swapped
 */
void GenericOutput::swap(GenericOutput & x) noexcept
{
    std::swap(m_dir, x.m_dir);
    std::swap(m_filename, x.m_filename);
    std::swap(m_csv, x.m_csv);
    std::swap(m_binary, x.m_binary);
    std::swap(m_input, x.m_input);
    std::swap(m_result, x.m_result);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
GenericOutput::buildPorts(){

    PortManager::instance().addPort(M_COORDS, MC_VECARR3, MD_FLOAT);
    PortManager::instance().addPort(M_DISPLS, MC_VECARR3, MD_FLOAT);
    PortManager::instance().addPort(M_FILTER, MC_VECTOR, MD_FLOAT);
    PortManager::instance().addPort(M_POINT,  MC_ARRAY3, MD_FLOAT);
    PortManager::instance().addPort(M_DIMENSION, MC_ARRAY3, MD_INT);
    PortManager::instance().addPort(M_VALUED, MC_SCALAR, MD_FLOAT);
    PortManager::instance().addPort(M_VALUEI, MC_SCALAR, MD_INT);
    PortManager::instance().addPort(M_VALUEB, MC_SCALAR, MD_BOOL);
    PortManager::instance().addPort(M_DEG, MC_SCALAR, MD_INT);
    PortManager::instance().addPort(M_FILENAME, MC_SCALAR, MD_STRING);
    PortManager::instance().addPort(M_DIR, MC_SCALAR, MD_STRING);
    
    
    bool built = true;
    built = (built && createPortIn<dvecarr3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvecarr3E>, PortType::M_COORDS, true, 1));
    built = (built && createPortIn<dvecarr3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvecarr3E>, PortType::M_DISPLS, true, 1));
    built = (built && createPortIn<dvector1D, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvector1D>, PortType::M_FILTER, true, 1));
    built = (built && createPortIn<darray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<darray3E>, PortType::M_POINT, true, 1));
    built = (built && createPortIn<iarray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<iarray3E>, PortType::M_DIMENSION, true, 1));
    built = (built && createPortIn<double, GenericOutput>(this, &mimmo::GenericOutput::setInput<double>, PortType::M_VALUED, true, 1));
    built = (built && createPortIn<int, GenericOutput>(this, &mimmo::GenericOutput::setInput<int>, PortType::M_VALUEI, true, 1));
    built = (built && createPortIn<bool, GenericOutput>(this, &mimmo::GenericOutput::setInput<bool>, PortType::M_VALUEB, true, 1));
    built = (built && createPortIn<iarray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<iarray3E>, PortType::M_DEG, true, 1));
    built = (built && createPortIn<dmpvector1D, GenericOutput>(this, &mimmo::GenericOutput::setResult<dmpvector1D>, PortType::M_SCALARFIELD));
    built = (built && createPortIn<dmpvecarr3E, GenericOutput>(this, &mimmo::GenericOutput::setResult<dmpvecarr3E>, PortType::M_VECTORFIELD));
//    built = (built && createPortIn<string, GenericOutput>(this, &mimmo::GenericOutput::setInput<string>, PortType::M_FILENAME, true, 1));
//    built = (built && createPortIn<string, GenericOutput>(this, &mimmo::GenericOutput::setInput<string>, PortType::M_DIR, true, 1));

    m_arePortsBuilt = built;
}

/*!
 * It sets the name of the directory to write the output file.
 * \param[in] dir Name of the directory of the output file.
 */
void
GenericOutput::setWriteDir(std::string dir){
    m_dir = dir;
};

/*!
 * It sets the name of the output file.
 * \param[in] filename Name of the output file.
 */
void
GenericOutput::setFilename(std::string filename){
    m_filename = filename;
};

/*!
 * It sets if the output file has to be written in csv format.
 * \param[in] csv Write the output file in csv format.
 */
void
GenericOutput::setCSV(bool csv){
    m_csv = csv;
};

/*!It sets if the output file is in Binary format.
 * \param[in] binary Is the output file in Binary format?
 */
void
GenericOutput::setBinary(bool binary){
    m_binary = binary;
};


/*!
 * It clear the input member of the object
 */
void
GenericOutput::clearInput(){
    m_input.reset(nullptr);
}

/*!
 * It clear the result member of the object
 */
void
GenericOutput::clearResult(){
    m_result.reset(nullptr);
}

/*!Execution command.
 * It does nothing, the real execution of the object
 * happens in setInput.
 */
void
GenericOutput::execute(){

    bool ok = rename("output.txt", m_filename.c_str());
    if (!ok)
        (*m_log)<<"Generic Output execution failed renaming..."<<std::endl;
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
GenericOutput::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Filename")){
        std::string input = slotXML.get("Filename");
        input = bitpit::utils::string::trim(input);
        setFilename(input);
    };

    if(slotXML.hasOption("WriteDir")){
        std::string input = slotXML.get("WriteDir");
        input = bitpit::utils::string::trim(input);
        setWriteDir(input);
    };

    if(slotXML.hasOption("CSV")){
        std::string input = slotXML.get("CSV");
        input = bitpit::utils::string::trim(input);
        bool temp = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setCSV(temp);
    };

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
GenericOutput::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    
    slotXML.set("Filename", m_filename);
    slotXML.set("WriteDir", m_dir);
    slotXML.set("CSV", std::to_string((int)m_csv));
};


//specializations of setResult

/*!
 * Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Data to be written and to be used to set the input/result.
 */
template<>
void
GenericOutput::setResult(dmpvector1D data){
    _setInput(data);
    _setResult(data);
    std::fstream file;
    file.open(m_dir+"/"+m_filename);
    if (file.is_open()){
        if (m_binary){
            bitpit::genericIO::flushBINARY(file, data, true);
        }
        else{
            bitpit::genericIO::flushASCII(file, 1, data, true);
        }
        file.close();
    }
}

/*!
 * Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Data to be written and to be used to set the input/result.
 */
template<>
void
GenericOutput::setResult(dmpvecarr3E data){
    _setInput(data);
    _setResult(data);
    std::fstream file;
    file.open(m_dir+"/"+m_filename);
    if (file.is_open()){
        if (m_binary){
            bitpit::genericIO::flushBINARY(file, data, true);
        }
        else{
            bitpit::genericIO::flushASCII(file, 3, data, true);
        }
        file.close();
    }
}

}
