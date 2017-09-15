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

#include "GenericInput.hpp"
#include "Operators.hpp"
#include <fstream>

using namespace std;
namespace mimmo {

////////GENERICINPUT///////////////////////////////////////
    
/*!
 * Default constructor of GenericInput.
 * \param[in] readFromFile True if the object reads the values from file (default value false).
 * \param[in] csv True if the input file is a csv format file (default value false).
 */
GenericInput::GenericInput(bool readFromFile, bool csv){
    m_readFromFile  = readFromFile;
    m_csv           = csv;
    m_portsType     = BaseManipulation::ConnectionType::BOTH;
    m_name          = "mimmo.GenericInput";
    m_dir           = "./";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
GenericInput::GenericInput(const bitpit::Config::Section & rootXML){

    m_readFromFile  = false;
    m_csv           = false;
    m_portsType     = BaseManipulation::ConnectionType::BOTH;
    m_name             = "mimmo.GenericInput";
    m_dir       = "./";
    m_filename  = "input.txt";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.GenericInput"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom constructor of GenericInput.
 * \param[in] dir Directory of the input file.
 * \param[in] filename Name of the input file.
 * \param[in] csv True if the input file is a csv format file (default value false).
 * The m_readFromFile flag is set to true.
 */
GenericInput::GenericInput(std::string dir, std::string filename, bool csv){
    m_name          = "mimmo.GenericInput";
    m_readFromFile  = true;
    m_csv           = csv;
    m_dir           = dir;
    m_filename      = filename;
    m_portsType     = BaseManipulation::ConnectionType::BOTH;
};

GenericInput::~GenericInput(){};

/*!
 * Copy constructor of GenericInput. m_input and  m_result members are not copied.
 */
GenericInput::GenericInput(const GenericInput & other):BaseManipulation(other){
    m_readFromFile  = other.m_readFromFile;
    m_csv           = other.m_csv;
    m_dir           = other.m_dir;
    m_filename      = other.m_filename;
};

/*!
 * Assignment operator. m_input and  m_result members are not copied. 
 */
GenericInput & GenericInput::operator=(GenericInput other){
    swap(other);
    return *this;
}

/*!
 *Swap method.
 *\param[in] x objet to be swapped.
 */
void GenericInput::swap(GenericInput & x) noexcept
{
    std::swap(m_readFromFile, x.m_readFromFile);
    std::swap(m_csv         , x.m_csv);
    std::swap(m_dir         , x.m_dir);
    std::swap(m_filename    , x.m_filename);
    std::swap(m_input       , x.m_input);
    std::swap(m_result      , x.m_result);
    BaseManipulation::swap(x);
}

/*!
 * It sets if the object imports the displacements from an input file.
 * \param[in] readFromFile True if the object reads the values from file.
 */
void
GenericInput::setReadFromFile(bool readFromFile){
    m_readFromFile = readFromFile;
};

/*!It sets the name of the input file.
 * \param[in] filename Name of the input file.
 */
void
GenericInput::setFilename(std::string filename){
    m_filename = filename;
};

/*!It sets the name of the directory of input file.
 * \param[in] dir Name of the directory of input file.
 */
void
GenericInput::setReadDir(std::string dir){
    m_dir = dir;
};

/*! It builds the input/output ports of the object
 */
void
GenericInput::buildPorts(){

    bool built = true;

    built = (built && createPortOut<dvecarr3E, GenericInput>(this, &mimmo::GenericInput::getResult<dvecarr3E>, M_COORDS));
    built = (built && createPortOut<dvecarr3E, GenericInput>(this, &mimmo::GenericInput::getResult<dvecarr3E>, M_DISPLS));
    built = (built && createPortOut<dvector1D, GenericInput>(this, &mimmo::GenericInput::getResult<dvector1D>, M_DATAFIELD));
    built = (built && createPortOut<darray3E, GenericInput>(this, &mimmo::GenericInput::getResult<darray3E>, M_POINT));
    built = (built && createPortOut<darray3E, GenericInput>(this, &mimmo::GenericInput::getResult<darray3E>, M_SPAN));
    built = (built && createPortOut<iarray3E, GenericInput>(this, &mimmo::GenericInput::getResult<iarray3E>, M_DIMENSION));;
    built = (built && createPortOut<double, GenericInput>(this, &mimmo::GenericInput::getResult<double>, M_VALUED));
    built = (built && createPortOut<int, GenericInput>(this, &mimmo::GenericInput::getResult<int>, M_VALUEI));
    built = (built && createPortOut<bool, GenericInput>(this, &mimmo::GenericInput::getResult<bool>, M_VALUEB));
    built = (built && createPortOut<iarray3E, GenericInput>(this, &mimmo::GenericInput::getResult<iarray3E>, M_DEG));

    m_arePortsBuilt = built;
}

/*!It clear the input member of the object
 */
void
GenericInput::clearInput(){
    m_input.reset(nullptr);
}

/*!It clear the result member of the object
 */
void
GenericInput::clearResult(){
    m_result.reset(nullptr);
}

/*!Execution command.
 * It does nothing, the real execution of the object
 * happens in setInput/getResult.
 */
void
GenericInput::execute(){};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
GenericInput::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("ReadFromFile")){
        std::string input = slotXML.get("ReadFromFile");
        input = bitpit::utils::string::trim(input);
        bool temp = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setReadFromFile(temp);
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

    if(slotXML.hasOption("Filename") && m_readFromFile){
        std::string input = slotXML.get("Filename");
        input = bitpit::utils::string::trim(input);
        setFilename(input);
    };

    if(slotXML.hasOption("ReadDir") && m_readFromFile){
        std::string input = slotXML.get("ReadDir");
        input = bitpit::utils::string::trim(input);
        setReadDir(input);
    };

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
GenericInput::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    
    slotXML.set("ReadFromFile", std::to_string((int)m_readFromFile));
    slotXML.set("CSV", std::to_string((int)m_csv));
    slotXML.set("ReadDir", m_dir);
    slotXML.set("Filename", m_filename);
};


/*!
 * Default constructor of GenericInputMPVData.
 * \param[in] readFromFile True if the object reads the values from file (default value false).
 * \param[in] csv True if the input file is a csv format file (default value false).
 */
GenericInputMPVData::GenericInputMPVData(bool csv){
    m_csv           = csv;
    m_portsType     = BaseManipulation::ConnectionType::BOTH;
    m_name          = "mimmo.GenericInputMPVData";
    m_dir           = "./";
    m_binary        = false;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
GenericInputMPVData::GenericInputMPVData(const bitpit::Config::Section & rootXML){
    
    m_csv           = false;
    m_portsType     = BaseManipulation::ConnectionType::BOTH;
    m_name             = "mimmo.GenericInputMPVData";
    m_dir       = "./";
    m_filename  = "input.txt";
    m_binary        = false;
    
    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.GenericInputMPVData"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom constructor of GenericInputMPVData.
 * \param[in] dir Directory of the input file.
 * \param[in] filename Name of the input file.
 * \param[in] csv True if the input file is a csv format file (default value false).
 * The m_readFromFile flag is set to true.
 */
GenericInputMPVData::GenericInputMPVData(std::string dir, std::string filename, bool csv){
    m_csv           = csv;
    m_dir           = dir;
    m_filename      = filename;
    m_name          = "mimmo.GenericInputMPVData";
    m_binary        = false;
    m_portsType     = BaseManipulation::ConnectionType::BOTH;
};

GenericInputMPVData::~GenericInputMPVData(){};

/*!
 * Copy constructor of GenericInputMPVData. m_input and  m_result members are not copied.
 */
GenericInputMPVData::GenericInputMPVData(const GenericInputMPVData & other):BaseManipulation(other){
    m_csv           = other.m_csv;
    m_dir           = other.m_dir;
    m_filename      = other.m_filename;
    m_binary        = other.m_binary;
};

/*!
 * Assignment operator. m_input and  m_result members are not copied. 
 */
GenericInputMPVData & GenericInputMPVData::operator=(GenericInputMPVData other){
    swap(other);
    return *this;
}

/*!
 * Swap method.
 *\param[in] x objet to be swapped.
 */
void GenericInputMPVData::swap(GenericInputMPVData & x) noexcept
{
    std::swap(m_csv         , x.m_csv);
    std::swap(m_dir         , x.m_dir);
    std::swap(m_filename    , x.m_filename);
    std::swap(m_binary      , x.m_binary);
    std::swap(m_input       , x.m_input);
    std::swap(m_result      , x.m_result);
    BaseManipulation::swap(x);
}


/*!It sets if the input file is in csv format.
 * \param[in] csv Is the input file write in comma separated value format?
 */
void
GenericInputMPVData::setCSV(bool csv){
    m_csv = csv;
};

/*!It sets if the input file is in Binary format.
 * \param[in] binary Is the input file in Binary format?
 */
void
GenericInputMPVData::setBinary(bool binary){
    m_binary = binary;
};

/*!It sets the name of the input file.
 * \param[in] filename Name of the input file.
 */
void
GenericInputMPVData::setFilename(std::string filename){
    m_filename = filename;
};

/*!It sets the name of the directory of input file.
 * \param[in] dir Name of the directory of input file.
 */
void
GenericInputMPVData::setReadDir(std::string dir){
    m_dir = dir;
};

/*! It builds the input/output ports of the object
 */
void
GenericInputMPVData::buildPorts(){
    
    bool built = true;
    
    built = (built && createPortIn<MimmoObject*, GenericInputMPVData>(this, &mimmo::GenericInputMPVData::setGeometry, M_GEOM, true));
    built = (built && createPortOut<dmpvector1D, GenericInputMPVData>(this, &mimmo::GenericInputMPVData::getResult<dmpvector1D>, M_SCALARFIELD));
    built = (built && createPortOut<dmpvecarr3E, GenericInputMPVData>(this, &mimmo::GenericInputMPVData::getResult<dmpvecarr3E>, M_VECTORFIELD));
    
    
    m_arePortsBuilt = built;
}

/*!It clear the input member of the object
 */
void
GenericInputMPVData::clearInput(){
    m_input.reset(nullptr);
}

/*!It clear the result member of the object
 */
void
GenericInputMPVData::clearResult(){
    m_result.reset(nullptr);
}

/*!Execution command.
 * It does nothing, the real execution of the object
 * happens in setInput/getResult.
 */
void
GenericInputMPVData::execute(){};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
GenericInputMPVData::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
    
    BITPIT_UNUSED(name);
    
    std::string input;
    
    BaseManipulation::absorbSectionXML(slotXML, name);
    
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
    
    if(slotXML.hasOption("Filename") && m_readFromFile){
        std::string input = slotXML.get("Filename");
        input = bitpit::utils::string::trim(input);
        setFilename(input);
    };
    
    if(slotXML.hasOption("ReadDir") && m_readFromFile){
        std::string input = slotXML.get("ReadDir");
        input = bitpit::utils::string::trim(input);
        setReadDir(input);
    };
    
    if(slotXML.hasOption("Binary")){
        std::string input = slotXML.get("Binary");
        input = bitpit::utils::string::trim(input);
        bool temp = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setBinary(temp);
    };
    
}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
GenericInputMPVData::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
    
    BITPIT_UNUSED(name);
    
    BaseManipulation::flushSectionXML(slotXML, name);
    
    slotXML.set("CSV", std::to_string((int)m_csv));
    slotXML.set("ReadDir", m_dir);
    slotXML.set("Filename", m_filename);
    slotXML.set("Binary", std::to_string((int)m_binary));
};



//specializations of getResult
/*!
 * Overloaded function of base class getResult.
 * It gets the result of the object, equal to the input.
 * In the case it reads the input from file before to set and to get the result.
 * \return Pointer to data stored in result member.
 */
template <>
dmpvector1D
GenericInputMPVData::getResult(){
    dmpvector1D data;
    if (getGeometry() == NULL) return data;
    int nEle;
    if(m_csv)   m_binary = false;
    std::fstream file;
    file.open(m_dir+"/"+m_filename);
        bitpit::PiercedVector<double> pvdata;
        if (file.is_open()){
            if (m_binary){
                bitpit::genericIO::absorbBINARY(file, nEle);
                bitpit::genericIO::absorbBINARY(file, pvdata, nEle);
            }
            else if (m_csv){
                ifstreamcsv(file, data);
            }else{
                bitpit::genericIO::absorbASCII(file, nEle);
                bitpit::genericIO::absorbASCII(file, pvdata, nEle);
            }
            file.close();
        }else{
            (*m_log) << "file not open --> exit" << std::endl;
            throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested");
        }
        
        data->setGeometry(getGeometry());
        if(!m_csv){
            data = pvdata;
            data->setDataLocation(data.recoverGeometryReferenceLocation());
        }    
        if(data.checkDataIdsCoherence()) _setResult(data);
    
    return *static_cast<IODataT<dmpvector1D>*>(m_result.get())->getData();
}


//specializations of getResult
/*!
 * Overloaded function of base class getResult.
 * It gets the result of the object, equal to the input.
 * In the case it reads the input from file before to set and to get the result.
 * \return Pointer to data stored in result member.
 */
template <>
dmpvecarr3E
GenericInputMPVData::getResult(){
    dmpvecarr3E data;
    if (getGeometry() == NULL) return data;
    int nEle;
    if(m_csv)   m_binary = false;
    std::fstream file;
    file.open(m_dir+"/"+m_filename);
    bitpit::PiercedVector<std::array<double,3> > pvdata;
    if (file.is_open()){
        if (m_binary){
            bitpit::genericIO::absorbBINARY(file, nEle);
            bitpit::genericIO::absorbBINARY(file, pvdata, nEle);
        }
        else if (m_csv){
            ifstreamcsv(file, data);
        }else{
            bitpit::genericIO::absorbASCII(file, nEle);
            bitpit::genericIO::absorbASCII(file, pvdata, nEle);
        }
        file.close();
    }else{
        (*m_log) << "file not open --> exit" << std::endl;
        throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested");
    }
    
    data->setGeometry(getGeometry());
    if(!m_csv){
        data = pvdata;
        data->setDataLocation(data.recoverGeometryReferenceLocation());
    }    
    if(data.checkDataIdsCoherence()) _setResult(data);

    return *static_cast<IODataT<dmpvecarr3E>*>(m_result.get())->getData();
}





}

