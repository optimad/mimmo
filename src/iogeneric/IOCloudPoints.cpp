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

#include "IOCloudPoints.hpp"

namespace mimmo {

/*!
 * Default constructor of IOCloudPoints.
 * \param[in] readMode True if the object is in read mode, false if in Write mode.
 */
IOCloudPoints::IOCloudPoints(bool readMode){
    m_name         = "mimmo.IOCloudPoints";
    m_read         = readMode;
    m_template     = false;
    m_dir       = ".";
    m_filename     = m_name+"_source.dat";
    m_isInternal = true;
    m_intgeo.reset(nullptr);
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOCloudPoints::IOCloudPoints(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.IOCloudPoints";
    m_template     = false;
    m_dir       = ".";
    m_filename     = m_name+"_source.dat";
    m_read         = true;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);

    std::string fallback_name2 = "1";
    std::string input2 = rootXML.get("IOmode", fallback_name2);
    input2 = bitpit::utils::string::trim(input2);

    m_read = bool(std::atoi(input2.c_str()));

    if(input == "mimmo.IOCloudPoints"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };

    m_isInternal = true;
    m_intgeo.reset(nullptr);

}


IOCloudPoints::~IOCloudPoints(){};

/*!
 * Copy constructor of IOCloudPoints. Labels, points and data attached are no copied.
 */
IOCloudPoints::IOCloudPoints(const IOCloudPoints & other):BaseManipulation(other){
    m_read         = other.m_read;
    m_dir       = other.m_dir;
    m_filename     = other.m_filename;
    m_template = other.m_template;

    if(other.m_isInternal){
        m_geometry = other.m_intgeo.get();
    }
    m_isInternal = false;

};

/*!Assignement operator of IOCloudPoints. Labels, points and data attached are no copied.
 */
IOCloudPoints & IOCloudPoints::operator=(IOCloudPoints other){
    swap(other);
    return *this;
};

/*!
 * Swap method
 * \param[in] x object to be swapped
 */
void IOCloudPoints::swap(IOCloudPoints & x) noexcept
{
    std::swap(m_read        , x.m_read);
    std::swap(m_dir         , x.m_dir);
    std::swap(m_filename    , x.m_filename);
    std::swap(m_template    , x.m_template);
    std::swap(m_labels      , x.m_labels);
    std::swap(m_points      , x.m_points);
    std::swap(m_scalarfield , x.m_scalarfield);
    std::swap(m_vectorfield , x.m_vectorfield);
    m_scalarfield.swap(x.m_scalarfield);
    m_vectorfield.swap(x.m_vectorfield);
    std::swap(m_isInternal, x.m_isInternal);
    std::swap(m_intgeo, x.m_intgeo);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
IOCloudPoints::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dmpvector1D*, IOCloudPoints>(this, &mimmo::IOCloudPoints::setScalarField, M_SCALARFIELD));
    built = (built && createPortIn<dmpvecarr3E*, IOCloudPoints>(this, &mimmo::IOCloudPoints::setVectorField, M_VECTORFIELD));
    built = (built && createPortIn<MimmoObject*, IOCloudPoints>(this, &mimmo::IOCloudPoints::setGeometry, M_GEOM));

    built = (built && createPortOut<MimmoObject*, IOCloudPoints>(this, &mimmo::IOCloudPoints::getGeometry, M_GEOM));
    built = (built && createPortOut<dmpvector1D*, IOCloudPoints>(this, &mimmo::IOCloudPoints::getScalarField, M_SCALARFIELD));
    built = (built && createPortOut<dmpvecarr3E*, IOCloudPoints>(this, &mimmo::IOCloudPoints::getVectorField, M_VECTORFIELD));

    m_arePortsBuilt = built;
}

/*!
 * It gets the point cloud as a linked MimmoObject.
 * \return Pointer to point cloud geometry.
 */
MimmoObject*
IOCloudPoints::getGeometry(){
    // Return internal geometry
    if (m_isInternal) return m_intgeo.get();
    return m_geometry;
};

/*!
 * Return the scalar field stored in the class as pointer to MimmoPiercedVector object.
 * \return scalar field stored in the class
 */
dmpvector1D*
IOCloudPoints::getScalarField(){
    return &m_scalarfield;
};

/*!
 * Return the vector field stored in the class as pointer to MimmoPiercedVector object.
 * \return vector field stored in the class
 */
dmpvecarr3E*
IOCloudPoints::getVectorField(){
    return &m_vectorfield;
};

/*!
 * Return if template option is active. This method is only meant in class working in Write mode
 * \return template flag
 */
bool
IOCloudPoints::isTemplate(){
    return m_template;
};

/*!
 * It sets the name of the input directory. Active only in read mode.
 * \param[in] dir directory path
 */
void
IOCloudPoints::setReadDir(std::string dir){
    if(!m_read)    return;
    m_dir = dir;
};

/*!
 * It sets the name of the input file. Active only in read mode.
 * \param[in] filename filename with tag extension included.
 */
void
IOCloudPoints::setReadFilename(std::string filename){
    if(!m_read)    return;
    m_filename = filename;
};

/*!
 * It sets the name of the output directory. Active only in write mode.
 * \param[in] dir directory path
 */
void
IOCloudPoints::setWriteDir(std::string dir){
    if(m_read)    return;
    m_dir = dir;
};

/*!
 * It sets the point cloud as a linked MimmoObject.
 * \param[in] geometry Pointer to point cloud geometry.
 */
void
IOCloudPoints::setGeometry(MimmoObject* geometry){

    // Check if external geometry is a point cloud
    if (geometry->getType() != 3) return;

    // Set external geometry
    m_isInternal = false;
    m_intgeo.reset(nullptr);
    m_geometry = geometry;

};

/*!
 * It sets the name of the output file. Active only in write mode.
 * \param[in] filename filename with tag extension included.
 */
void
IOCloudPoints::setWriteFilename(std::string filename){
    if(m_read)    return;
    m_filename = filename;
};

/*!
 * It sets the scalar field associated to point.
 * The method is not active in Read mode.
 * \param[in] scalarfield pointer to scalar field
 */
void
IOCloudPoints::setScalarField(dmpvector1D* scalarfield){
    if(m_read) return;
    m_scalarfield = *scalarfield;
};

/*!
 * It sets the vector field associated to point.
 * The method is not active in Read mode.
 * \param[in] vectorfield pointer to vector field
 */
void
IOCloudPoints::setVectorField(dmpvecarr3E* vectorfield){
    if(m_read) return;
    m_vectorfield = *vectorfield;
};


/*!
 * Enables the template writing mode. The method is not active in Read mode.
 * \param[in] flag true to enable the template writing
 */
void
IOCloudPoints::setTemplate(bool flag){
    if(m_read) return;
    m_template = flag;
};

/*!
 * Clear all data stored in the class
 */
void
IOCloudPoints::clear(){
    m_vectorfield.clear();
    m_scalarfield.clear();
    m_template     = false;
    m_filename     = m_name+"_source.dat";
    m_isInternal      = true;
    m_intgeo.reset(nullptr);
    BaseManipulation::clear();
}

/*!
 * Execution command.
 * Read data from or Write data on linked filename
 */
void
IOCloudPoints::execute(){
    if(m_read){
        read();
    }else{
        write();
    }
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOCloudPoints::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    //checking IOmode

    if(slotXML.hasOption("IOmode")){
        input = slotXML.get("IOmode");
        bool value = true;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss>>value;
        }
        if (m_read != value){
            (*m_log)<< "warning in class "<<m_name<<": cannot absorb xml parameters for class IOmode mismatching"<<std::endl;
            throw std::runtime_error (m_name + " : xml absorbing failed ");
        }
    };

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(m_read){
        if(slotXML.hasOption("ReadDir")){
            std::string input = slotXML.get("ReadDir");
            input = bitpit::utils::string::trim(input);
            setReadDir(input);
        };

        if(slotXML.hasOption("ReadFilename")){
            std::string input = slotXML.get("ReadFilename");
            input = bitpit::utils::string::trim(input);
            setReadFilename(input);
        };
    }else{
        if(slotXML.hasOption("WriteDir")){
            std::string input = slotXML.get("WriteDir");
            input = bitpit::utils::string::trim(input);
            setWriteDir(input);
        };

        if(slotXML.hasOption("WriteFilename")){
            std::string input = slotXML.get("WriteFilename");
            input = bitpit::utils::string::trim(input);
            setWriteFilename(input);
        };
    }


    if(slotXML.hasOption("Template")){
        std::string input = slotXML.get("Template");
        input = bitpit::utils::string::trim(input);
        bool temp = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setTemplate(temp);
    };

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOCloudPoints::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("IOmode", std::to_string(int(m_read)));

    if(m_read){
        slotXML.set("ReadDir", m_dir);
        slotXML.set("ReadFilename", m_filename);
    }else{
        slotXML.set("WriteDir", m_dir);
        slotXML.set("WriteFilename", m_filename);
    }

    slotXML.set("Template", std::to_string(int(m_template)));
};


/*!
 * Plot cloud point and store it in *.vtu file
 */
void
IOCloudPoints::plotOptionalResults(){

    if (getGeometry() == nullptr) return;
    //Assessing data;
    MimmoPiercedVector<double> f1 = m_scalarfield;
    MimmoPiercedVector<std::array<double,3>> f2 = m_vectorfield;

    if (!f1.getGeometry()) f1.setGeometry(getGeometry());
    if (!f2.getGeometry()) f2.setGeometry(getGeometry());
    //force both data to have location on POINT, otherwise all these class is pointless
    f1.setDataLocation(MPVLocation::POINT);
    f2.setDataLocation(MPVLocation::POINT);
    f1.completeMissingData(0.0);
    f2.completeMissingData({{0.0,0.0,0.0}});

    if(f1.getName() == f2.getName()){
        std::string temp = f1.getName();
        f1.setName(temp+"_scalar");

        temp = f2.getName();
        f2.setName(temp+"_vector");

    }
    //Write
    BaseManipulation::write(getGeometry(), f1, f2);
}

/*!
 * Read displacement data from file
 */
void
IOCloudPoints::read(){

    std::unordered_map<long, int> mapP;
    std::ifstream reading;
    std::string source = m_dir+"/"+m_filename;
    reading.open(source.c_str());
    if(reading.is_open()){

        m_points.clear();
        m_labels.clear();
        m_scalarfield.clear();
        m_vectorfield.clear();

        std::string keyword, line;
        long label;
        darray3E dtrial;
        double temp;

        do{
            line.clear();
            keyword.clear();
            std::getline(reading,line);
            line = bitpit::utils::string::trim(line);
            std::stringstream ss(line);
            ss>>keyword;
            keyword = bitpit::utils::string::trim(keyword);
            if(keyword == "$POINT")    {

                ss>>label;
                m_labels.push_back(label);

                dtrial.fill(0.0);
                ss>>dtrial[0]>>dtrial[1]>>dtrial[2];
                m_points.push_back(dtrial);
            }

        }while(!reading.eof());

        int counter = 0;
        for(auto &lab :m_labels){
            mapP[lab] = counter;
            ++counter;
        }

        m_scalarfield.reserve(m_points.size());
        m_vectorfield.reserve(m_points.size());

        reading.clear();
        reading.seekg(0, std::ios::beg);

        do{
            line.clear();
            keyword.clear();
            std::getline(reading,line);
            line = bitpit::utils::string::trim(line);
            std::stringstream ss(line);
            ss>>keyword;
            keyword = bitpit::utils::string::trim(keyword);
            if(keyword == "$SCALARF")    {

                ss>>label>>temp;
                m_scalarfield.insert(mapP[label],temp);
            }
            if(keyword == "$VECTORF")    {

                dtrial.fill(0.0);
                ss>>label>>dtrial[0]>>dtrial[1]>>dtrial[2];
                m_vectorfield.insert(mapP[label],dtrial);
            }

        }while(!reading.eof());

    }else{
        (*m_log)<<"error of "<<m_name<<" : cannot open "<<m_filename<< " requested. Exiting... "<<std::endl;
        throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested. Exiting... ");
    }

    reading.close();


    // Fill geometry

    // Force external geometry pointer to null pointer
    m_geometry = nullptr;

    // Instantiate a new MimmoObject of type point cloud and fill it
    m_isInternal = true;
    m_intgeo.reset(nullptr);
    std::unique_ptr<MimmoObject> dum(new MimmoObject(3));
    m_intgeo = std::move(dum);

    // Fill geometry with point cloud list
    std::size_t ind = 0;
    for (darray3E & coords : m_points){
            long label = m_labels[ind];
            m_intgeo->addVertex(coords, label);
            m_intgeo->addConnectedCell(livector1D(1,label), bitpit::ElementType::VERTEX, label);
            ind++;
    }

    // Clear points and label structure
    dvecarr3E().swap(m_points);
    livector1D().swap(m_labels);

    // Set the mimmo pierced vectors
    m_scalarfield.setGeometry(m_intgeo.get());
    m_scalarfield.setDataLocation(MPVLocation::POINT);
    m_scalarfield.setName("scalarfield");
    m_vectorfield.setGeometry(m_intgeo.get());
    m_vectorfield.setDataLocation(MPVLocation::POINT);
    m_vectorfield.setName("vectorfield");

};

/*!
 * Write displacement data to file
 */
void
IOCloudPoints::write(){

    //assessing data;
    if (!m_scalarfield.getGeometry()) m_scalarfield.setGeometry(getGeometry());
    if (!m_vectorfield.getGeometry()) m_vectorfield.setGeometry(getGeometry());
    //force both data to have location on POINT, otherwise all these class is pointless
    m_scalarfield.setDataLocation(MPVLocation::POINT);
    m_vectorfield.setDataLocation(MPVLocation::POINT);

    if(!m_scalarfield.empty())  m_scalarfield.completeMissingData(0.0);
    if(!m_vectorfield.empty()) m_vectorfield.completeMissingData({{0.0,0.0,0.0}});

    std::string filename;
    if(m_template){
        filename = "TEMPLATE_"+m_filename;
    }else{
        filename = m_filename;
    }

    std::ofstream writing;
    std::string source = m_dir+"/"+m_filename;
    writing.open(source.c_str());
    std::string keyT1 = "{", keyT2 = "}";
    if(writing.is_open()){

        MimmoObject* geometry = getGeometry();
        darray3E coords;
        for(const long & label : geometry->getVerticesIds()){
            coords = geometry->getVertexCoords(label);
            writing<<"$POINT"<<'\t'<<label<<'\t'<<coords[0]<<'\t'<<coords[1]<<'\t'<<coords[2]<<std::endl;
        }
        writing<<""<<std::endl;

        for(const long & label : geometry->getVerticesIds()){
            if(m_template){
                std::string str1 = keyT1+"s"+std::to_string(label)+keyT2;
                writing<<"$SCALARF"<<'\t'<<label<<'\t'<<str1<<std::endl;
            }else{
                if(!m_scalarfield.exists(label)) continue;
                double val = m_scalarfield[label];
                writing<<"$SCALARF"<<'\t'<<label<<'\t'<<val<<std::endl;
            }
        }
        writing<<""<<std::endl;

        for(const long & label : geometry->getVerticesIds()){
            if(m_template){
                std::string str1 = keyT1+"x"+std::to_string(label)+keyT2;
                std::string str2 = keyT1+"y"+std::to_string(label)+keyT2;
                std::string str3 = keyT1+"z"+std::to_string(label)+keyT2;

                writing<<"$VECTORF"<<'\t'<<label<<'\t'<<str1<<'\t'<<str2<<'\t'<<str3<<std::endl;
            }else{
                if(!m_vectorfield.exists(label))    continue;
                std::array<double,3> val = m_vectorfield[label];
                writing<<"$VECTORF"<<'\t'<<label<<'\t'<<val[0]<<'\t'<<val[1]<<'\t'<<val[2]<<std::endl;
            }
        }
        writing<<""<<std::endl;
    }else{
        (*m_log)<<"error of "<<m_name<<" : cannot open "<<m_filename<< " requested. Exiting... "<<std::endl;
        throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested. Exiting... ");
    }

    writing.close();
};


}
