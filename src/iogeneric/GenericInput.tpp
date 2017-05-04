#include <fstream>
#include "Operators.hpp"

namespace mimmo{
/*!Overloaded function of base class setInput.
 * It sets the input of the object, but at the same time it sets even the result.
 * \param[in] data Pointer to data to be used to set the input/result.
 */
template<typename T>
void
GenericInput::setInput(T* data){
    _setInput(data);
    _setResult(data);
}

/*!Overloaded function of base class setInput.
 * It sets the input of the object, but at the same time it sets even the result.
 * \param[in] data Data to be used to set the input/result.
 */
template<typename T>
void
GenericInput::setInput(T& data){
    _setInput(data);
    _setResult(data);
}

/*!Overloaded function of base class getResult.
 * It gets the result of the object, equal to the input.
 * In the case it reads the input from file before to set and to get the result.
 * \return Pointer to data stored in result member.
 */
template<typename T>
T
GenericInput::getResult(){
    if (m_readFromFile){
        T data;
        std::ifstream file;
        file.open(m_filename);
        if (file.is_open()){
            if (m_csv){
                ifstreamcsv(file, data);
            }
            else{
                file >> data;
            }
            file.close();
        }else{
            std::cout << "file not open --> exit" << std::endl;
            exit(1);
        }
        _setResult(data);
    }
    T temp = (*static_cast<IODataT<T>*>(m_result.get())->getData());

    return(temp);
}


// OLD BASEMANIPULATION CLASS TEMPLATED INPUT/RESULT METHODS	//


/*!It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
GenericInput::getInput(){
    return(static_cast<IODataT<T>*>(m_input.get())->getData());
}

/*!It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericInput::setResult(T* data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
    m_result = std::move(dummy);
}

/*!It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericInput::setResult(T& data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(data));
    m_result = std::move(dummy);
}

//==================================================   //
// OLD BASEMANIPULATION CLASS TEMPLATED INPUT METHODS  //
//==================================================   //

/*!It sets the input member of the object.
 * \param[in] data Pointer to data to be stored in the input member.
 */
template<typename T>
void
GenericInput::_setInput(T* data){
    clearInput();
    std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
    m_input = std::move(dummy);
}

/*!It sets the input member of the object.
 * \param[in] data Data to be stored in the input member.
 */
template<typename T>
void
GenericInput::_setInput(T& data){
    clearInput();
    std::unique_ptr<IOData> dummy(new IODataT<T>(data));
    m_input = std::move(dummy);
}

/*!It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
GenericInput::_getInput(){
    return(static_cast<IODataT<T>*>(m_input.get())->getData());
}

//==================================================   //
// OLD BASEMANIPULATION CLASS TEMPLATED RESULT METHODS //
//==================================================   //

/*!It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericInput::_setResult(T* data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
    m_result = std::move(dummy);
}

/*!It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericInput::_setResult(T& data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(data));
    m_result = std::move(dummy);
}

/*!It gets the result member of the object.
 * \return Pointer to data stored in the result member.
 */
template<typename T>
T*
GenericInput::_getResult(){
    return(static_cast<IODataT<T>*>(m_result.get())->getData());
}


template <typename T>
std::ifstream&  GenericInput::ifstreamcsv(std::ifstream &in, T &x){

T                   dummy;
char                delim;

if ((in.good())) {
    if (in >> dummy && in >> delim) { x = dummy;}
}
return(in);
}

template <typename T>
std::ifstream&  GenericInput::ifstreamcsvend(std::ifstream &in, T &x){

T                   dummy;
char                delim;

if ((in.good())) {
    if (in >> dummy) { x = dummy;}
}
return(in);
}

template <typename T>
std::ifstream&  GenericInput::ifstreamcsv(std::ifstream &in, std::vector< T > &x){

T       dummy;
char    delim;
int     i, n = x.size();

//i = 0;
//while ((in.good()) && (i < n)) {
//    if (ifstreamcsv(in,dummy)) { x[i] = dummy; i++; }
//} //next i
while (in.good()) {
    if (ifstreamcsv(in,dummy)) { x.push_back(dummy); }
}
return(in);
}

template <typename T, size_t d>
std::ifstream&  GenericInput::ifstreamcsv(std::ifstream &in, std::array< T, d> &x){

T       dummy;
char    delim;
int     i;

i = 0;
while ((in.good()) && (i < d-1)) {
    if (ifstreamcsv(in,dummy)) { x[i] = dummy;}
    i++;
} //next i
if (ifstreamcsvend(in,dummy)) { x[i] = dummy;}
return(in);
};





}
