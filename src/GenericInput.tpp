#include <fstream>
#include "Operators.hpp"

/*!Overloaded function of base class setInput.
 * It sets the input of the object, but at the same time it sets even the result.
 * \param[in] data Pointer to data to be used to set the input/result.
 */
template<typename T>
void
GenericInput::setInput(T* data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
}

/*!Overloaded function of base class setInput.
 * It sets the input of the object, but at the same time it sets even the result.
 * \param[in] data Data to be used to set the input/result.
 */
template<typename T>
void
GenericInput::setInput(T& data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
}

/*!Overloaded function of base class getResult.
 * It gets the result of the object, equal to the input.
 * In the case it reads the input from file before to set and to get the result.
 * \return Pointer to data stored in result member.
 */
template<typename T>
T*
GenericInput::getResult(){
	if (m_readFromFile){
		T data;
		std::ifstream file;
		file.open(m_filename);
		if (file.is_open()){
			file >> data;
			file.close();
		}else{
			std::cout << "file not open --> exit" << std::endl;
			exit(1);
		}
		BaseManipulation::setResult(data);
	}
	return(static_cast<IODataT<T>*>(m_result.get())->getData());
}




