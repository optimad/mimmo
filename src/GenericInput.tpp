#include <fstream>
#include "Operators.hpp"

template<typename T>
void
GenericInput::setInput(T* data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
}

template<typename T>
void
GenericInput::setInput(T& data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
}


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
	return BaseManipulation::getResult<T>();
}

