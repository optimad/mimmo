
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
