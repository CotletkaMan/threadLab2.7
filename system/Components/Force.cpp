#include "Component.h"

class Force : public Component{
	private:
	double x;
	double y;
	double (*forceFunction)(double, double, double);
	public:
		Force(double x , double y , double (*function)(double, double, double)){
			this -> x = x;
			this -> y = y;
			forceFunction = function;
		}

		double* include(double param){
			matrix[0][0] = 0 ; matrix[0][1] = 0 ; matrix[0][2] = -forceFunction(x , y , ModelingTime::getTime() + ModelingTime::getStep());
			matrix[1][0] = 0 ; matrix[1][1] = 0 ; matrix[1][2] = forceFunction(x , y , ModelingTime::getTime() + ModelingTime::getStep());
			return &matrix[0][0];
		}
};