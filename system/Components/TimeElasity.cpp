#include "Component.h"

class TimeElasity : public Component{
	private:
	double lastValue;
	public:
		TimeElasity(double last){
			lastValue = last;
		}

		double* include(double param){
			double dt = ModelingTime::getStep() * ModelingTime::getStep();
			matrix[0][0] = 1./ dt ; matrix[0][1] = 0; matrix[0][2] = (2 * param - lastValue) / dt;
			matrix[1][0] = 0 ; matrix[1][1] = 1./dt ; matrix[1][2] =  (2 * param - lastValue) / dt;
			lastValue = param;
			return &matrix[0][0];
		}
};