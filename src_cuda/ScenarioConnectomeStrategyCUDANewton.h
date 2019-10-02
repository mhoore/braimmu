#pragma once

#include "ScenarioConnectomeAbstractStrategy.h"

class ScenarioConnectomeStrategyCUDANewton
	: public ScenarioConnectomeAbstractStrategy
{

  public:

  struct array_properties {
    double *Dtau;
  };
  
  protected:
  array_properties arr_prop;
  double *agent, *deriv;
  int *type;

	public:

		ScenarioConnectomeStrategyCUDANewton(ScenarioConnectome* pthis);
		~ScenarioConnectomeStrategyCUDANewton() override;

	void update() override;
	void derivatives() override;

	void push() override;
	void pop() override;
};
