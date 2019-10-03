#pragma once

#include "ScenarioConnectomeAbstractStrategy.h"

class ScenarioConnectomeStrategyCUDA
	: public ScenarioConnectomeAbstractStrategy
{

  public:

  struct array_properties {
    double *Dtau;
  };
  
  protected:
  array_properties arr_prop;
  double *agent, *agent2, *deriv;
  int *type;

	public:

		ScenarioConnectomeStrategyCUDA(ScenarioConnectome* pthis);
		~ScenarioConnectomeStrategyCUDA() override;

	void update() override;
	void derivatives() override;

	void push() override;
	void pop() override;
};
