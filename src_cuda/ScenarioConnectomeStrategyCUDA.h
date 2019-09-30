#pragma once

#include "ScenarioConnectomeAbstractStrategy.h"

class ScenarioConnectomeStrategyCUDA
	: public ScenarioConnectomeAbstractStrategy
{
  protected:

  struct array_properties {
    double *Dtau;
  } arr_prop;

  double *agent, *deriv;
  int *type;

	public:

		ScenarioConnectomeStrategyCUDA(ScenarioConnectome* pthis);
		~ScenarioConnectomeStrategyCUDA() override;

	void update() override;
	void derivatives() override;

	void push() override;
	void pop() override;
};
