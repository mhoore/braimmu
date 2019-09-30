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

		ScenarioConnectomeStrategyCPU(ScenarioConnectome* pthis)
		 	: ScenarioConnectomeAbstractStrategy(pthis) {}

	void update() override;
	void derivatives() override;
};
