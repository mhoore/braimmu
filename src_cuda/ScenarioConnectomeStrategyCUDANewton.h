#pragma once

#include "ScenarioConnectomeAbstractStrategy.h"

class ScenarioConnectomeStrategyCUDANewton
	: public ScenarioConnectomeStrategyCUDA
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

		using ScenarioConnectomeStrategyCUDA::ScenarioConnectomeStrategyCUDA;

    void derivatives() override;

};
