#pragma once

#include "ScenarioConnectomeAbstractStrategy.h"

class ScenarioConnectomeStrategyCUDA
	: public ScenarioConnectomeAbstractStrategy
{

  public:

  struct array_properties {
    double *Dtau;
  };

  struct AllocPitch {
	  size_t pDouble, pInt;
  };
  
  protected:
  array_properties arr_prop;
  double *agent, *deriv;
  int *type;

  AllocPitch m_allocPitch;

	inline size_t nvx();
	inline size_t nvyz();

	public:

		ScenarioConnectomeStrategyCUDA(ScenarioConnectome* pthis);
		~ScenarioConnectomeStrategyCUDA() override;

	void update() override;
	void derivatives() override;

	void push() override;
	void pop() override;
};
