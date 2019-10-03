#pragma once

#include "ScenarioConnectomeAbstractStrategy.h"

class ScenarioConnectomeStrategyCUDANewton
	: public ScenarioConnectomeStrategyCUDA
{
	public:

		using ScenarioConnectomeStrategyCUDA::ScenarioConnectomeStrategyCUDA;

    void derivatives() override;

};
