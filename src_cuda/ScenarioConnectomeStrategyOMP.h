#pragma once

#include "ScenarioConnectomeAbstractStrategy.h"

class ScenarioConnectomeStrategyOMP
	: public ScenarioConnectomeAbstractStrategy
{
	public:

		using ScenarioConnectomeAbstractStrategy::ScenarioConnectomeAbstractStrategy;
		// ScenarioConnectomeStrategyCPU(ScenarioConnectome* pthis)
		// 	: ScenarioConnectomeAbstractStrategy(pthis) {}

	void update() override;
	void derivatives() override;
};
