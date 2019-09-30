#pragma once

class ScenarioConnectome;

class ScenarioConnectomeAbstractStrategy
{
	protected:
	ScenarioConnectome* m_this;

	public:

		ScenarioConnectomeAbstractStrategy(ScenarioConnectome* pthis)
			: m_this(pthis) {}
	
	virtual void update() = 0;
	virtual void derivatives() = 0;

	virtual ~ScenarioConnectomeAbstractStrategy() {}
};
