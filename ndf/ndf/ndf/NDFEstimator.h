#pragma once

#include "globals.h"

class NDFEstimator
{
public:
	NDFEstimator() {}
	virtual ~NDFEstimator() {}
	virtual float Estimate(const glm::vec3& w) = 0;
};

