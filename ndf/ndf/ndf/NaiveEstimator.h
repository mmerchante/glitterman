#pragma once
#include "NDFEstimator.h"
#include "globals.h"
#include <string>

class NaiveEstimator : public NDFEstimator
{
public:
	NaiveEstimator(std::string normalMapFilename);
	virtual ~NaiveEstimator();

	virtual float Estimate(const glm::vec3& w);

private:

	void Initialize(std::string normalMapFilename);

	Vector3f * normalMap;
	int normalMapWidth;
	int normalMapHeight;

	Vector2f from;
	Vector2f to;
};

