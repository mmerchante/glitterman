#include "NaiveEstimator.h"
#include <SFML\Graphics.hpp>
#include "globals.h"

#define BINS_COUNT 8

// zero mean, small sigma
float Gaussian(float x)
{
	float sigma2 = .01f * .01f;
	return glm::exp(-x * x / (2.f * sigma2)) / glm::sqrt(2.f * 3.141592f * sigma2);
}

NaiveEstimator::NaiveEstimator(std::string normalMapFilename) : normalMap(nullptr), normalMapWidth(100), normalMapHeight(100)
{
	this->Initialize(normalMapFilename);
}

NaiveEstimator::~NaiveEstimator()
{
}

float NaiveEstimator::Estimate(const glm::vec3 & w)
{
	Vector2i from(0, 0);
	Vector2i to(normalMapWidth, normalMapHeight);

	int pixels = glm::abs((to.x - from.x) * (to.y - from.y));
	float ndf = 0.f;

	for (int x = from.x; x < to.x; x++)
	{
		for (int y = from.y; y < to.y; y++)
		{
			int i = (y * normalMapWidth) + x;
			Vector3f wh = this->normalMap[i];
/*
			float d = glm::clamp(glm::dot(wh, w), 0.f, 0.9999f);
			ndf += 1.f / (1.f - d);*/
			ndf += Gaussian(glm::length(wh - w));
		}
	}

	return ndf / pixels;
}

void NaiveEstimator::BinningMethod(float * pixels, int width, int height)
{
	Vector2i from(0, 0);
	Vector2i to(normalMapWidth, normalMapHeight);

	for (int x = from.x; x < to.x; x++)
	{
		for (int y = from.y; y < to.y; y++)
		{
			int i = (y * normalMapWidth) + x;
			Vector3f wh = this->normalMap[i];

			int pX = (wh.x * .5f + .5f) * width;
			int pY = (wh.y * .5f + .5f) * width;

			pixels[pY * width + pX]++;
		}
	}
}

void NaiveEstimator::Initialize(std::string normalMapFilename)
{
	sf::Image normalMapImage;

	if (normalMapImage.loadFromFile(normalMapFilename))
	{
		this->normalMapWidth = normalMapImage.getSize().x;
		this->normalMapHeight = normalMapImage.getSize().y;

		int pixels = normalMapHeight * normalMapWidth;
		this->normalMap = new Vector3f[pixels];

		const sf::Uint8 * rawImage = normalMapImage.getPixelsPtr();

		for (int p = 0; p < pixels; p++)
		{
			sf::Uint8 r = rawImage[p * 4];
			sf::Uint8 g = rawImage[p * 4 + 1];
			sf::Uint8 b = rawImage[p * 4 + 2];

			this->normalMap[p] = glm::normalize(Vector3f(r / 255.f, g / 255.f, b / 255.f) * 2.f - Vector3f(1.f));
		}
	}
}