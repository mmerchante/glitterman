#include "NaiveEstimator.h"
#include <SFML\Graphics.hpp>
#include "globals.h"

#define BINS_COUNT 8

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
	Vector2i to(128, 128);

	int pixels = glm::abs((to.x - from.x) * (to.y - from.y));
	float ndf = 0.f;

	for (int x = from.x; x < to.x; x++)
	{
		for (int y = from.y; y < to.y; y++)
		{
			int i = (y * normalMapWidth) + x;
			Vector3f wh = this->normalMap[i];

			float d = glm::clamp(glm::dot(wh, w), 0.f, 0.9999f);
			ndf += 1.f / (1.f - d);
		}
	}

	return ndf / pixels;
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