#include <SFML/Graphics.hpp>

#include "pcg32.h"
#include "globals.h"
#include "NaiveEstimator.h"

void GenerateNDF(sf::Image& image, NaiveEstimator& estimator, int width, int height)
{
	image.create(width, height, sf::Color::Black);
	pcg32 random(0);

	int sampleCount = 8;

	float * data = new float[width * height];

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			float ndf = 0.f;

			for (int sample = 0; sample < sampleCount; sample++)
			{
				float fX = x + random.nextFloat();
				float fY = y + random.nextFloat();

				float s = (fX / width) * 2.f - 1.f;
				float t = (fY / height) * 2.f - 1.f;
				float v2 = 1.f - (s * s) - (t * t);

				if (v2 > 0.f)
				{
					float v = std::sqrt(v2);
					Vector3f w = glm::normalize(Vector3f(s, t, v));
					ndf += estimator.Estimate(w);
				}
			}

			data[y * width + x] = ndf / sampleCount;
		}
	}

	float maxValue = 0.f;

	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++)
			maxValue = glm::max(data[y * width + x], maxValue);

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			int g = glm::clamp(data[y  * width + x] / maxValue, 0.f, 1.f) * 255;
			image.setPixel(x, y, sf::Color(g, g, g));
		}
	}

	delete[] data;
}

int main()
{
	int width = 512;
	int height = 512;
	sf::RenderWindow window(sf::VideoMode(width, height), "NDF Estimator");
	sf::RectangleShape background(sf::Vector2f(width, height));

	NaiveEstimator estimator("normal5.png");

	sf::Image image;
	GenerateNDF(image, estimator, width, height);

	image.saveToFile("output.png");

	sf::Texture tx;
	tx.loadFromImage(image);

	background.setTexture(&tx);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear();
		window.draw(background);
		window.display();
	}

	return 0;
}