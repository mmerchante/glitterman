#include <SFML/Graphics.hpp>

#include "pcg32.h"
#include "globals.h"
#include "NaiveEstimator.h"

// zero mean, small sigma
float Gaussian1D(float x)
{
	float sigma2 = .01f * .01f;
	return glm::exp(-x * x / (2.f * sigma2)) / glm::sqrt(2.f * 3.141592f * sigma2);
}

inline
float Gaussian2D(float x, float y, float sigma)
{
	float sigma2 = sigma * sigma;
	float den = 1.f / glm::sqrt(sigma2 * 2.f * glm::pi<float>());
	float quad = ((x*x) + (y*y)) / (2.f * sigma2);
	return den * glm::exp(-quad);
}

// Arbitrary function to test flatland gaussian reconstruction
float GenerateFlatlandNormal(float x)
{
	return glm::cos(x*15.f) * .15f + .5f + glm::cos(x * 41.f) * .05f ;
}

float GenerateFlatlandNormalDerivative(float x)
{
	// Central differences, similar to what Sobel will do to a normal map
	float delta = .001f;
	float x1 = GenerateFlatlandNormal(x - delta);
	float x2 = GenerateFlatlandNormal(x + delta);
	return (x2 - x1) / (2.f * delta);
}

void GenerateFlatlandNormalImage(sf::Image& image, int width, int height)
{
	image.create(width, height, sf::Color::Black);

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			float fX = x / (float)width;

			float h = GenerateFlatlandNormal(fX) * height;
			float d = glm::abs(h - y) / 5.f;

			if (d < 1.f)
			{
				d = 1.f - d;
				image.setPixel(x, y, sf::Color(255 * d, 255 * d, 255 * d, 255));
			}
		}
	}
}

void CurvedElementsFlatlandNormalImage(sf::Image& image, int width, int height, float sigmaR)
{
	image.create(width, height, sf::Color::Black);

	// Amount of gaussians
	int m = 60;
	glm::vec3 * gaussianSeeds = new glm::vec3[m];
	glm::mat2 * invCovarianceMatrices = new glm::mat2[m];

	float h = 1.f / m;
	float sigmaH = h / glm::sqrt(8.f * glm::log(2.f));

	float sigmaH2 = sigmaH * sigmaH;
	float sigmaR2 = sigmaR * sigmaR;

	glm::mat2 invCov;
	for (int i = 0; i < m; i++)
	{
		float x = i * h; // ui
		float y = GenerateFlatlandNormal(x); // n(ui)
		float derivative = GenerateFlatlandNormalDerivative(x); // n'(ui)
		gaussianSeeds[i] = glm::vec3(x, y, derivative);

		invCov[0][0] = (1.f / sigmaH2) + (derivative * derivative / sigmaR2);
		invCov[1][0] = (-derivative / sigmaR2);
		invCov[1][1] = (1.f / sigmaR2);
		invCov[0][1] = (-derivative / sigmaR2);
		invCovarianceMatrices[i] = invCov;
	}

	sigmaH2 *= 2.f;
	sigmaR2 *= 2.f;

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			float accum = 0.f;

			for (int i = 0; i < m; i++)
			{
				glm::vec3 gaussianSeed = gaussianSeeds[i]; // Equivalent to (ui, n(ui))

				float fX = x / (float)width; // Equivalent to u
				float fY = y / (float)height; // Equivalent to s

				float Ci = 1.f;

				// Explicit approach
				//float positionBandGi = glm::pow(fX - gaussianSeed.x, 2.f) / sigmaH2;
				//float normalQuad = fY - gaussianSeed.y - (gaussianSeed.z * (fX - gaussianSeed.x));
				//float normalBandGi = glm::pow(-normalQuad, 2.f) / sigmaR2;
				//float Gi = ci * glm::exp(-positionBandGi) * glm::exp(-normalBandGi);

				// Matrix approach
				glm::vec2 dX(fX - gaussianSeed.x, fY - gaussianSeed.y);
				float quad = -0.5f * glm::dot(dX, invCovarianceMatrices[i] * dX);
				float Gi = Ci * glm::exp(quad);

				accum += Gi;
			}

			int c = glm::clamp(accum, 0.f, 1.f) * 255;
			image.setPixel(x, y, sf::Color(c, c, c, 255));
		}
	}

	delete[] gaussianSeeds;
}

void FlatElementsFlatlandNormalImage(sf::Image& image, int width, int height, float sigmaR)
{
	image.create(width, height, sf::Color::Black);

	// Amount of gaussians
	int m = 80;
	glm::vec2 * gaussianSeeds = new glm::vec2[m];

	float h = 1.f / m;
	float sigmaH = h / glm::sqrt(8.f * glm::log(2.f));

	for (int i = 0; i < m; i++)
	{
		float x = i * h;
		float y = GenerateFlatlandNormal(x);
		gaussianSeeds[i] = glm::vec2(x, y);
	}
	
	float sigmaH2 = 2.f * sigmaH * sigmaH;
	float sigmaR2 = 2.f * sigmaR * sigmaR;

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			float accum = 0.f;

			for (int i = 0; i < m; i++)
			{
				glm::vec2 gaussianSeed = gaussianSeeds[i]; // Equivalent to (ui, n(ui))

				float fX = x / (float)width; // Equivalent to u
				float fY = y / (float)height; // Equivalent to s

				float ci = 1.f;
				
				float positionBandGi = glm::pow(fX - gaussianSeed.x, 2.f) / sigmaH2;
				float normalBandGi = glm::pow(fY - gaussianSeed.y, 2.f) / sigmaR2;
				float Gi = ci * glm::exp(-positionBandGi) * glm::exp(-normalBandGi);

				accum += Gi;
			}

			int c = glm::clamp(accum, 0.f, 1.f) * 255;
			image.setPixel(x, y, sf::Color(c, c, c, 255));
		}
	}

	delete[] gaussianSeeds;
}

void BinningMethod(sf::Image& image, NaiveEstimator& estimator, int width, int height)
{
	image.create(width, height, sf::Color::Black);

	float * data = new float[width * height];
	estimator.BinningMethod(data, width, height);
	
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

sf::Image TestFlatland()
{
	int width = 512;
	int height = 512;


	// Intrinsic roughness
	float sigmaR = .005f;

	sf::Image image;
	GenerateFlatlandNormalImage(image, width, height);
	image.saveToFile("output/flatland_source.png");

	FlatElementsFlatlandNormalImage(image, width, height, sigmaR);
	image.saveToFile("output/flatland_flat_elements.png");

	CurvedElementsFlatlandNormalImage(image, width, height, sigmaR);
	image.saveToFile("output/flatland_curved_elements.png");
	
	return image;
}

int main()
{
	int width = 512;
	int height = 512;
	sf::RenderWindow window(sf::VideoMode(width, height), "NDF Estimator");
	sf::RectangleShape background(sf::Vector2f(width, height));

	std::string filename("cutlery.png");

	NaiveEstimator estimator("images/" + filename);

	//sf::Image image;
	sf::Image image = TestFlatland();
	//GenerateNDF(image, estimator, width, height);
	//BinningMethod(image, estimator, width, height);


	image.saveToFile("output/" + filename);

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