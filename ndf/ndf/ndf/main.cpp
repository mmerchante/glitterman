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
	// Central differences
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

inline
glm::vec2 SampleNormalMap(const sf::Image& normalMap, glm::vec2 uv)
{
	uv = glm::clamp(uv, glm::vec2(0.f), glm::vec2(1.f - glm::epsilon<float>()));
	sf::Color rawNormal = normalMap.getPixel(uv.x * normalMap.getSize().x, uv.y * normalMap.getSize().y);
	return glm::vec2(rawNormal.r / 255.f, rawNormal.g / 255.f) * 2.0f - glm::vec2(1.0);
}

inline
glm::mat2 SampleNormalMapJacobian(const sf::Image& normalMap, glm::vec2 uv, float h)
{
	glm::vec2 xDiff = SampleNormalMap(normalMap, uv + glm::vec2(h, 0.f)) - SampleNormalMap(normalMap, uv - glm::vec2(h, 0.f));
	glm::vec2 yDiff = SampleNormalMap(normalMap, uv + glm::vec2(0.f, h)) - SampleNormalMap(normalMap, uv - glm::vec2(0.f, h));

	// Encoded as 
	// fx/dx fx/dy
	// fy/dx fy/dy
	// Remember glm is column major
	return glm::mat2(xDiff.x, yDiff.x, xDiff.y, yDiff.y) * (.5f / h);
}

struct GaussianElementData
{
	glm::vec4 seed; // (u, n(u))

	// The matrices required for the 2d slice
	glm::mat2 A;
	glm::mat2 B;
	glm::mat2 C;
};

void CurvedElements4DNDF(sf::Image& image, sf::Image& normalMap, int width, int height, float sigmaR)
{
	image.create(width, height, sf::Color::Black);

	// Amount of gaussians
	int mX = normalMap.getSize().x; // Texel density
	int mY = normalMap.getSize().y; // Texel density
	GaussianElementData * gaussians = new GaussianElementData[mX * mY];

	float h = 1.f / mX; // Texel density
	float sigmaH = h / glm::sqrt(8.f * glm::log(2.f));

	float sigmaH2 = sigmaH * sigmaH;
	float sigmaR2 = sigmaR * sigmaR;

	float invSigmaH2 = 1.f / sigmaH2;
	float invSigmaR2 = 1.f / sigmaR2;

	glm::mat4 invCov;
	for (int i = 0; i < mX * mY; i++)
	{
		float x = (i % mX);
		float y = (i / mX);

		glm::vec2 Xi(x / (float)mX, y / (float)mY);

		GaussianElementData data;

		// If texel density is bigger, we need bilinear filtering! TODO
		glm::vec2 normal = SampleNormalMap(normalMap, Xi);
		data.seed = glm::vec4(Xi.x, Xi.y, normal.x, normal.y);

		glm::mat2 jacobian = SampleNormalMapJacobian(normalMap, Xi, h);
		glm::mat2 trJacobian = glm::transpose(jacobian);

		data.A = glm::mat2(invSigmaH2) + trJacobian * jacobian * invSigmaR2;
		data.B = -trJacobian * invSigmaR2;
		data.C = glm::mat2(invSigmaR2);

		gaussians[i] = data;
	}

	sigmaH2 *= 2.f;
	sigmaR2 *= 2.f;

	glm::vec2 regionSize = glm::vec2(32, 32);
	glm::vec2 from = glm::vec2(100, 100);
	glm::vec2 scalingFactor = glm::vec2((float)mX / regionSize.x, (float)mY / regionSize.y);


	float * ndf = new float[width * height];
	float cumulativeImage = 0.f;

	// For each direction, S
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			float accum = 0.f;
			float s = x / (float)width;
			float t = y / (float)height;

			// For each gaussian in the region...
			for (int gX = from.x; gX < regionSize.x + from.x; gX++)
			{
				for (int gY = from.y; gY < regionSize.y + from.y; gY++)
				{
					GaussianElementData data = gaussians[gY * mX + gX];
					glm::vec4 gaussianSeed = data.seed;

					glm::vec2 S((s * 2.f) - 1.f, (t * 2.f) - 1.f);

					//S = glm::vec2(gaussianSeed.z, gaussianSeed.w) - S;
					S -= glm::vec2(gaussianSeed.z, gaussianSeed.w);

					// We reduce the 4D gaussian into 2D by fixing S, see appendix
					glm::mat2 invCov = data.A;
					glm::vec2 u0 = -((glm::inverse(data.A)) * data.B) * S;

					float inner = glm::dot(S, data.C * S) - glm::dot(u0, data.A * u0);
					float c = 1.f * glm::exp(-0.5f * inner);

					// TODO: Calculate the resulting gaussian by multiplying Gp * Gi...

					// Analytically integrate the 2D gaussian
					// Remember, dimensions = 2
					glm::mat2 covariance = glm::inverse(invCov);
					float integral = c * 2.f * glm::pi<float>() * glm::sqrt(glm::determinant(covariance));

					accum += integral;
				}
			}

			cumulativeImage = glm::max(cumulativeImage, accum);			
			ndf[y * width + x] = accum;
		}
	}


	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			float n = ndf[y * width + x];
			int c = glm::clamp(n / cumulativeImage, 0.f, 1.f) * 255;
			image.setPixel(x, y, sf::Color(c, c, c, 255));
		}
	}

	delete[] ndf;
	delete[] gaussians;
}

void FlatElements4DNDF(sf::Image& image, sf::Image& normalMap, int width, int height, float sigmaR)
{
	image.create(width, height, sf::Color::Black);

	// Amount of gaussians
	int mX = width; // Texel density
	int mY = height; // Texel density
	glm::vec4 * gaussianSeeds = new glm::vec4[mX * mY];
	glm::mat4 * invCovarianceMatrices = new glm::mat4[mX * mY];

	float h = 1.f / mX; // Texel density
	float sigmaH = h / glm::sqrt(8.f * glm::log(2.f));

	float sigmaH2 = sigmaH * sigmaH;
	float sigmaR2 = sigmaR * sigmaR;
	
	glm::mat4 invCov;
	for (int i = 0; i < mX * mY; i++)
	{
		float x = (i % mX);
		float y = (i / mX);

		glm::vec2 Xi(x / (float)mX, y / (float)mY);

		// If texel density is bigger, we need bilinear filtering!
		glm::vec2 normal = SampleNormalMap(normalMap, Xi);
		gaussianSeeds[i] = glm::vec4(Xi.x + h * .5f, Xi.y + h * .5f, normal.x, normal.y);
		
		//gaussianSeeds[i] = glm::vec4(x, y, derivative);

		//invCov[0][0] = (1.f / sigmaH2) + (derivative * derivative / sigmaR2);
		//invCov[1][0] = (-derivative / sigmaR2);
		//invCov[1][1] = (1.f / sigmaR2);
		//invCov[0][1] = (-derivative / sigmaR2);
		//invCovarianceMatrices[i] = invCov;
	}

	sigmaH2 *= 2.f;
	sigmaR2 *= 2.f;

	glm::vec2 region = glm::vec2(45, 45);
	glm::vec2 scalingFactor = glm::vec2((float)mX / region.x, (float)mY / region.y);

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			float accum = 0.f;
			float fX = x / (float)width;
			float fY = y / (float)height;

			glm::vec2 Xi = glm::vec2(fX, fY) / scalingFactor;

			for (int jX = 0; jX < 30; jX++)
			{
				for (int jY = 0; jY < 30; jY++)
				{
					float s = (jX / 15.f);// *2.f - 1.f;
					float t = (jY / 15.f);// *2.f - 1.f;
					glm::vec2 S(s, t);

					for (int gX = 0; gX < region.x; gX++)
					{
						for (int gY = 0; gY < region.y; gY++)
						{
							glm::vec4 gaussianSeed = gaussianSeeds[gY * mX + gX]; // Equivalent to (ui, n(ui))

							glm::vec2 gaussianCenter = glm::vec2(gaussianSeed.x, gaussianSeed.y);
							glm::vec2 gaussianNormal(gaussianSeed.z, gaussianSeed.w);

							float positionDistance = glm::pow(glm::distance(Xi, gaussianCenter), 2.f);
							float normalDistance = glm::pow(glm::distance(S, gaussianNormal), 2.f);

							float Ci = 1.f;

							// Explicit approach
							float positionBandGi = positionDistance / sigmaH2;
							float normalBandGi = normalDistance / sigmaR2;
							float Gi = Ci * glm::exp(-positionBandGi) * glm::exp(-normalBandGi);

							accum += Gi;
						}
					}
				}
			}

			int c = glm::clamp(accum, 0.f, 1.f) * 255;
			image.setPixel(x, y, sf::Color(c, c, c, 255));
		}
	}

	delete[] gaussianSeeds;
	delete[] invCovarianceMatrices;
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
	int width = 128;
	int height = 128;


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


sf::Image Test4D(std::string normalMapFilename)
{
	int width = 256;
	int height = 256;

	// Intrinsic roughness
	float sigmaR = .005f;
	sf::Image image;

	sf::Image normalMapImage;

	if (normalMapImage.loadFromFile(normalMapFilename))
	{
		CurvedElements4DNDF(image, normalMapImage, width, height, sigmaR);
		image.saveToFile("output/4D_flat_elements.png");
	}

	return image;
}

int main()
{
	int width = 256;
	int height = 256;
	sf::RenderWindow window(sf::VideoMode(width, height), "NDF Estimator");
	sf::RectangleShape background(sf::Vector2f(width, height));

	std::string filename("blender.png");

	NaiveEstimator estimator("images/" + filename);

	//sf::Image image;
	sf::Image image = Test4D("images/" + filename);
	//sf::Image image = TestFlatland();
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