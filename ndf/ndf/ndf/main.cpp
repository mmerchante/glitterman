#include <SFML/Graphics.hpp>
#include <iostream>
#include <thread>

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

glm::vec2 SampleRawNormal(const sf::Image& normal, int x, int y)
{
	int w = normal.getSize().x;
	int h = normal.getSize().y;

	x = glm::clamp(x, 0, w - 1);
	y = glm::clamp(y, 0, h - 1);

	sf::Color c = normal.getPixel(x, y);
	return glm::vec2(c.r / 255.f, c.g / 255.f) * 2.0f - glm::vec2(1.0);
}

inline
glm::vec2 SampleNormalMap(const sf::Image& normalMap, glm::vec2 uv)
{
	uv = glm::fract(uv);
	int w = normalMap.getSize().x;
	int h = normalMap.getSize().y;
	
	int x = uv.x * w;
	int y = uv.y * h;

	glm::vec2 st = glm::vec2((uv.x * w) - x, (uv.y * h) - y);

	glm::vec2 p1 = SampleRawNormal(normalMap, x, y);
	glm::vec2 p2 = SampleRawNormal(normalMap, x + 1, y);
	glm::vec2 l1 = glm::mix(p1, p2, st.x);

	glm::vec2 p3 = SampleRawNormal(normalMap, x, y + 1);
	glm::vec2 p4 = SampleRawNormal(normalMap, x + 1, y + 1);
	glm::vec2 l2 = glm::mix(p3, p4, st.x);

	return glm::mix(l1, l2, st.y);
}

inline
glm::mat2 SampleNormalMapJacobian(const sf::Image& normalMap, glm::vec2 uv)
{
	float hX = 1.f / normalMap.getSize().x;
	float hY = 1.f / normalMap.getSize().y;

	glm::vec2 xDiff = glm::vec2(0.f);
	xDiff += SampleNormalMap(normalMap, uv + glm::vec2(hX, hY)) * 1.f;
	xDiff += SampleNormalMap(normalMap, uv + glm::vec2(hX, 0)) * 2.f;
	xDiff += SampleNormalMap(normalMap, uv + glm::vec2(hX, -hY)) * 1.f;

	xDiff += SampleNormalMap(normalMap, uv - glm::vec2(hX, hY)) * -1.f;
	xDiff += SampleNormalMap(normalMap, uv - glm::vec2(hX, 0)) * -2.f;
	xDiff += SampleNormalMap(normalMap, uv - glm::vec2(hX, -hY)) * -1.f;

	glm::vec2 yDiff;
	yDiff += SampleNormalMap(normalMap, uv + glm::vec2(hX, hY)) * 1.f;
	yDiff += SampleNormalMap(normalMap, uv + glm::vec2(0 , hY)) * 2.f;
	yDiff += SampleNormalMap(normalMap, uv + glm::vec2(-hX, hY)) * 1.f;

	yDiff += SampleNormalMap(normalMap, uv - glm::vec2(hX, hY)) * -1.f;
	yDiff += SampleNormalMap(normalMap, uv - glm::vec2(0 , hY)) * -2.f;
	yDiff += SampleNormalMap(normalMap, uv - glm::vec2(-hX, hY)) * -1.f;

	xDiff /= 8.f * hX;
	yDiff /= 8.f * hY;
	
	/*
	glm::vec2 xDiff = SampleNormalMap(normalMap, uv + glm::vec2(hX, 0.f)) - SampleNormalMap(normalMap, uv - glm::vec2(hX, 0.f));
	glm::vec2 yDiff = SampleNormalMap(normalMap, uv + glm::vec2(0.f, hY)) - SampleNormalMap(normalMap, uv - glm::vec2(0.f, hY));

	xDiff /= 2.f * hX;
	yDiff /= 2.f * hY;
	*/

	// Encoded as 
	// fx/dx fx/dy
	// fy/dx fy/dy
	// Remember glm is column major
	glm::mat2 J;
	J[0][0] = xDiff.x;
	J[0][1] = yDiff.x;

	J[1][0] = xDiff.y;
	J[1][1] = yDiff.y;
	return J;
}

struct GaussianElementData
{
	glm::vec4 seed; // (u, n(u))

	// The matrices required for the 2d slice
	glm::mat2 A;
	glm::mat2 B;
	glm::mat2 C;

	float coeff;
};

inline
float EvaluateGaussian(float c, const glm::vec2& x, const glm::vec2& u, const glm::mat2& InvCov)
{
	float inner = glm::dot(x - u, InvCov * (x - u));
	return c * glm::exp(-.5f * inner);
}

inline
float GetGaussianCoefficient(const glm::mat2& InvCov)
{
	float det = glm::determinant(2.f * glm::pi<float>() * glm::inverse(InvCov));

	if(det > 0.f)
		return 1.f / glm::sqrt(det);

	return 0.f;
}

void RenderJacobian(sf::Image& image, sf::Image& normalMap, int width, int height)
{
	width = normalMap.getSize().x;
	height = normalMap.getSize().y;
	image.create(width, height, sf::Color::Black);

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			glm::vec2 uv = glm::vec2(x / (float)width, y / (float)height);
			glm::mat2 J = SampleNormalMapJacobian(normalMap, uv);

			float det = glm::abs(glm::determinant(J));
			sf::Uint8 d = glm::clamp(det * 255.f, 0.f, 255.f);
			image.setPixel(x, y, sf::Color(d, d, d, 255));
		}
	}
}

void CurvedElementsIntegrationThread(int threadNumber, int bucketHeight, int width, int height, int mX, int mY, GaussianElementData * gaussians, float * ndf)
{
	glm::vec2 regionCenter = glm::vec2(256, 256);
	glm::vec2 regionSize = glm::vec2(256, 256);
	glm::vec2 from = regionCenter - regionSize * .5f;

	float footprintRadius = regionSize.x * .5f / (float)mX;
	float sigmaP = footprintRadius * .5f;

	glm::mat2 footprintCovarianceInv = glm::inverse(glm::mat2(sigmaP * sigmaP));
	glm::vec2 footprintMean = (from + regionSize * .5f) * glm::vec2(1.f / mX, 1.f / mY);

	pcg32 generator;
	generator.seed(14041956 + threadNumber * 127361);

	int samplesPerPixel = 8;

	int fromY = bucketHeight * threadNumber;
	int toY = glm::min(bucketHeight * (threadNumber + 1), height);

	float invW = 1.f / (float) width;
	float invH = 1.f / (float) height;

	// For each direction, S
	for (int y = fromY; y < toY; y++)
	{
		for (int x = 0; x < width; x++)
		{
			float accum = 0.f;
			float s = x * invW;
			float t = y * invH;

			glm::vec2 imageS = glm::vec2((s * 2.f) - 1.f, (t * 2.f) - 1.f);

			// Outside our disk
			if (glm::length(imageS) > .975f)
			{
				ndf[y * width + x] = 0.f;
				continue;
			}

			for (int sample = 0; sample < samplesPerPixel; sample++)
			{
				s = (x + generator.nextFloat()) * invW;
				t = (y + generator.nextFloat()) * invH;

				// For each gaussian in the region...
				for (int gX = from.x; gX < regionSize.x + from.x; gX++)
				{
					for (int gY = from.y; gY < regionSize.y + from.y; gY++)
					{
						GaussianElementData data = gaussians[gY * mX + gX];
						glm::vec4 gaussianSeed = data.seed;

						// Direction, S - N(Xi)
						glm::vec2 S((s * 2.f) - 1.f, (t * 2.f) - 1.f);
						S = S - glm::vec2(gaussianSeed.z, gaussianSeed.w);

						// We reduce the 4D gaussian into 2D by fixing S, see appendix
						glm::mat2 invCov = data.A;
						glm::vec2 u0 = -((glm::inverse(data.A)) * data.B) * S;
						float inner = glm::dot(S, data.C * S) - glm::dot(u0, data.A * u0);
						float c = data.coeff * glm::exp(-0.5f * inner);

						// Calculate the resulting gaussian by multiplying Gp * Gi
						glm::mat2 resultInvCovariance = invCov + footprintCovarianceInv;
						glm::mat2 resultCovariance = glm::inverse(resultInvCovariance);
						glm::vec2 resultMean = resultCovariance * (invCov * u0 + footprintCovarianceInv * (footprintMean - glm::vec2(gaussianSeed.x, gaussianSeed.y)));

						float resultC = EvaluateGaussian(c, resultMean, u0, invCov) *
							EvaluateGaussian(GetGaussianCoefficient(footprintCovarianceInv), resultMean, footprintMean - glm::vec2(gaussianSeed.x, gaussianSeed.y), footprintCovarianceInv);

						float det = (glm::determinant(resultCovariance * 2.f * glm::pi<float>()));

						if (det > 0.f)
							accum += resultC * glm::sqrt(det);
					}
				}
			}

			accum /= (mX / (float)regionSize.x) * .8f;
			accum /= samplesPerPixel;

			ndf[y * width + x] = accum;
		}
	}
}

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

	for (int i = 0; i < mX * mY; i++)
	{
		float x = (i % mX);
		float y = (i / mX);

		glm::vec2 Xi(x / (float)mX, y / (float)mY);

		GaussianElementData data;

		glm::vec2 normal = SampleNormalMap(normalMap, Xi);
		data.seed = glm::vec4(Xi.x, Xi.y, normal.x, normal.y);

		glm::mat2 jacobian = SampleNormalMapJacobian(normalMap, Xi);
		glm::mat2 trJacobian = glm::transpose(jacobian);

		data.A = ((trJacobian * jacobian) * invSigmaR2) + glm::mat2(invSigmaH2);
		data.B = -trJacobian * invSigmaR2;
		data.C = glm::mat2(invSigmaR2);

		glm::mat4 invCov4D;

		// Upper left
		invCov4D[0][0] = data.A[0][0];
		invCov4D[0][1] = data.A[0][1];
		invCov4D[1][0] = data.A[1][0];
		invCov4D[1][1] = data.A[1][1];

		// Upper right
		invCov4D[2][0] = data.B[0][0];
		invCov4D[2][1] = data.B[0][1];
		invCov4D[3][0] = data.B[1][0];
		invCov4D[3][1] = data.B[1][1];

		// Lower left
		glm::mat2 trB = -jacobian * invSigmaR2;
		invCov4D[0][2] = trB[0][0];
		invCov4D[0][3] = trB[0][1];
		invCov4D[1][2] = trB[1][0];
		invCov4D[1][3] = trB[1][1];

		// Lower right
		invCov4D[2][2] = data.C[0][0];
		invCov4D[2][3] = data.C[0][1];
		invCov4D[3][2] = data.C[1][0];
		invCov4D[3][3] = data.C[1][1];

		float det = glm::determinant(glm::inverse(invCov4D) * 2.f * glm::pi<float>());

		if (det <= 0.f)
			data.coeff = 0.f;
		else
			data.coeff = h * h / glm::sqrt(det);
		
		gaussians[i] = data;
	}

	int threadCount = 8;
	std::vector<std::thread> threads;
	int bucketHeight = height / threadCount;

	float * ndf = new float[width * height];

	for (int i = 0; i < threadCount; i++)
		threads.push_back(std::thread(CurvedElementsIntegrationThread, i, bucketHeight, width, height, mX, mY, gaussians, ndf));

	for (int i = 0; i < threadCount; i++)
		threads[i].join();

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			float n = ndf[y * width + x];
			glm::uint8 c = glm::clamp(n * 255.f, 0.f, 255.f);
			image.setPixel(x, y, sf::Color(c, c, c, 255));
		}
	}

	//pcg32 generator;
	//generator.seed(14041956);

	//// For each direction, S
	//for (int x = 0; x < width; x++)
	//{
	//	for (int y = 0; y < height; y++)
	//	{
	//		float accum = 0.f;

	//		for (int sample = 0; sample < samplesPerPixel; sample++)
	//		{
	//			float s = (x + generator.nextFloat()) / (float)width;
	//			float t = (y + generator.nextFloat()) / (float)height;

	//			glm::vec2 imageS = glm::vec2((s * 2.f) - 1.f, (t * 2.f) - 1.f);
	//			if (glm::length(imageS) > .975f)
	//			{
	//				image.setPixel(x, y, sf::Color(0, 0, 0, 255));
	//				continue;
	//			}

	//			// For each gaussian in the region...
	//			for (int gX = from.x; gX < regionSize.x + from.x; gX++)
	//			{
	//				for (int gY = from.y; gY < regionSize.y + from.y; gY++)
	//				{
	//					GaussianElementData data = gaussians[gY * mX + gX];
	//					glm::vec4 gaussianSeed = data.seed;

	//					// Direction, S - N(Xi)
	//					glm::vec2 S((s * 2.f) - 1.f, (t * 2.f) - 1.f);
	//					S = S - glm::vec2(gaussianSeed.z, gaussianSeed.w);

	//					// We reduce the 4D gaussian into 2D by fixing S, see appendix
	//					glm::mat2 invCov = data.A;
	//					glm::vec2 u0 = -((glm::inverse(data.A)) * data.B) * S;
	//					float inner = glm::dot(S, data.C * S) - glm::dot(u0, data.A * u0);
	//					float c = data.coeff * glm::exp(-0.5f * inner);

	//					// Calculate the resulting gaussian by multiplying Gp * Gi
	//					glm::mat2 resultInvCovariance = invCov + footprintCovarianceInv;
	//					glm::mat2 resultCovariance = glm::inverse(resultInvCovariance);
	//					glm::vec2 resultMean = resultCovariance * (invCov * u0 + footprintCovarianceInv * (footprintMean - glm::vec2(gaussianSeed.x, gaussianSeed.y)));

	//					float resultC = EvaluateGaussian(c, resultMean, u0, invCov) *
	//						EvaluateGaussian(GetGaussianCoefficient(footprintCovarianceInv), resultMean, footprintMean - glm::vec2(gaussianSeed.x, gaussianSeed.y), footprintCovarianceInv);

	//					float det = glm::determinant(resultCovariance * 2.f * glm::pi<float>());

	//					if (det > 0.f)
	//						accum += resultC * glm::sqrt(det);
	//				}
	//			}
	//		}

	//		accum /= (mX / (float)regionSize.x) * .85f;
	//		accum /= samplesPerPixel;

	//		int c = glm::clamp(accum * 255.f, 0.f, 255.f);
	//		image.setPixel(x, y, sf::Color(c, c, c, 255));
	//	}
	//}

	delete[] ndf;
	delete[] gaussians;
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
	int width = 256;
	int height = 256;

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
		image.saveToFile("output/4D_curved.png");

		//RenderJacobian(image, normalMapImage, width, height);
		//image.saveToFile("output/jacobian.png");
	}

	return image;
}

int main()
{
	int width = 512;
	int height = 512;
	sf::RenderWindow window(sf::VideoMode(width, height), "NDF Estimator");
	sf::RectangleShape background(sf::Vector2f(width, height));

	std::string filename("normal5.png");

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