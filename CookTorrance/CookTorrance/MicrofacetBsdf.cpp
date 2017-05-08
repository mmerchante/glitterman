#include "MicrofacetBsdf.h"
#include <vector>
#include <string>

#include "glm\glm.hpp"
#include "glm\common.hpp"
#include "glm\matrix.hpp"
#include "glm\geometric.hpp"
#include "glm\trigonometric.hpp"
#include "glm\gtc\constants.hpp"

struct GaussianElementData
{
	// (u, n(u))
	glm::vec4 seed;

	// The matrices required for the 2d slice
	glm::mat2 A;
	glm::mat2 B;
	glm::mat2 C;

	// Coefficient that integrates to 1
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

	if (det > 0.f)
		return 1.f / glm::sqrt(det);

	return 0.f;
}

class Image
{
public:
	void Initialize(const char *filename)
	{
		unsigned char* redata = (unsigned char*)read_tga(filename, &width, &height);

		this->data = new RtColorRGB[width * height];

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int index = (y * width + x) * 3;

				RtColorRGB color;
				color.b = redata[index] / 255.0f;
				color.g = redata[index + 1] / 255.0f;
				color.r = redata[index + 2] / 255.0f;

				data[y * width + x] = color;
			}
		}
	}

	short le_short(unsigned char *bytes)
	{
		return bytes[0] | ((char)bytes[1] << 8);
	}

	void *read_tga(const char *filename, int *width, int *height)
	{
		struct tga_header {
			char  id_length;
			char  color_map_type;
			char  data_type_code;
			unsigned char  color_map_origin[2];
			unsigned char  color_map_length[2];
			char  color_map_depth;
			unsigned char  x_origin[2];
			unsigned char  y_origin[2];
			unsigned char  width[2];
			unsigned char  height[2];
			char  bits_per_pixel;
			char  image_descriptor;
		} header;
		int i, color_map_size, pixels_size;
		FILE *f;
		size_t read;
		void *pixels;
		/*errno_t err;
		err = */
		fopen_s(&f, filename, "rb");

		if (!f) {
			fprintf(stderr, "Unable to open %s for reading\n", filename);
			return NULL;
		}

		read = fread(&header, 1, sizeof(header), f);

		if (read != sizeof(header)) {
			fprintf(stderr, "%s has incomplete tga header\n", filename);
			fclose(f);
			return NULL;
		}
		if (header.data_type_code != 2) {
			fprintf(stderr, "%s is not an uncompressed RGB tga file\n", filename);
			fclose(f);
			return NULL;
		}

		for (i = 0; i < header.id_length; ++i)
			if (getc(f) == EOF) {
				fprintf(stderr, "%s has incomplete id string\n", filename);
				fclose(f);
				return NULL;
			}

		color_map_size = le_short(header.color_map_length) * (header.color_map_depth / 8);
		for (i = 0; i < color_map_size; ++i)
			if (getc(f) == EOF) {
				fprintf(stderr, "%s has incomplete color map\n", filename);
				fclose(f);
				return NULL;
			}

		*width = le_short(header.width); *height = le_short(header.height);
		pixels_size = *width * *height * (header.bits_per_pixel / 8);
		pixels = malloc(pixels_size);

		read = fread(pixels, 1, pixels_size, f);
		fclose(f);

		if (read != (unsigned int)pixels_size) {
			fprintf(stderr, "%s has incomplete image\n", filename);
			free(pixels);
			return NULL;
		}

		return pixels;
	}

	// Sobel
	glm::mat2 SampleNormalMapJacobian(glm::vec2 uv)
	{
		float hX = 1.f / (float)width;
		float hY = 1.f / (float)height;

		glm::vec2 xDiff = glm::vec2(0.f);
		xDiff += SampleNormalMap(uv + glm::vec2(hX, hY)) * 1.f;
		xDiff += SampleNormalMap(uv + glm::vec2(hX, 0)) * 2.f;
		xDiff += SampleNormalMap(uv + glm::vec2(hX, -hY)) * 1.f;

		xDiff += SampleNormalMap(uv - glm::vec2(hX, hY)) * -1.f;
		xDiff += SampleNormalMap(uv - glm::vec2(hX, 0)) * -2.f;
		xDiff += SampleNormalMap(uv - glm::vec2(hX, -hY)) * -1.f;

		glm::vec2 yDiff;
		yDiff += SampleNormalMap(uv + glm::vec2(hX, hY)) * 1.f;
		yDiff += SampleNormalMap(uv + glm::vec2(0, hY)) * 2.f;
		yDiff += SampleNormalMap(uv + glm::vec2(-hX, hY)) * 1.f;

		yDiff += SampleNormalMap(uv - glm::vec2(hX, hY)) * -1.f;
		yDiff += SampleNormalMap(uv - glm::vec2(0, hY)) * -2.f;
		yDiff += SampleNormalMap(uv - glm::vec2(-hX, hY)) * -1.f;

		xDiff /= 8.f * hX;
		yDiff /= 8.f * hY;

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

	// Bilinear filtering
	glm::vec2 SampleNormalMap(glm::vec2 uv)
	{
		uv = glm::fract(uv);

		int x = uv.x * width;
		int y = uv.y * height;

		glm::vec2 st = glm::vec2((uv.x * width) - x, (uv.y * height) - y);

		glm::vec2 p1 = SampleRawNormal(x, y);
		glm::vec2 p2 = SampleRawNormal(x + 1, y);
		glm::vec2 l1 = glm::mix(p1, p2, st.x);

		glm::vec2 p3 = SampleRawNormal(x, y + 1);
		glm::vec2 p4 = SampleRawNormal(x + 1, y + 1);
		glm::vec2 l2 = glm::mix(p3, p4, st.x);

		return glm::mix(l1, l2, st.y);
	}

	glm::vec2 SampleRawNormal(int x, int y)
	{
		RtColorRGB c = GetColor(x, y);
		return glm::vec2(c.r, c.g) * 2.0f - glm::vec2(1.0);
	}

	RtColorRGB GetColor(int x, int y)
	{
		x = glm::clamp(x, 0, width - 1);
		y = glm::clamp(y, 0, height - 1);
		return data[x + y * width];
	}

	int width;
	int height;
	RtColorRGB * data;
};

bool initialized = false;
Image normalMap;
GaussianElementData * gaussians = 0;

static const RtFloat k_minfacing = .0001f; // NdV < k_minfacing is invalid
static const unsigned char k_reflBlinnLobeId = 0;

static RixBXLobeSampled s_reflBlinnLobe;
static RixBXLobeTraits s_reflBlinnLobeTraits;

// BSDF Helper Inline Functions, from pbrt
inline RtFloat CosTheta(const RtVector3 &w) { return w.z; }
inline RtFloat Cos2Theta(const RtVector3 &w) { return w.z * w.z; }
inline RtFloat AbsCosTheta(const RtVector3 &w) { return std::abs(w.z); }
inline RtFloat Sin2Theta(const RtVector3 &w) { return std::max(0.f, 1.f - Cos2Theta(w)); }
inline RtFloat SinTheta(const RtVector3 &w) { return std::sqrt(Sin2Theta(w)); }
inline RtFloat TanTheta(const RtVector3 &w) { return SinTheta(w) / CosTheta(w); }
inline RtFloat Tan2Theta(const RtVector3 &w) { return Sin2Theta(w) / Cos2Theta(w); }

RtFloat RoughnessToAlpha(RtFloat roughness)
{
	// This seems to be closer to PRMan's current implementation...
	return std::max(roughness * roughness, .00001f);
}

MicrofacetBsdf::MicrofacetBsdf(RixShadingContext const * sc, RixBxdfFactory * bx, RixBXLobeTraits const & lobesWanted,
	RtColorRGB const * DiffuseColor, RtColorRGB const * SpecularColor, RtFloat const * Roughness, RtInt const * distributionUsed,
	RtColorRGB const * indexOfRefraction, RtColorRGB const * absorptionCoefficient) :
	RixBsdf(sc, bx),
	m_lobesWanted(lobesWanted),
	diffuseColor(DiffuseColor),
	specularColor(SpecularColor),
	microfacetRoughness(Roughness),
	distributionUsed(distributionUsed),
	indexOfRefraction(indexOfRefraction),
	extinctionCoefficient(absorptionCoefficient)
{
	RixBXLobeTraits lobes = s_reflBlinnLobeTraits;

	m_lobesWanted &= lobes;

	sc->GetBuiltinVar(RixShadingContext::k_P, &m_P);
	sc->GetBuiltinVar(RixShadingContext::k_Nn, &m_Nn);
	sc->GetBuiltinVar(RixShadingContext::k_Ngn, &m_Ngn);
	sc->GetBuiltinVar(RixShadingContext::k_Tn, &m_Tn);
	sc->GetBuiltinVar(RixShadingContext::k_Vn, &m_Vn);
	sc->GetBuiltinVar(RixShadingContext::k_u, &m_U);
	sc->GetBuiltinVar(RixShadingContext::k_v, &m_V);
	// TODO: get outside IOR
}

RixBXEvaluateDomain MicrofacetBsdf::GetEvaluateDomain()
{
	return k_RixBXReflect;
}

void MicrofacetBsdf::GetAggregateLobeTraits(RixBXLobeTraits * traits)
{
	*traits = m_lobesWanted;
}

// Samples in tangent space
RtVector3 MicrofacetBsdf::SampleGGXDistribution(const RtFloat2 & xi, const RtFloat alpha)
{
	RtFloat phi = xi.y * F_TWOPI;
	
	RtFloat tanTheta2 = alpha * alpha * xi.x / (1.0f - xi.x);
	RtFloat cosTheta = 1.f / std::sqrt(std::max(0.0001f, 1.f + tanTheta2));

	RtFloat sinTheta = std::sqrt(std::max(0.f, 1.f - cosTheta * cosTheta));
	return RixSphericalDirection(sinTheta, cosTheta, phi);
}

// Samples in tangent space
RtVector3 MicrofacetBsdf::SampleBeckmannDistribution(const RtFloat2& xi, const RtFloat alpha)
{
	RtFloat logSample = std::log(xi.x);

	if (std::isinf(logSample))
		logSample = 0.f;

	RtFloat tan2Theta = - (alpha * alpha * logSample);
	RtFloat phi = xi.y * F_TWOPI;

	RtFloat cosTheta = (1.f / std::sqrt(std::max(0.0001f, 1.f + tan2Theta)));
	RtFloat sinTheta = std::sqrt(std::max(0.f, 1.f - (cosTheta * cosTheta)));

	return RixSphericalDirection(sinTheta, cosTheta, phi);
}

RtVector3 MicrofacetBsdf::SampleDistribution(const RtFloat2 & xi, const RtFloat3 & normal, const RtFloat3 & tangent, 
											 const RtVector3& Wo, const RtFloat alpha, const RtInt distribution)
{
	RtFloat3 TY = Cross(normal, tangent);

	RtNormal3 localWh = distribution == 0 ? SampleBeckmannDistribution(xi, alpha) : SampleGGXDistribution(xi, alpha);
	RtVector3 Wh = RixChangeBasisFrom(localWh, tangent, TY, normal);

	if (Dot(Wh, Wo) < 0.f) 
		Wh *= -1.f;

	Wh.Normalize();
	RtVector3 Wi = RixReflect(Wo, Wh);
	Wi.Normalize();
	return Wi;
}

RtFloat MicrofacetBsdf::EvaluateMaskingShadow(const RtVector3 & Wi, const RtVector3 & Wo, const RtNormal3& normal, 
											   const RtNormal3& tangent, const RtFloat alpha, const RtInt& distribution)
{
	// TODO: Rethink lambdas so we dont need tangent space vectors
	RtFloat3 TY = Cross(normal, tangent);
	RtVector3 localWi = RixChangeBasisTo(Wi, tangent, TY, normal);
	RtVector3 localWo = RixChangeBasisTo(Wo, tangent, TY, normal);

	if(distribution == 0)
		return 1.f / (1.f + BeckmannLambda(localWi, alpha) + BeckmannLambda(localWo, alpha));
	else
		return 1.f / (1.f + GGXLambda(localWi, alpha) + GGXLambda(localWo, alpha));
}

RtFloat MicrofacetBsdf::BeckmannLambda(const RtVector3 & w, const RtFloat alpha)
{
	RtFloat absTanTheta = std::abs(TanTheta(w));
	
	if (std::isinf(absTanTheta))
		return 0.f;
	
	// Isotropic
	RtFloat a = 1.f / (alpha * absTanTheta);
	
	if (a >= 1.6f) 
		return 0.f;

	return (1.f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

RtFloat MicrofacetBsdf::GGXLambda(const RtVector3 & w, const RtFloat alpha)
{
	RtFloat absTanTheta = std::abs(TanTheta(w));

	if (std::isinf(absTanTheta))
		return 0.f;

	// Isotropic
	RtFloat alpha2Tan2Theta = alpha * absTanTheta;
	alpha2Tan2Theta *= alpha2Tan2Theta;

	return (-1.f + std::sqrt(std::max(0.0001f, 1.f + alpha2Tan2Theta))) * .5f;
}

void MicrofacetBsdf::GenerateSample(RixBXTransportTrait transportTrait, RixBXLobeTraits const * lobesWanted, RixRNG * rng, 
	RixBXLobeSampled * lobeSampled, RtVector3 * Ln, RixBXLobeWeights & W, RtFloat * FPdf, RtFloat * RPdf, RtColorRGB * compTrans)
{
	RtInt nPts = shadingCtx->numPts;
	RixBXLobeTraits all = GetAllLobeTraits();

	RtFloat2 * xi = (RtFloat2 *) RixAlloca(sizeof(RtFloat2) * nPts);
	rng->DrawSamples2D(nPts, xi);

	RtFloat * lobeR = (RtFloat *) RixAlloca(sizeof(RtFloat) * nPts);
	rng->DrawSamples1D(nPts, lobeR);

	RtColorRGB *reflDiffuseWgt = NULL;

	RtNormal3 facingNormal;
	RtFloat NdV;

	for (int i = 0; i < nPts; i++)
	{
		lobeSampled[i].SetValid(false);

		RixBXLobeTraits lobes = (all & lobesWanted[i]);
		bool doDiff = (lobes & s_reflBlinnLobeTraits).HasAny();

		if (!reflDiffuseWgt && doDiff)
			reflDiffuseWgt = W.AddActiveLobe(s_reflBlinnLobe);

		RtFloat const &uv_U = m_U[i];
		RtFloat const &uv_V = m_V[i];
		RtFloat2 uv = RtFloat2(uv_U, uv_V);

		if (doDiff)
		{
			facingNormal = RixGetForwardFacingNormal(m_Vn[i], m_Nn[i], &NdV);

			if (NdV > k_minfacing)
			{
				RtFloat alpha = RoughnessToAlpha(microfacetRoughness[i]);
				RtFloat cosTheta = 0.f;

				// TODO: Translate this to lobes...
				if(lobeR[i] > .5f || diffuseColor[i].IsBlack())
					Ln[i] = SampleDistribution(xi[i], m_Nn[i], m_Tn[i], m_Vn[i], alpha, distributionUsed[i]);
				else
					RixCosDirectionalDistribution(xi[i], m_Nn[i], Ln[i], cosTheta);

				RtFloat NdL = facingNormal.Dot(Ln[i]);

				if (NdL > 0.f)
				{
					EvaluateMicrofacetBRDF(uv, facingNormal, m_Tn[i], diffuseColor[i], specularColor[i], microfacetRoughness[i], distributionUsed[i],
						indexOfRefraction[i], extinctionCoefficient[i], Ln[i], m_Vn[i], reflDiffuseWgt[i], FPdf[i], RPdf[i]);

					lobeSampled[i] = s_reflBlinnLobe;
				}
			}
			// else invalid.. NullTrait
		}
	}
}

void MicrofacetBsdf::EvaluateSample(RixBXTransportTrait transportTrait, RixBXLobeTraits const * lobesWanted, RixRNG * rng, 
	RixBXLobeTraits * lobesEvaluated, RtVector3 const * Ln, RixBXLobeWeights & W, RtFloat * FPdf, RtFloat * RPdf)
{
	RtInt nPts = shadingCtx->numPts;
	RtNormal3 facingNormal;
	RtFloat NdV;

	RixBXLobeTraits all = GetAllLobeTraits();
	RtColorRGB * reflDiffuseWgt = NULL;

	for (int i = 0; i < nPts; i++)
	{
		lobesEvaluated[i].SetNone();
		RixBXLobeTraits lobes = (all & lobesWanted[i]);
		bool doDiff = (lobes & s_reflBlinnLobeTraits).HasAny();

		if (!reflDiffuseWgt && doDiff)
			reflDiffuseWgt = W.AddActiveLobe(s_reflBlinnLobe);

		RtFloat const &uv_U = m_U[i];
		RtFloat const &uv_V = m_V[i];
		RtFloat2 uv = RtFloat2(uv_U, uv_V);

		if (doDiff)
		{
			facingNormal = RixGetForwardFacingNormal(m_Vn[i], m_Nn[i], &NdV);

			if (NdV > k_minfacing)
			{
				RtFloat NdL = facingNormal.Dot(Ln[i]);

				if (NdL > 0.f)
				{
					EvaluateMicrofacetBRDF(uv, facingNormal, m_Tn[i], diffuseColor[i], specularColor[i], microfacetRoughness[i], distributionUsed[i],
						indexOfRefraction[i], extinctionCoefficient[i], Ln[i], m_Vn[i], reflDiffuseWgt[i], FPdf[i], RPdf[i]);

					lobesEvaluated[i] |= s_reflBlinnLobeTraits;
				}
			}
		}
	}
}

void MicrofacetBsdf::EvaluateSamplesAtIndex(RixBXTransportTrait transportTrait, RixBXLobeTraits const & lobesWanted, RixRNG * rng, 
	RtInt index, RtInt nsamps, RixBXLobeTraits * lobesEvaluated, RtVector3 const * Ln, RixBXLobeWeights & W, RtFloat * FPdf, RtFloat * RPdf)
{
	for (int i = 0; i < nsamps; i++)
		lobesEvaluated[i].SetNone();

	RixBXLobeTraits lobes = lobesWanted & GetAllLobeTraits();
	bool doDiff = (lobes & s_reflBlinnLobeTraits).HasAny();

	if (!doDiff)
		return;

	RtNormal3 const &Nn = m_Nn[index];
	RtNormal3 const &Ngn = m_Ngn[index];
	RtVector3 const &Vn = m_Vn[index];
	RtColorRGB const &diff = diffuseColor[index];
	RtColorRGB const &spec = specularColor[index];
	RtFloat const &roughness = microfacetRoughness[index];
	RtInt const &distribution = distributionUsed[index];
	RtFloat const &uv_U = m_U[index];
	RtFloat const &uv_V = m_V[index];

	RtFloat2 uv = RtFloat2(uv_U, uv_V);

	RtColorRGB const &ior = indexOfRefraction[index];
	RtColorRGB const &absorption = extinctionCoefficient[index];

	// Make any lobes that we may evaluate or write to active lobes,
	// initialize their lobe weights to zero and fetch a pointer to the
	// lobe weight arrays.
	RtColorRGB *reflDiffuseWgt = doDiff ? W.AddActiveLobe(s_reflBlinnLobe) : NULL;

	RtNormal3 Nf;
	RtFloat NdV;
	Nf = RixGetForwardFacingNormal(Vn, Nn, &NdV);
	
	if (NdV > k_minfacing)
	{
		for (int i = 0; i < nsamps; ++i)
		{
			RtFloat NdL = Nf.Dot(Ln[i]);

			if (NdL > 0.f)
			{
				EvaluateMicrofacetBRDF(uv, Nf, m_Tn[i], diff, spec, roughness, distribution, ior, absorption, Ln[i], Vn, reflDiffuseWgt[i], FPdf[i], RPdf[i]);
				lobesEvaluated[i] |= s_reflBlinnLobeTraits;
			}
		}
	}
}

float EvaluatePNDF(glm::vec2 uv, glm::vec2 st)
{
	int mX = normalMap.width;
	int mY = normalMap.height;

	glm::vec2 regionCenter = uv;
	glm::vec2 regionSize = glm::vec2(12, 12);
	glm::vec2 from = glm::clamp(regionCenter - regionSize * .5f, glm::vec2(), glm::vec2(mX - 1, mY - 1));
	glm::vec2 to = glm::clamp(regionCenter + regionSize * .5f, glm::vec2(), glm::vec2(mX - 1, mY - 1));

	float footprintRadius = regionSize.x * .5f / (float)mX;
	float sigmaP = footprintRadius * 0.0625f;

	glm::mat2 footprintCovarianceInv = glm::inverse(glm::mat2(sigmaP * sigmaP));
	glm::vec2 footprintMean = (from + regionSize * .5f) * glm::vec2(1.f / mX, 1.f / mY);

	// Outside our disk
	if (glm::length(st) > .975f)
		return 0.f;

	float accum = 0.f;

	// For each gaussian in the region...
	for (int gY = from.y; gY < to.y; gY++)
	{
		for (int gX = from.x; gX < to.x; gX++)
		{
			GaussianElementData data = gaussians[gY * mX + gX];
			glm::vec4 gaussianSeed = data.seed;

			// Direction, S - N(Xi)
			glm::vec2 S(st);
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

	accum /= (mX / (float)regionSize.x);

	return accum;
}

RtFloat MicrofacetBsdf::EvaluateBeckmann(RtFloat cosTheta, RtFloat roughness)
{
	RtFloat cos2Theta = cosTheta * cosTheta;
	RtFloat tan2Theta = std::max(0.f, 1.f - cos2Theta) / cos2Theta;

	if (std::isinf(tan2Theta))
		return 0.f;
	
	RtFloat cos4Theta = cos2Theta * cos2Theta;
	RtFloat alpha2 = RoughnessToAlpha(roughness);
	alpha2 *= alpha2;
	return std::exp(-tan2Theta / alpha2) / std::max(0.0001f, F_PI * alpha2 * cos4Theta);
}

RtFloat MicrofacetBsdf::EvaluateGGX(RtFloat cosTheta, RtFloat roughness)
{
	RtFloat cos2Theta = cosTheta * cosTheta;
	RtFloat tan2Theta = std::max(0.f, 1.f - cos2Theta) / cos2Theta;

	if (std::isinf(tan2Theta))
		return 0.f;

	RtFloat cos4Theta = cos2Theta * cos2Theta;
	RtFloat alpha2 = RoughnessToAlpha(roughness);
	alpha2 *= alpha2;

	RtFloat e = 1.f + (tan2Theta / alpha2);
	return F_INVPI / std::max(0.0001f, alpha2 * cos4Theta * e * e);
}

RtFloat MicrofacetBsdf::EvaluateDistribution(RtFloat cosTheta, RtFloat roughness, RtInt distribution)
{
	if (cosTheta > 0.f)
	{
		if (distribution == 0)
			return EvaluateBeckmann(cosTheta, roughness);
		else
			return EvaluateGGX(cosTheta, roughness);
	}

	return 0.f;
}

void MicrofacetBsdf::EvaluateMicrofacetBRDF(const RtFloat2& uv, const RtNormal3 & normal, const RtNormal3 & tangent, const RtColorRGB & diffuseColor,
	const RtColorRGB & specularColor, const RtFloat & roughness, const RtInt & distribution, const RtColorRGB & ior, const RtColorRGB & extinction, const RtVector3& Wi, const RtVector3& Wo, RtColorRGB & outRadiance,
	RtFloat & FPdf, RtFloat & RPdf)
{
	RtFloat cosThetaO = AbsDot(normal, Wo);
	RtFloat cosThetaI = AbsDot(normal, Wi);
	RtVector3 Wh = (Wi + Wo);

	if (cosThetaI == 0.f || cosThetaO == 0.f || (Wh.x == 0.f && Wh.y == 0.f && Wh.z == 0.f))
	{
		outRadiance = RtColorRGB(0.f);
		FPdf = 0.f;
		RPdf = 0.f;
		return;
	}

	// Evaluate our normal distribution function
	Wh.Normalize();
	RtFloat WhCosTheta = Dot(Wh, normal);
	//RtFloat D = EvaluateDistribution(WhCosTheta, roughness, distribution);

	RtFloat3 TY = Cross(normal, tangent);
	RtVector3 localWh = RixChangeBasisTo(Wh, tangent, TY, normal);

	RtFloat D = EvaluatePNDF(glm::vec2(uv.x, uv.y), glm::vec2(localWh.x, localWh.y));
	
	// Fresnel
	RtColorRGB F = RtColorRGB(0.f);
	RtColorRGB eta = ior.IsBlack() ? ior : (RtColorRGB(1.f) / ior);
	RtColorRGB kappa = extinction.IsBlack() ? extinction : (ior / extinction); // TODO: Find the actual mapping to kappa
	RixFresnelConductor(Dot(Wh, Wi), eta, kappa, &F);

	// Shadowing
	RtFloat G = EvaluateMaskingShadow(Wi, Wo, normal, tangent, RoughnessToAlpha(roughness), distribution);

	// The reference material does not conserve specular + diffuse energy
	RtColorRGB specularRadiance = specularColor * F * (D * G / (4.f * cosThetaO));
	RtColorRGB diffuseRadiance = diffuseColor * RixClamp(cosThetaI, 0.f, 1.f) / M_PI;
	outRadiance = specularRadiance + diffuseRadiance;

	// PDFs for bekcmann and ggx are the same with this approach
	RtFloat absCosThetaW = std::abs(WhCosTheta);
	FPdf = D * absCosThetaW / (4.f * Dot(Wo, Wh));
	RPdf = D * absCosThetaW / (4.f * Dot(Wi, Wh));

	// TODO: separate this into proper lobes...
	if (!diffuseRadiance.IsBlack())
	{
		FPdf = FPdf * .5f + F_INVTWOPI;
		RPdf = RPdf * .5f + F_INVTWOPI;
	}
}

extern "C" PRMANEXPORT RixBxdfFactory *CreateRixBxdfFactory(const char *hint)
{
	return new MicrofacetBsdfFactory();
}

extern "C" PRMANEXPORT void DestroyRixBxdfFactory(RixBxdfFactory *bxdf)
{
	delete (MicrofacetBsdfFactory *)bxdf;
}

/*-----------------------------------------------------------------------*/
MicrofacetBsdfFactory::MicrofacetBsdfFactory()
{
	m_DiffuseColorDflt = RtColorRGB(.5f);
	m_SpecularColorDefault = RtColorRGB(1.0f);
	m_RoughnessDefault = .35f;
	indexOfRefractionDefault = RtColorRGB(1.5f);
	extinctionCoefficientDefault = RtColorRGB(0.f); // Dielectric by default
	distributionUsedDefault = 0;
}

MicrofacetBsdfFactory::~MicrofacetBsdfFactory()
{
}


void InitializeGaussianData()
{
	if(gaussians != 0)
		delete [] gaussians;

	// Amount of gaussians
	int mX = normalMap.width; // Texel density
	int mY = normalMap.height; // Texel density
	gaussians = new GaussianElementData[mX * mY];

	float h = 1.f / mX; // Texel density
	float sigmaH = h / glm::sqrt(8.f * glm::log(2.f));
	float sigmaR = .0065f;

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

		glm::vec2 normal = normalMap.SampleNormalMap(Xi);
		data.seed = glm::vec4(Xi.x, Xi.y, normal.x, normal.y);

		glm::mat2 jacobian = normalMap.SampleNormalMapJacobian(Xi);
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
}

// Init
//  should be called once per RIB-instance. We look for parameter name
//  errors, and "cache" an understanding of our graph-evaluation requirements
//  in the form of allocation sizes.
int
MicrofacetBsdfFactory::Init(RixContext &ctx, char const *pluginpath)
{
	RixMessages* msgs = (RixMessages*)ctx.GetRixInterface(k_RixMessages);

	normalMap.Initialize("D:\\Penn\\cis660\\AuthoringTool\\glitterman\\ndf\\ndf\\ndf\\images\\snail.tga");
	msgs->Info("Normal map initialized.");

	InitializeGaussianData();
	msgs->Info("Gaussian data initialized.");

	return 0;
}

// Synchronize: delivers occasional status information
// from the renderer. Parameterlist contents depend upon the SyncMsg.
// This method is optional and the default implementation ignores all
// events.
void
MicrofacetBsdfFactory::Synchronize(RixContext &ctx, RixSCSyncMsg syncMsg,
	RixParameterList const *parameterList)
{
	if (syncMsg == k_RixSCRenderBegin)
	{
		s_reflBlinnLobe = RixBXLookupLobeByName(ctx, false, true, true, false,
			k_reflBlinnLobeId,
			"Specular");

		s_reflBlinnLobeTraits = RixBXLobeTraits(s_reflBlinnLobe);
	}
}

enum paramIds
{
	k_DiffuseColor,
	k_SpecularColor,
	k_Roughness,
	k_Distribution,
	k_IndexOfRefraction,
	k_ExtinctionCoefficient,
	k_numParams
};

RixSCParamInfo const *
MicrofacetBsdfFactory::GetParamTable()
{
	// see .args file for comments, etc...
	static RixSCParamInfo s_ptable[] =
	{
		RixSCParamInfo("DiffuseColor", k_RixSCColor),
		RixSCParamInfo("SpecularColor", k_RixSCColor),
		RixSCParamInfo("Roughness", k_RixSCFloat),
		RixSCParamInfo("Distribution", k_RixSCInteger),
		RixSCParamInfo("IndexOfRefraction", k_RixSCColor),
		RixSCParamInfo("ExtinctionCoefficient", k_RixSCColor),
		RixSCParamInfo() // end of table
	};
	return &s_ptable[0];
}

// CreateInstanceData:
//    analyze plist to determine our response to GetOpacityHints.
//    Checks these inputs:
//          transmissionBehavior (value),
//          presence (networked)
int
MicrofacetBsdfFactory::CreateInstanceData(RixContext &ctx,
	char const *handle,
	RixParameterList const *plist,
	InstanceData *idata)
{
	RtUInt64 req = k_TriviallyOpaque;
	idata->data = (void *)req; // no memory allocated, overload pointer
	idata->freefunc = NULL;
	return 0;
}

int MicrofacetBsdfFactory::GetInstanceHints(RtConstPointer instanceData) const
{
	InstanceHints const &hints = (InstanceHints const&)instanceData;
	return hints;
}

void MicrofacetBsdfFactory::Finalize(RixContext &)
{
}

RixBsdf * MicrofacetBsdfFactory::BeginScatter(RixShadingContext const *sCtx,
	RixBXLobeTraits const &lobesWanted,
	RixSCShadingMode sm,
	RtConstPointer instanceData)
{
	// Get all input data
	RtColorRGB const * diffuseColor;
	RtColorRGB const * specularColor;

	RtFloat const * roughness;

	RtInt const * distribution;

	RtColorRGB const * indexOfRefraction;
	RtColorRGB const * extinctionCoefficient;

	sCtx->EvalParam(k_DiffuseColor, -1, &diffuseColor, &m_DiffuseColorDflt, true);
	sCtx->EvalParam(k_SpecularColor, -1, &specularColor, &m_SpecularColorDefault, true);
	sCtx->EvalParam(k_Roughness, -1, &roughness, &m_RoughnessDefault, true);
	sCtx->EvalParam(k_Distribution, -1, &distribution, &distributionUsedDefault, true);
	sCtx->EvalParam(k_IndexOfRefraction, -1, &indexOfRefraction, &indexOfRefractionDefault, true);
	sCtx->EvalParam(k_ExtinctionCoefficient, -1, &extinctionCoefficient, &extinctionCoefficientDefault, true);

	RixShadingContext::Allocator pool(sCtx);
	void *mem = pool.AllocForBxdf<MicrofacetBsdf>(1);

	// Must use placement new to set up the vtable properly
	MicrofacetBsdf *eval = new (mem) MicrofacetBsdf(sCtx, this, lobesWanted, diffuseColor, specularColor, roughness, distribution,
		indexOfRefraction, extinctionCoefficient);

	return eval;
}

void
MicrofacetBsdfFactory::EndScatter(RixBsdf *)
{
}