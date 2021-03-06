#pragma once
#define RENDERMAN21
#define _USE_MATH_DEFINES

#include "RixBxdf.h"
#include "RixRNG.h"
#include "RixShadingUtils.h"
#include <cstring>

class MicrofacetBsdf : public RixBsdf
{
public:
	MicrofacetBsdf(RixShadingContext const *sc, RixBxdfFactory *bx, RixBXLobeTraits const &lobesWanted,
		RtColorRGB const *DiffuseColor, RtColorRGB const *SpecularColor, RtFloat const *Roughness, RtInt const * distributionUsed,
		RtColorRGB const * indexOfRefraction, RtColorRGB const * extinctionCoefficient);

	virtual RixBXEvaluateDomain GetEvaluateDomain();

	virtual void GetAggregateLobeTraits(RixBXLobeTraits * traits);

	virtual void GenerateSample(RixBXTransportTrait transportTrait, RixBXLobeTraits const *lobesWanted, RixRNG *rng, RixBXLobeSampled *lobeSampled,
		RtVector3 * Ln, RixBXLobeWeights &W, RtFloat *FPdf, RtFloat *RPdf, RtColorRGB* compTrans);

	virtual void EvaluateSample(RixBXTransportTrait transportTrait, RixBXLobeTraits const *lobesWanted, RixRNG *rng, 
		RixBXLobeTraits *lobesEvaluated, RtVector3 const *Ln, RixBXLobeWeights &W, RtFloat *FPdf, RtFloat *RPdf);

	virtual void EvaluateSamplesAtIndex(RixBXTransportTrait transportTrait, RixBXLobeTraits const &lobesWanted,
		RixRNG *rng, RtInt index, RtInt nsamps, RixBXLobeTraits *lobesEvaluated, RtVector3 const *Ln,
		RixBXLobeWeights &W, RtFloat *FPdf, RtFloat *RPdf);

private:

	RtVector3 SampleGGXDistribution(const RtFloat2& xi, const RtFloat alpha);
	RtVector3 SampleBeckmannDistribution(const RtFloat2& xi, const RtFloat alpha);
	RtVector3 SampleDistribution(const RtFloat2& xi, const RtFloat3& normal, const RtFloat3& tangent, 
								 const RtVector3& Wo, const RtFloat alpha, const RtInt distribution);

	RtFloat EvaluateMaskingShadow(const RtVector3& Wi, const RtVector3& Wo, const RtNormal3& normal, const RtNormal3& tangent, 
								  const RtFloat alpha, const RtInt& distribution);
	RtFloat BeckmannLambda(const RtVector3& W, const RtFloat alpha);
	RtFloat GGXLambda(const RtVector3& W, const RtFloat alpha);

	RtFloat EvaluateBeckmann(RtFloat cosTheta, RtFloat roughness);
	RtFloat EvaluateGGX(RtFloat cosTheta, RtFloat roughness);
	RtFloat EvaluateDistribution(RtFloat cosTheta, RtFloat roughness, RtInt distribution);

	PRMAN_INLINE
	void EvaluateMicrofacetBRDF(const RtFloat2& uv, const RtNormal3 & normal, const RtNormal3 & tangent, const RtColorRGB & diffuseColor,
		const RtColorRGB & specularColor, const RtFloat & roughness, const RtInt & distribution, const RtColorRGB & ior, const RtColorRGB & absorption,
		const RtVector3& Wi, const RtVector3& Wo, RtColorRGB & outRadiance, RtFloat & FPdf, RtFloat & RPdf);

private:
	RixBXLobeTraits m_lobesWanted;
	
	RtColorRGB const * diffuseColor;
	RtColorRGB const * specularColor;
	
	RtFloat const * microfacetRoughness;

	// 0 Beckmann, 1 GGX
	RtInt const * distributionUsed;
	
	RtColorRGB const * indexOfRefraction;
	RtColorRGB const * extinctionCoefficient;

	RtPoint3 const * m_P;
	RtFloat const * m_U;
	RtFloat const * m_V;
	RtVector3 const * m_Vn;
	RtVector3 const * m_Tn;
	RtNormal3 const * m_Nn;
	RtNormal3 const * m_Ngn;
};

// PxrBlinnFactory Implementation
class MicrofacetBsdfFactory : public RixBxdfFactory
{
public:

	MicrofacetBsdfFactory();
	~MicrofacetBsdfFactory();

	virtual int Init(RixContext &, char const *pluginpath);
	RixSCParamInfo const *GetParamTable();
	virtual void Finalize(RixContext &);

	virtual void Synchronize(RixContext &ctx, RixSCSyncMsg syncMsg,
		RixParameterList const *parameterList);

	virtual int CreateInstanceData(RixContext &,
		char const *handle,
		RixParameterList const *,
		InstanceData *id);

	virtual int GetInstanceHints(RtConstPointer instanceData) const;

	virtual RixBsdf *BeginScatter(RixShadingContext const *,
		RixBXLobeTraits const &lobesWanted,
		RixSCShadingMode sm,
		RtConstPointer instanceData);
	virtual void EndScatter(RixBsdf *);

private:
	// these hold the default (def) values
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief Defualt colour of our blinn BRDF
	//----------------------------------------------------------------------------------------------------------------------
	RtColorRGB m_DiffuseColorDflt;
	RtColorRGB m_SpecularColorDefault;

	RtFloat m_RoughnessDefault;

	RtInt distributionUsedDefault;

	RtColorRGB indexOfRefractionDefault;
	RtColorRGB extinctionCoefficientDefault;
};