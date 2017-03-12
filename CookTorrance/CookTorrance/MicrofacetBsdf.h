#pragma once
#define RENDERMAN21
#define _USE_MATH_DEFINES

#include "RixBxdf.h"
#include "RixRNG.h"
#include "RixShadingUtils.h"
#include <cstring> // memset

class MicrofacetBsdf : public RixBsdf
{
public:
	MicrofacetBsdf(RixShadingContext const *sc, RixBxdfFactory *bx, RixBXLobeTraits const &lobesWanted,
		RtColorRGB const *DiffuseColor, RtColorRGB const *SpecularColor, RtFloat const *SpecularHardness);

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

	RtFloat EvaluateBeckmann(RtFloat cosTheta, RtFloat roughness);
	RtFloat EvaluateDistribution(RtFloat cosTheta, RtFloat roughness);

	PRMAN_INLINE
	void EvaluateMicrofacetBRDF(RtFloat NdV, RtFloat NdL, const RtNormal3 &normal, const RtColorRGB &diffuseColor,
			const RtColorRGB &specularColor, const RtFloat &roughness, const RtVector3 &Wi, const RtVector3 &Wo,
			RtColorRGB  &outRadiance, RtFloat &FPdf, RtFloat &RPdf);

private:
	RixBXLobeTraits m_lobesWanted;
	RtColorRGB const * m_DiffuseColor;
	RtColorRGB const * m_SpecularColor;
	RtFloat const * m_SpecularHardness;
	RtPoint3 const * m_P;
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
	RtColorRGB m_SpecualrColorDflt;

	//----------------------------------------------------------------------------------------------------------------------
	/// @brief Defualt exponent of our blinn BRDF
	//----------------------------------------------------------------------------------------------------------------------
	RtFloat m_SpecularHardnessDflt;
	//----------------------------------------------------------------------------------------------------------------------
};