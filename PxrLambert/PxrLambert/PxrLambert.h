#pragma once
#define RENDERMAN21
#define _USE_MATH_DEFINES

#include "RixBxdf.h"
#ifdef RENDERMAN21
#include "RixRNG.h"
#endif
#include "RixShadingUtils.h"
#include <cstring>

static const unsigned char k_reflLambertLobeId = 0;

static const RtFloat k_minfacing = .0001f; // NdV < k_minfacing is invalid

static RixBXLobeSampled s_reflLambertLobe;

static RixBXLobeTraits s_reflLambertLobeTraits;

class PxrLambert : public RixBsdf
{
	public:
	
		PxrLambert(RixShadingContext const *sc, RixBxdfFactory *bx, RixBXLobeTraits const &lobesWanted, RtColorRGB const *DiffuseColor);
		virtual RixBXEvaluateDomain GetEvaluateDomain();
		virtual void GetAggregateLobeTraits(RixBXLobeTraits *t);

		virtual void GenerateSample(RixBXTransportTrait transportTrait,
			RixBXLobeTraits const *lobesWanted,
			RixRNG *rng,
			RixBXLobeSampled *lobeSampled,
			RtVector3   *Ln,
			RixBXLobeWeights &W,
			RtFloat *FPdf, RtFloat *RPdf,
			RtColorRGB* compTrans);

		virtual void EvaluateSample(RixBXTransportTrait transportTrait,
			RixBXLobeTraits const *lobesWanted,
			RixRNG *rng,
			RixBXLobeTraits *lobesEvaluated,
			RtVector3 const *Ln, RixBXLobeWeights &W,
			RtFloat *FPdf, RtFloat *RPdf);

		virtual void EvaluateSamplesAtIndex(RixBXTransportTrait transportTrait,
			RixBXLobeTraits const &lobesWanted,
			RixRNG *rng,
			RtInt index, RtInt nsamps,
			RixBXLobeTraits *lobesEvaluated,
			RtVector3 const *Ln,
			RixBXLobeWeights &W,
			RtFloat *FPdf, RtFloat *RPdf);

	private:

		PRMAN_INLINE
		void generate(RtFloat NdL, const RtColorRGB &DiffuseColor, RtColorRGB  &W)
		{
			W = DiffuseColor * NdL;
		}

		PRMAN_INLINE
		void evaluate(RtFloat NdL, const RtColorRGB &DiffuseColor, RtColorRGB  &W)
		{
			W = DiffuseColor * NdL;
		}

	private:
		RixBXLobeTraits m_lobesWanted;
		RtColorRGB const *m_DiffuseColor;
		RtNormal3 const *m_Nn;
		RtVector3 const* m_Vn;
};

// PxrLambertFactory Implementation
class PxrLambertFactory : public RixBxdfFactory
{
public:


	PxrLambertFactory();
	~PxrLambertFactory();

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
	RtColorRGB m_DiffuseColorDflt;
};

extern "C" PRMANEXPORT RixBxdfFactory *CreateRixBxdfFactory(const char *hint)
{
	return new PxrLambertFactory();
}

extern "C" PRMANEXPORT void DestroyRixBxdfFactory(RixBxdfFactory *bxdf)
{
	delete (PxrLambertFactory *)bxdf;
}

