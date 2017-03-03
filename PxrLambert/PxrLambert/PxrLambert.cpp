#pragma once
#include "PxrLambert.h"



	PxrLambert::PxrLambert(RixShadingContext const *sc, RixBxdfFactory *bx, RixBXLobeTraits const &lobesWanted, RtColorRGB const *DiffuseColor) :
		RixBsdf(sc, bx),
		m_lobesWanted(lobesWanted),
		m_DiffuseColor(DiffuseColor)
	{
		RixBXLobeTraits lobes = s_reflLambertLobeTraits;

		m_lobesWanted &= lobes;

		sc->GetBuiltinVar(RixShadingContext::k_Nn, &m_Nn);
		sc->GetBuiltinVar(RixShadingContext::k_Vn, &m_Vn);
	}

	RixBXEvaluateDomain PxrLambert::GetEvaluateDomain()
	{
		return k_RixBXReflect;   // Same as v19/20 kRixBXFront
	}

	void PxrLambert::GetAggregateLobeTraits(RixBXLobeTraits *t)
	{
		*t = m_lobesWanted;
	}

	void PxrLambert::GenerateSample(RixBXTransportTrait transportTrait,
		RixBXLobeTraits const *lobesWanted,
		RixRNG *rng,
		RixBXLobeSampled *lobeSampled,
		RtVector3   *Ln,
		RixBXLobeWeights &W,
		RtFloat *FPdf, RtFloat *RPdf,
		RtColorRGB* compTrans)
	{
		RtInt nPts = shadingCtx->numPts;
		RixBXLobeTraits all = GetAllLobeTraits();
		RtFloat2 *xi = (RtFloat2 *)RixAlloca(sizeof(RtFloat2) * nPts);
		rng->DrawSamples2D(nPts, xi);

		RtColorRGB *reflDiffuseWgt = NULL;

		RtNormal3 Nf;

		for (int i = 0; i < nPts; i++)
		{
			lobeSampled[i].SetValid(false);

			RixBXLobeTraits lobes = (all & lobesWanted[i]);
			bool doDiff = (lobes & s_reflLambertLobeTraits).HasAny();

			if (!reflDiffuseWgt && doDiff)
				reflDiffuseWgt = W.AddActiveLobe(s_reflLambertLobe);
			if (doDiff)
			{
				// we generate samples on the (front) side of Vn since
				// we have no translucence effects.
				RtFloat NdV;
				NdV = m_Nn[i].Dot(m_Vn[i]);
				if (NdV >= 0.f)
				{
					Nf = m_Nn[i];
				}
				else
				{
					Nf = -m_Nn[i];
					NdV = -NdV;
				}

				if (NdV > k_minfacing)
				{
					RtFloat NdL = Nf.Dot(Ln[i]);

					generate(NdL, m_DiffuseColor[i], reflDiffuseWgt[i]);

					FPdf[i] = 1.0f / Nf.Dot(m_Vn[i]);
					RPdf[i] = 1.0f / NdL;

					lobeSampled[i] = s_reflLambertLobe;
				}				
			}
		}

	}

	void PxrLambert::EvaluateSample(RixBXTransportTrait transportTrait,
		RixBXLobeTraits const *lobesWanted,
		RixRNG *rng,
		RixBXLobeTraits *lobesEvaluated,
		RtVector3 const *Ln, RixBXLobeWeights &W,
		RtFloat *FPdf, RtFloat *RPdf)

	{
		RtNormal3 Nf;
		RtInt nPts = shadingCtx->numPts;
		RixBXLobeTraits all = GetAllLobeTraits();

		RtColorRGB *reflDiffuseWgt = NULL;

		for (int i = 0; i < nPts; i++)
		{
			lobesEvaluated[i].SetNone();
			RixBXLobeTraits lobes = (all & lobesWanted[i]);
			bool doDiff = (lobes & s_reflLambertLobeTraits).HasAny();

			if (!reflDiffuseWgt && doDiff)
				reflDiffuseWgt = W.AddActiveLobe(s_reflLambertLobe);

			if (doDiff)
			{
				RtFloat NdV;
				NdV = m_Nn[i].Dot(m_Vn[i]);
				if (NdV >= 0.f)
					Nf = m_Nn[i];
				else
				{
					Nf = -m_Nn[i];
					NdV = -NdV;
				}
				if (NdV > k_minfacing)
				{
					RtFloat NdL = Nf.Dot(Ln[i]);
					if (NdL > 0.f)
					{
						evaluate(NdL, m_DiffuseColor[i], reflDiffuseWgt[i]);
						lobesEvaluated[i] |= s_reflLambertLobeTraits;

						FPdf[i] = 1.0f / Nf.Dot(m_Vn[i]);
						RPdf[i] = 1.0f / NdL;
					}
				}
			}
		}


	}

	void PxrLambert::EvaluateSamplesAtIndex(RixBXTransportTrait transportTrait,
		RixBXLobeTraits const &lobesWanted,
		RixRNG *rng,
		RtInt index, RtInt nsamps,
		RixBXLobeTraits *lobesEvaluated,
		RtVector3 const *Ln,
		RixBXLobeWeights &W,
		RtFloat *FPdf, RtFloat *RPdf)

	{
			for (int i = 0; i < nsamps; i++)
			lobesEvaluated[i].SetNone();

			RixBXLobeTraits lobes = lobesWanted & GetAllLobeTraits();
			bool doDiff = (lobes & s_reflLambertLobeTraits).HasAny();

			if (!doDiff)
				return;

			RtNormal3 const &Nn = m_Nn[index];
			RtVector3 const &Vn = m_Vn[index];
			RtColorRGB const &diff = m_DiffuseColor[index];

			// Make any lobes that we may evaluate or write to active lobes,
			// initialize their lobe weights to zero and fetch a pointer to the
			// lobe weight arrays.

			RtColorRGB *reflDiffuseWgt = doDiff
				? W.AddActiveLobe(s_reflLambertLobe) : NULL;

			RtNormal3 Nf;
			RtFloat NdV;

			NdV = Nn.Dot(Vn);
			if (NdV >= .0f)
				Nf = Nn;
			else
			{
				Nf = -Nn;
				NdV = -NdV;
			}
			RtFloat NfdV;
			NfdV = NdV;
			if (NdV > k_minfacing)
			{
				for (int i = 0; i < nsamps; ++i)
				{
					RtFloat NdL = Nf.Dot(Ln[i]);
					if (NdL > 0.f)
					{
						evaluate(NdL, diff, reflDiffuseWgt[i]);
						FPdf[i] = 1.0f / Nf.Dot(m_Vn[i]);
						RPdf[i] = 1.0f / NdL;
						lobesEvaluated[i] |= s_reflLambertLobeTraits;
					}
				}
			}		
	}



PxrLambertFactory::PxrLambertFactory()
{
	m_DiffuseColorDflt = RtColorRGB(1.0f, 0.0f, 0.0f);
}

PxrLambertFactory::~PxrLambertFactory()
{
}

// Init
//  should be called once per RIB-instance. We look for parameter name
//  errors, and "cache" an understanding of our graph-evaluation requirements
//  in the form of allocation sizes.
int
PxrLambertFactory::Init(RixContext &ctx, char const *pluginpath)
{
	return 0;
}

// Synchronize: delivers occasional status information
// from the renderer. Parameterlist contents depend upon the SyncMsg.
// This method is optional and the default implementation ignores all
// events.
void
PxrLambertFactory::Synchronize(RixContext &ctx, RixSCSyncMsg syncMsg,
	RixParameterList const *parameterList)
{
	if (syncMsg == k_RixSCRenderBegin)
	{
		s_reflLambertLobe = RixBXLookupLobeByName(ctx, false, true, true, false,
			k_reflLambertLobeId, "Diffuse");

		s_reflLambertLobeTraits = RixBXLobeTraits(s_reflLambertLobe);
	}
}

enum paramIds
{
	k_DiffuseColor,
	k_numParams
};

RixSCParamInfo const *
PxrLambertFactory::GetParamTable()
{
	// see .args file for comments, etc...
	static RixSCParamInfo s_ptable[] =
	{
		RixSCParamInfo("DiffuseColor", k_RixSCColor),
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
PxrLambertFactory::CreateInstanceData(RixContext &ctx,
	char const *handle,
	RixParameterList const *plist,
	InstanceData *idata)
{
	RtUInt64 req = k_TriviallyOpaque;
	idata->data = (void *)req; // no memory allocated, overload pointer
	idata->freefunc = NULL;
	return 0;
}

int
PxrLambertFactory::GetInstanceHints(RtConstPointer instanceData) const
{
	// our instance data is the RixBxdfFactory::InstanceHints bitfield.
	InstanceHints const &hints = (InstanceHints const&)instanceData;
	return hints;
}

// Finalize:
//  companion to Init, called with the expectation that any data
//  allocated there will be released here.
void
PxrLambertFactory::Finalize(RixContext &)
{
}
RixBsdf *
PxrLambertFactory::BeginScatter(RixShadingContext const *sCtx,
	RixBXLobeTraits const &lobesWanted,
	RixSCShadingMode sm,
	RtConstPointer instanceData)
{
	// Get all input data
	RtColorRGB const * DiffuseColor;
	sCtx->EvalParam(k_DiffuseColor, -1, &DiffuseColor, &m_DiffuseColorDflt, true);

	RixShadingContext::Allocator pool(sCtx);
	void *mem = pool.AllocForBxdf<PxrLambert>(1);

	// Must use placement new to set up the vtable properly
	PxrLambert *eval = new (mem) PxrLambert(sCtx, this, lobesWanted, DiffuseColor);

	return eval;
}

void
PxrLambertFactory::EndScatter(RixBsdf *)
{
}