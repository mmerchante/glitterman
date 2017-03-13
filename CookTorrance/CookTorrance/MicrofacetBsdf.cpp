#include "MicrofacetBsdf.h"

static const RtFloat k_minfacing = .0001f; // NdV < k_minfacing is invalid
static const unsigned char k_reflBlinnLobeId = 0;

static RixBXLobeSampled s_reflBlinnLobe;
static RixBXLobeTraits s_reflBlinnLobeTraits;

MicrofacetBsdf::MicrofacetBsdf(RixShadingContext const * sc, RixBxdfFactory * bx, RixBXLobeTraits const & lobesWanted,
	RtColorRGB const * DiffuseColor, RtColorRGB const * SpecularColor, RtFloat const * Roughness, RtFloat const * IOR) :
	RixBsdf(sc, bx),
	m_lobesWanted(lobesWanted),
	m_DiffuseColor(DiffuseColor),
	m_SpecularColor(SpecularColor),
	m_Roughness(Roughness),
	m_IOR(IOR)
{
	RixBXLobeTraits lobes = s_reflBlinnLobeTraits;

	m_lobesWanted &= lobes;

	sc->GetBuiltinVar(RixShadingContext::k_P, &m_P);
	sc->GetBuiltinVar(RixShadingContext::k_Nn, &m_Nn);
	sc->GetBuiltinVar(RixShadingContext::k_Ngn, &m_Ngn);
	sc->GetBuiltinVar(RixShadingContext::k_Tn, &m_Tn);
	sc->GetBuiltinVar(RixShadingContext::k_Vn, &m_Vn);
}

RixBXEvaluateDomain MicrofacetBsdf::GetEvaluateDomain()
{
	return k_RixBXReflect;
}

void MicrofacetBsdf::GetAggregateLobeTraits(RixBXLobeTraits * traits)
{
	*traits = m_lobesWanted;
}

void MicrofacetBsdf::GenerateSample(RixBXTransportTrait transportTrait, RixBXLobeTraits const * lobesWanted, RixRNG * rng, RixBXLobeSampled * lobeSampled, RtVector3 * Ln, RixBXLobeWeights & W, RtFloat * FPdf, RtFloat * RPdf, RtColorRGB * compTrans)
{
	RtInt nPts = shadingCtx->numPts;
	RixBXLobeTraits all = GetAllLobeTraits();

	RtFloat2 * xi = (RtFloat2 *) RixAlloca(sizeof(RtFloat2) * nPts);
	rng->DrawSamples2D(nPts, xi);

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

		if (doDiff)
		{
			facingNormal = -RixGetBackwardFacingNormal(m_Vn[i], m_Nn[i], &NdV);
			NdV *= -1.f;

			if (NdV > k_minfacing)
			{
				RtFloat cosTheta = 0.f;
				RixCosDirectionalDistribution(xi[i], facingNormal, Ln[i], cosTheta);

				RtFloat NdL = facingNormal.Dot(Ln[i]);

				if (NdL > 0.f)
				{
					EvaluateMicrofacetBRDF(NdV, NdL, facingNormal, m_DiffuseColor[i], m_SpecularColor[i], m_Roughness[i], m_IOR[i],
						Ln[i], m_Vn[i], reflDiffuseWgt[i], FPdf[i], RPdf[i]);

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

		if (doDiff)
		{
			facingNormal = -RixGetBackwardFacingNormal(m_Vn[i], m_Nn[i], &NdV);
			NdV *= -1.f;

			if (NdV > k_minfacing)
			{
				RtFloat NdL = facingNormal.Dot(Ln[i]);

				if (NdL > 0.f)
				{
					EvaluateMicrofacetBRDF(NdV, NdL, facingNormal, m_DiffuseColor[i], m_SpecularColor[i], m_Roughness[i], m_IOR[i],
						Ln[i], m_Vn[i], reflDiffuseWgt[i], FPdf[i], RPdf[i]);

					lobesEvaluated[i] |= s_reflBlinnLobeTraits;
				}
			}
		}
	}
}

void MicrofacetBsdf::EvaluateSamplesAtIndex(RixBXTransportTrait transportTrait, RixBXLobeTraits const & lobesWanted, RixRNG * rng, RtInt index, RtInt nsamps, RixBXLobeTraits * lobesEvaluated, RtVector3 const * Ln, RixBXLobeWeights & W, RtFloat * FPdf, RtFloat * RPdf)
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
	RtColorRGB const &diff = m_DiffuseColor[index];
	RtColorRGB const &spec = m_SpecularColor[index];
	RtFloat const &roughness = m_Roughness[index];
	RtFloat const &ior = m_IOR[index];

	// Make any lobes that we may evaluate or write to active lobes,
	// initialize their lobe weights to zero and fetch a pointer to the
	// lobe weight arrays.

	RtColorRGB *reflDiffuseWgt = doDiff ? W.AddActiveLobe(s_reflBlinnLobe) : NULL;

	RtNormal3 Nf;
	RtFloat NdV;
	Nf = -RixGetBackwardFacingNormal(Vn, Nn, &NdV);
	NdV *= -1.f;

	if (NdV > k_minfacing)
	{
		for (int i = 0; i < nsamps; ++i)
		{
			RtFloat NdL = Nf.Dot(Ln[i]);

			if (NdL > 0.f)
			{
				EvaluateMicrofacetBRDF(NdV, NdL, Nf, diff, spec, roughness, ior, Ln[i], Vn, reflDiffuseWgt[i], FPdf[i], RPdf[i]);
				lobesEvaluated[i] |= s_reflBlinnLobeTraits;
			}
		}
	}
}

// Reference: pbrtv3
RtFloat RoughnessToAlpha(RtFloat roughness) 
{
	roughness = std::max(roughness, .0001f);
	RtFloat x = std::log(roughness);
	return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
}

RtFloat MicrofacetBsdf::EvaluateBeckmann(RtFloat cosTheta, RtFloat roughness)
{
	RtFloat cos2Theta = cosTheta * cosTheta;
	RtFloat tan2Theta = std::max(0.f, 1.f - cos2Theta) / cos2Theta;

	if (std::isinf(tan2Theta))
		return 0.f;
	
	RtFloat cos4Theta = cos2Theta * cos2Theta;
	RtFloat alpha2 = (roughness);
	alpha2 *= alpha2;
	return std::exp(-tan2Theta / alpha2) / (M_PI * alpha2 * cos4Theta);
}

RtFloat MicrofacetBsdf::EvaluateDistribution(RtFloat cosTheta, RtFloat roughness)
{
	if (cosTheta > 0.f)
	{
		return EvaluateBeckmann(cosTheta, roughness);
	}

	return 0.f;
}

void MicrofacetBsdf::EvaluateMicrofacetBRDF(RtFloat NdV, RtFloat NdL, const RtNormal3 & normal, const RtColorRGB & diffuseColor, 
	const RtColorRGB & specularColor, const RtFloat & roughness, const RtFloat & ior, const RtVector3& Wi, const RtVector3& Wo, RtColorRGB & outRadiance,
	RtFloat & FPdf, RtFloat & RPdf)
{
	RtVector3 halfVector = (Wi + Wo);
	halfVector.Normalize();
	float halfVectorCosTheta = fabs(halfVector.Dot(normal));

	// Evaluate our normal distribution function
	float D = EvaluateDistribution(halfVectorCosTheta, .35f);

	float F = 0.f;
	RixFresnelDielectric(NdV, 1.f / ior, &F);

	// Shadow function
	float IdN = normal.Dot(Wo);
	float OdN = normal.Dot(Wi);

	float G1, G2;
	if (IdN <= 0.f)
		G1 = 0.f;
	else
	{
		float sinVSqrd = 1.f - IdN*IdN;
		float tanV = sqrtf(sinVSqrd / (IdN*IdN));
		float a = sqrtf(0.5f*roughness + 1.f) / tanV;
		if (a<1.6f)
		{
			G1 = (3.535f*a + 2.181f*a*a) / (1.f + 2.276f*a + 2.577f*a*a);
		}
		else
		{
			G1 = 1.f;
		}
	}
	if (OdN <= 0.f)
		G2 = 0.f;
	else
	{
		float sinVSqrd = 1.f - OdN*OdN;
		float tanV = sqrtf(sinVSqrd / (OdN*OdN));
		float a = sqrtf(0.5f*roughness + 1.f) / tanV;
		if (a<1.6f)
		{
			G2 = (3.535f*a + 2.181f*a*a) / (1.f + 2.276f*a + 2.577f*a*a);
		}
		else
		{
			G2 = 1.f;
		}
	}

	float radiance = F;// (G1 * G2 * D * F) / (4.f * IdN);
	outRadiance = specularColor * radiance;// +diffuseColor * NdL / M_PI;
	FPdf = D * G1 / IdN;
	RPdf = D * G2 / OdN;
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
	m_SpecualrColorDflt = RtColorRGB(1.0f);
	m_RoughnessDefault = .35f;
	m_IORDefault = 1.5f;
}

MicrofacetBsdfFactory::~MicrofacetBsdfFactory()
{
}

// Init
//  should be called once per RIB-instance. We look for parameter name
//  errors, and "cache" an understanding of our graph-evaluation requirements
//  in the form of allocation sizes.
int
MicrofacetBsdfFactory::Init(RixContext &ctx, char const *pluginpath)
{
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
	k_IOR,
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
		RixSCParamInfo("IOR", k_RixSCFloat),
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
	RtColorRGB const * DiffuseColor;
	RtColorRGB const * SpecularColor;
	RtFloat const * Roughness;
	RtFloat const * IOR;
	sCtx->EvalParam(k_DiffuseColor, -1, &DiffuseColor, &m_DiffuseColorDflt, true);
	sCtx->EvalParam(k_SpecularColor, -1, &SpecularColor, &m_SpecualrColorDflt, true);
	sCtx->EvalParam(k_Roughness, -1, &Roughness, &m_RoughnessDefault, true);
	sCtx->EvalParam(k_IOR, -1, &IOR, &m_IORDefault, true);

	RixShadingContext::Allocator pool(sCtx);
	void *mem = pool.AllocForBxdf<MicrofacetBsdf>(1);

	// Must use placement new to set up the vtable properly
	MicrofacetBsdf *eval = new (mem) MicrofacetBsdf(sCtx, this, lobesWanted, DiffuseColor, SpecularColor, Roughness, IOR);

	return eval;
}

void
MicrofacetBsdfFactory::EndScatter(RixBsdf *)
{
}