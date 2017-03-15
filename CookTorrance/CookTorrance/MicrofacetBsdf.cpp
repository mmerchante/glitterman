#include "MicrofacetBsdf.h"

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
	RtColorRGB const * DiffuseColor, RtColorRGB const * SpecularColor, RtFloat const * Roughness,
	RtColorRGB const * indexOfRefraction, RtColorRGB const * absorptionCoefficient) :
	RixBsdf(sc, bx),
	m_lobesWanted(lobesWanted),
	diffuseColor(DiffuseColor),
	specularColor(SpecularColor),
	microfacetRoughness(Roughness),
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

RtVector3 MicrofacetBsdf::SampleDistribution(const RtFloat2 & xi, const RtFloat3 & normal, const RtFloat3 & tangent, const RtVector3& Wo, const RtFloat alpha)
{
	RtFloat3 TY = Cross(normal, tangent);

	RtNormal3 localWh = SampleBeckmannDistribution(xi, alpha);
	RtVector3 Wh = RixChangeBasisFrom(localWh, tangent, TY, normal);

	if (Dot(Wh, Wo) < 0.f) 
		Wh *= -1.f;

	Wh.Normalize();
	RtVector3 Wi = RixReflect(Wo, Wh);
	Wi.Normalize();
	return Wi;
}

RtFloat MicrofacetBsdf::EvaluateMaskingShadow(const RtVector3 & Wi, const RtVector3 & Wo, const RtNormal3& normal, const RtNormal3& tangent, const RtFloat alpha)
{
	// TODO: Rethink lambdas so we dont need tangent space vectors
	RtFloat3 TY = Cross(normal, tangent);
	RtVector3 localWi = RixChangeBasisTo(Wi, tangent, TY, normal);
	RtVector3 localWo = RixChangeBasisTo(Wo, tangent, TY, normal);

	return 1.f / (1.f + BeckmannLambda(localWi, alpha) + BeckmannLambda(localWo, alpha));
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

		if (doDiff)
		{
			facingNormal = RixGetForwardFacingNormal(m_Vn[i], m_Nn[i], &NdV);

			if (NdV > k_minfacing)
			{
				RtFloat alpha = RoughnessToAlpha(microfacetRoughness[i]);
				RtFloat cosTheta = 0.f;

				// TODO: Translate this to lobes...
				if(lobeR[i] > .5f || diffuseColor[i].IsBlack())
					Ln[i] = SampleDistribution(xi[i], m_Nn[i], m_Tn[i], m_Vn[i], alpha);
				else
					RixCosDirectionalDistribution(xi[i], m_Nn[i], Ln[i], cosTheta);

				RtFloat NdL = facingNormal.Dot(Ln[i]);

				if (NdL > 0.f)
				{
					EvaluateMicrofacetBRDF(facingNormal, m_Tn[i], diffuseColor[i], specularColor[i], microfacetRoughness[i], indexOfRefraction[i], extinctionCoefficient[i],
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
			facingNormal = RixGetForwardFacingNormal(m_Vn[i], m_Nn[i], &NdV);

			if (NdV > k_minfacing)
			{
				RtFloat NdL = facingNormal.Dot(Ln[i]);

				if (NdL > 0.f)
				{
					EvaluateMicrofacetBRDF(facingNormal, m_Tn[i], diffuseColor[i], specularColor[i], microfacetRoughness[i], indexOfRefraction[i], extinctionCoefficient[i],
						Ln[i], m_Vn[i], reflDiffuseWgt[i], FPdf[i], RPdf[i]);

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
				EvaluateMicrofacetBRDF(Nf, m_Tn[i], diff, spec, roughness, ior, absorption, Ln[i], Vn, reflDiffuseWgt[i], FPdf[i], RPdf[i]);
				lobesEvaluated[i] |= s_reflBlinnLobeTraits;
			}
		}
	}
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

RtFloat MicrofacetBsdf::EvaluateDistribution(RtFloat cosTheta, RtFloat roughness)
{
	if (cosTheta > 0.f)
	{
		// TODO: Choose GGX or Beckmann
		return EvaluateBeckmann(cosTheta, roughness);
	}

	return 0.f;
}

void MicrofacetBsdf::EvaluateMicrofacetBRDF(const RtNormal3 & normal, const RtNormal3 & tangent, const RtColorRGB & diffuseColor,
	const RtColorRGB & specularColor, const RtFloat & roughness, const RtColorRGB & ior, const RtColorRGB & extinction, const RtVector3& Wi, const RtVector3& Wo, RtColorRGB & outRadiance,
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
	RtFloat D = EvaluateDistribution(WhCosTheta, roughness);
	
	// Fresnel
	RtColorRGB F = RtColorRGB(0.f);
	RtColorRGB eta = RtColorRGB(1.f) / ior;
	RtColorRGB kappa = ior / extinction; // TODO: Find the actual mapping to kappa
	RixFresnelConductor(Dot(Wh, Wi), eta, kappa, &F);

	// Shadowing
	RtFloat G = EvaluateMaskingShadow(Wi, Wo, normal, tangent, RoughnessToAlpha(roughness));

	// The reference material does not conserve specular + diffuse energy
	RtFloat specularRadiance = D * G / (4.f * cosThetaO);
	outRadiance = specularColor * F * specularRadiance + (diffuseColor * cosThetaI / M_PI);

	RtFloat absCosThetaW = std::abs(WhCosTheta);
	FPdf = ((D * absCosThetaW / (4.f * Dot(Wo, Wh))) + (F_INVPI)) * .5f;
	RPdf = ((D * absCosThetaW / (4.f * Dot(Wi, Wh))) + (F_INVPI)) * .5f;
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

	RtColorRGB const * indexOfRefraction;
	RtColorRGB const * extinctionCoefficient;

	sCtx->EvalParam(k_DiffuseColor, -1, &diffuseColor, &m_DiffuseColorDflt, true);
	sCtx->EvalParam(k_SpecularColor, -1, &specularColor, &m_SpecularColorDefault, true);
	sCtx->EvalParam(k_Roughness, -1, &roughness, &m_RoughnessDefault, true);
	sCtx->EvalParam(k_IndexOfRefraction, -1, &indexOfRefraction, &indexOfRefractionDefault, true);
	sCtx->EvalParam(k_ExtinctionCoefficient, -1, &extinctionCoefficient, &extinctionCoefficientDefault, true);

	RixShadingContext::Allocator pool(sCtx);
	void *mem = pool.AllocForBxdf<MicrofacetBsdf>(1);

	// Must use placement new to set up the vtable properly
	MicrofacetBsdf *eval = new (mem) MicrofacetBsdf(sCtx, this, lobesWanted, diffuseColor, specularColor, roughness, 
		indexOfRefraction, extinctionCoefficient);

	return eval;
}

void
MicrofacetBsdfFactory::EndScatter(RixBsdf *)
{
}