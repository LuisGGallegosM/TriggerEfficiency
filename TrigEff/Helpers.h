
#ifndef HELPERS
#define HELPERS

#include"TLorentzVector.h"

#include"../OniaIO/Data.h"

const float dRthreshold = 0.3f;
const float dPtThreshold = 0.3f;
const bool dPtThresholdEnabled = false;

struct HltobjEntry
{
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> mass;
};

using AccFunction= bool(*)(float,float);

bool isMatched(const TLorentzVector* recoMuon, const HltobjEntry* onMuons);
bool isMuonInAcceptance(float pt, float abseta);
bool isJPsiInAcceptance(float pt, float abseta);
bool isUpsilonInAcceptance(float pt, float abseta);
bool isPassQualityCuts(const OniaInput* in, int index);

#endif