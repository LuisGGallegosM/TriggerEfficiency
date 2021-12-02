
#include"Helpers.h"

bool isMatched(const TLorentzVector* recoMuon, const HltobjEntry* onMuons)
{
    int onlineMuonSize= onMuons->eta.size();
    float eta=recoMuon->Eta();
    float phi=recoMuon->Phi();
    float pt=recoMuon->Pt();
    for (int i=0;i<onlineMuonSize;i++)
    {
        float deltaEta= eta - onMuons->eta[i];
        float deltaPhi= phi - onMuons->phi[i];
        if (sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi) < dRthreshold)
        {
            if (!dPtThresholdEnabled) 
                return true;
            if (abs(( onMuons->pt[i] - pt )/pt) < dPtThreshold) 
                return true;
        }
    }
    return false;
}

bool isMuonInAcceptance(float pt, float abseta)
{
   if (abseta > 2.4f ) return false;
   if (abseta < 1.2f ) return pt >3.5f;
   if (abseta < 2.1f ) return pt >= 5.47f-1.89f*abseta;
   return pt >1.5f;
}

bool isJPsiInAcceptance(float pt, float abseta)
{
    if (abseta < 1.6f ) return pt > 6.5f;
    return pt > 1.5f;
}

bool isUpsilonInAcceptance(float pt, float abseta)
{
    if (abseta > 2.4f) return false;
    return pt > 3.5f;
}

bool isPassQualityCuts(const OniaInput* in, int index)
{
    if (in->reco_mu_nPixWMea[index] <=0 ) return false;
    if (in->reco_mu_nTrkWMea[index] <=5 ) return false;
    if (in->reco_mu_dxy[index] >=0.3 ) return false;
    return (in->reco_mu_dz[index] < 20.0);
}
