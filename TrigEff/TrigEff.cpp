
#include"../OniaIO/OniaIO.h"
#include"../Utils/Utils.h"
#include"TrigEff.h"
#include"Helpers.h"

#include<iostream>
#include<unordered_map>

#include"TFile.h"
#include"TEfficiency.h"

const float ptMax=100.0f;

const bool muProcess=true;
const bool dimuProcess=true;

using std::cout;
using std::cerr;
using HltIndex=std::unordered_map<Long64_t,HltobjEntry>;

void ProcessMuon(Input* input, Output* output,const HltIndex& indexer);
void ProcessDimuon(Input* input, Output* output,const HltIndex& indexer);
HltIndex generateIndexer(Input* input);

Output allocateOutput();

void TrigEff(const char* oniaFilename, const char* triggerFilename, const char* triggerName, const char* outputFilename, const char* oniatype)
{
    //init inputs
    Input input;

    //check if is L1 using the trigger name
    input.isL1=std::string(triggerName).find("HLT_HIL1")!= std::string::npos;
    std::string oniat=std::string(oniatype);
    if (input.isL1) std::cout <<"L1 trigger detected\n";

    if (dimuProcess)
    {
        if (oniat == "Y" ) 
        {
            input.type = OniaType::Upsilon;
            cout << "Using Upsilon dimuons\n";
        }
        else if (oniat == "JPsi" )
        {
            input.type = OniaType::JPsi;
            cout << "Using JPsi dimuons\n";
        } 
        else 
        { 
            cout << "Error, onia type '" << oniat << "' not recognized.\n";
            cout << "Available types: Y , JPsi\n";
            return; 
            }
    }

    //open input files
    TFile* oniaFile = OpenFile(oniaFilename);
    if (oniaFile==nullptr) return;
    input.oniaTree =OpenTree(oniaFile,oniaTreeName);
    if(input.oniaTree ==nullptr) return;
    
    TFile* triggerFile = OpenFile(triggerFilename);
    if (triggerFile==nullptr) return;

    std::string triggerPath=std::string(hltobjDirectoryName)+triggerName;
    input.hltobjectTree=OpenTree(triggerFile,triggerPath.data());
    if (input.hltobjectTree==nullptr) return;
    input.hltanalysisTree=OpenTree(triggerFile,hltanalysisTreeName);
    if(input.hltanalysisTree==nullptr) return;

    std::string outFilename=outputFilename;
    TFile* outputFile = CreateFile(outFilename+"/output.root");
    if(outputFile==nullptr) return;

    //init outputs
    Output output=allocateOutput();

    const HltIndex indexer= generateIndexer(&input);

    if(muProcess)
    {
        ProcessMuon(&input,&output,indexer);
        output.muonpass->Write();
        output.muontotal->Write();
    }
    
    if(dimuProcess)
    {
        ProcessDimuon(&input,&output,indexer);
        output.dimuonpass->Write();
        output.dimuontotal->Write();
    }
    
    cout << "Success.\n";

    delete oniaFile;
    delete triggerFile;
    delete outputFile;
}

HltIndex generateIndexer(Input* input)
{
    HltIndex indexer;
    Reader<HltobjInput> hltObj(input->hltobjectTree);
    Reader<HltanalysisInput> hltAn(input->hltanalysisTree);

    Long64_t entryNum= input->hltanalysisTree->GetEntries();

    indexer.reserve(entryNum);

    std::cout << "indexing events..\n";
    long long repeatedEvents=0;

    for(Long64_t entry=0;entry < entryNum;entry++)
    {
        const HltanalysisInput* analysisInput = hltAn.readEntry(entry);
        const HltobjInput* hltobjInput = hltObj.readEntry(entry);

        auto item= indexer.emplace(analysisInput->event,HltobjEntry());
        
        if(!item.second)
        {
            repeatedEvents++;
        }

        HltobjEntry* hltentry= &(item.first->second);
        hltentry->eta=*(hltobjInput->eta);
        hltentry->mass=*(hltobjInput->mass);
        hltentry->phi=*(hltobjInput->phi);
        hltentry->pt=*(hltobjInput->pt);
    }

    std::cout << "Indexing finished, index generated :" << entryNum-repeatedEvents << "\n";
    std::cout << "repeated events: " << repeatedEvents << "\n"; 
    return indexer;

}

void ProcessMuon(Input* input, Output* output,const HltIndex& indexer)
{
    Reader<OniaInput> onia(input->oniaTree);

    Writer<OniaOutput> total(output->muontotal);
    Writer<OniaOutput> pass(output->muonpass);
    
    Long64_t oniaEntryNum= input->oniaTree->GetEntries();
    Long64_t entryStep= oniaEntryNum/50;

    cout << "Proccessing muons\n";

    for(Long64_t entry=0;entry< oniaEntryNum;entry++)
    {
        if ((entry % entryStep)==0)
        {
            cout << "processing entries : " << round((100.0f*entry)/oniaEntryNum)
                 << "% " << entry << " / " << oniaEntryNum << '\n';
        }
            
        const OniaInput* oniaInput = onia.readEntry(entry);

	    int size = oniaInput->reco_mu_size;
        for(int iMu=0;iMu< size;iMu++)
        {
            const TLorentzVector* mu= input->isL1 ? 
                    (TLorentzVector*) oniaInput->reco_mu_L1_mom4->At(iMu) 
                :   (TLorentzVector*) oniaInput->reco_mu_mom4->At(iMu);
            
            const float pt= mu->Pt();
            const float y= mu->Rapidity();
            const float eta= mu->Eta();
            const int cent= oniaInput->centrality;

            if (pt>ptMax) continue;
            if (!isMuonInAcceptance(pt,abs(eta))) continue;
            if (!isPassQualityCuts(oniaInput,iMu)) continue;

            //Passed acceptance and quality cuts
            total.output.pt=pt;
            total.output.y=y;
            total.output.cent = cent;
            total.output.m= mu->M();
            total.output.eta=eta;
            total.writeEntry();

            //read hltobj, if not found, continue
            const auto hltobjFound= indexer.find(oniaInput->event);
            if (hltobjFound== indexer.end()) continue;
            if (isMatched(mu,&(hltobjFound->second)))
            {
                pass.output= total.output;
                pass.writeEntry();
            }
        }
    }

    cout << "finished processing "<< oniaEntryNum << " entries\n";
}

void ProcessDimuon(Input* input, Output* output,const HltIndex& indexer)
{
    Reader<OniaInput> onia(input->oniaTree);

    Writer<OniaOutput> total(output->dimuontotal);
    Writer<OniaOutput> pass(output->dimuonpass);
    
    Long64_t oniaEntryNum= input->oniaTree->GetEntries();
    Long64_t entryStep= oniaEntryNum/50;

    cout << "Proccessing dimuons\n";

    AccFunction isInAcceptance = nullptr;
    switch(input->type)
    {
        case OniaType::JPsi:
        isInAcceptance = isJPsiInAcceptance;
        break;
        case OniaType::Upsilon:
        isInAcceptance = isUpsilonInAcceptance;
        break;
    }

    for(Long64_t entry=0;entry< oniaEntryNum;entry++)
    {
        if ((entry % entryStep)==0)
        {
            cout << "processing entries : " << round((100.0f*entry)/oniaEntryNum)
                 << "% " << entry << " / " << oniaEntryNum << '\n';
        }
            
        const OniaInput* oniaInput = onia.readEntry(entry);

	    int size = oniaInput->reco_QQ_size;
        for(int iQQ=0;iQQ< size;iQQ++)
        {
            const int idx_pl=oniaInput->reco_QQ_mupl_idx[iQQ];
            const int idx_mi=oniaInput->reco_QQ_mumi_idx[iQQ];

            if( (idx_pl<0) || (idx_mi<0) ) continue;

            const TLorentzVector* dimu  =    (TLorentzVector*) oniaInput->reco_QQ_mom4->At(iQQ);

            const TLorentzVector* simu_pl= input->isL1 ? 
                    (TLorentzVector*) oniaInput->reco_mu_L1_mom4->At(idx_pl) 
                :   (TLorentzVector*) oniaInput->reco_mu_mom4->At(idx_pl);
            
            const TLorentzVector* simu_mi= input->isL1 ? 
                    (TLorentzVector*) oniaInput->reco_mu_L1_mom4->At(idx_mi) 
                :   (TLorentzVector*) oniaInput->reco_mu_mom4->At(idx_mi);
            
            const float pt= dimu->Pt();
            const float y= dimu->Rapidity();
            const float eta= dimu->Eta();
            const int cent= oniaInput->centrality;

            const float pt_pl=simu_pl->Pt();
            const float abseta_pl=fabs(simu_pl->Eta());
            const float pt_mi=simu_mi->Pt();
            const float abseta_mi=fabs(simu_mi->Pt());

            if (!(isInAcceptance(pt_pl,abseta_pl) && isInAcceptance(pt_mi,abseta_mi))) continue;
            if (!(isPassQualityCuts(oniaInput,idx_mi) && isPassQualityCuts(oniaInput,idx_pl))) continue;

            //Passed acceptance and quality cuts
            total.output.pt=pt;
            total.output.y=y;
            total.output.cent = cent;
            total.output.m= dimu->M();
            total.output.eta=eta;
            total.writeEntry();

            //read hltobj, if not found, continue
            const auto hltobjFound= indexer.find(oniaInput->event);
            if (hltobjFound== indexer.end()) continue;

            if (isMatched(simu_mi,&(hltobjFound->second)) && isMatched(simu_pl,&(hltobjFound->second)) )
            {
                pass.output= total.output;
                pass.writeEntry();
            }
        }
    }

    cout << "finished processing "<< oniaEntryNum << " entries\n";
}

Output allocateOutput()
{
    Output out{nullptr,nullptr,nullptr,nullptr};

    if(muProcess)
    {
        out.muonpass = new TTree("muon_pass","passing muons");
        out.muontotal = new TTree("muon_total","total muons");
    }

    if(dimuProcess)
    {
        out.dimuonpass = new TTree("dimuon_pass","passing dimuons");
        out.dimuontotal = new TTree("dimuon_total","total dimuons");
    }

    return out;
}

#if !defined(__CLING__)

int main(int argc, char **argv)
{
    if (argc==6)
        TrigEff(argv[1],argv[2],argv[3],argv[4],argv[5]);
    else
        cerr << "Wrong number of parameters:\nUsage:\nOnia Filename , Trigger Filename, Trigger Name, OutputFilename, OniaType\n";
    return 0;
}

#endif
