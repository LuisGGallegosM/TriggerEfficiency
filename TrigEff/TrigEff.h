
#include"RtypesCore.h"

#include<array>
#include"TH1.h"
#include"TTree.h"

struct Output
{
    TTree* muontotal;
    TTree* muonpass;
    TTree* dimuontotal;
    TTree* dimuonpass;
};

enum class OniaType { Upsilon, JPsi };

struct Input
{
    TTree* oniaTree;
    TTree* hltanalysisTree;
    TTree* hltobjectTree;
    bool isL1;
    OniaType type;
};

const char oniaTreeName[]="hionia/myTree";
const char hltanalysisTreeName[] = "hltanalysis/HltTree";
const char hltobjDirectoryName[] = "hltobject/";

