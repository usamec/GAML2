#pragma once
//#include "common.h"
#include <array>
#include <vector>
#include "Sequence.h"
namespace dalign {
    extern "C" {
#include "DALIGN/align.h"
    }
}
using namespace std;

class Alignment {
    friend class DalignWrapper;
public:
    Alignment();
    Alignment(const Alignment& al);
    Alignment& operator=(const Alignment& al);
    virtual ~Alignment();
    int GetLengthOnA() const;
    int GetLengthOnB() const;
    pair<int, int> GetPosOnA() const;
    pair<int, int> GetPosOnB() const;
    pair<int, int> GetStartPos() const;
    pair<int, int> GetEndPos() const;
    bool ComputeTrace();
    bool GetAlignedPairs(vector<pair<int, int>>&pairs) const;
    bool PrintAlignment(const string& filename);
    bool GetCigarString(string &cigar) const;
    double GetSimilarity();
    int GetNumDifferences();
private:
    void Prepare(Sequence& _A, Sequence& _B, dalign::Work_Data* _workData, int _traceSpacing);
    void CopyData(const Alignment &al);
    dalign::Alignment alignment;
    dalign::Path path;
    dalign::Work_Data* workData;
    int traceSpacing;

    enum Status {
        AS_EMPTY,
        AS_PREPARED,
        AS_ALIGNMENT,
        AS_TRACE
    };

    Status status;
    bool deleteTrace;
};

class DalignWrapper {
public:
    virtual ~DalignWrapper();
    DalignWrapper();

    void SetAligningParameters(float corelation, int traceSpacing, const array<float, 4>& frequencies);
    bool ComputeAlignment(Sequence& A, Sequence& B, pair<int, int> seed, Alignment& al);
private:
    void FreeDalignData();
    dalign::Work_Data* workData;
    dalign::Align_Spec* alignSpec;
};