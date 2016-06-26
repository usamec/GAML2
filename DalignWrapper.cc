#include <string.h>

#include "DalignWrapper.h"

DalignWrapper::DalignWrapper()
: workData(NULL), alignSpec(NULL) {
}

DalignWrapper::~DalignWrapper() {
    FreeDalignData();
}

void DalignWrapper::SetAligningParameters(float corelation, int traceSpacing, const array<float, 4>& frequencies) {
    FreeDalignData();
    float freq[4] = {frequencies[0], frequencies[1], frequencies[2], frequencies[3]};
    alignSpec = dalign::New_Align_Spec(corelation, traceSpacing, freq);
    workData = dalign::New_Work_Data();
}

bool DalignWrapper::ComputeAlignment(Sequence& A, Sequence& B, pair<int, int> seed, Alignment& al) {
    if (workData == NULL || alignSpec == NULL) {
        return false;
    }
    al.Prepare(A, B, workData, dalign::Trace_Spacing(alignSpec));
    dalign::Local_Alignment(&(al.alignment), workData, alignSpec, seed.first, seed.first, seed.second);
    al.status = Alignment::AS_ALIGNMENT;
    return true;
}

void DalignWrapper::FreeDalignData() {
    if (workData != NULL) {
        dalign::Free_Work_Data(workData);
        workData = NULL;
    }
    if (alignSpec != NULL) {
        dalign::Free_Align_Spec(alignSpec);
        alignSpec = NULL;
    }
}

Alignment::Alignment()
: status(AS_EMPTY), deleteTrace(false) {
}

Alignment::Alignment(const Alignment& al) {
    CopyData(al);
}

Alignment& Alignment::operator=(const Alignment& al) {
    CopyData(al);
    return *this;
}

void Alignment::CopyData(const Alignment& al) {
    status = al.status;
    deleteTrace = false;
    if (al.status >= AS_PREPARED) {
        workData = al.workData;
        traceSpacing = al.traceSpacing;
        alignment.path = &path;
        alignment.aseq = al.alignment.aseq;
        alignment.alen = al.alignment.alen;
        alignment.bseq = al.alignment.bseq;
        alignment.blen = al.alignment.blen;
        alignment.flags = al.alignment.flags;
    }
    
    if (al.status == AS_TRACE) {
        path.abpos = al.path.abpos;
        path.aepos = al.path.aepos;
        path.bbpos = al.path.bbpos;
        path.bepos = al.path.bepos;
        path.diffs = al.path.diffs;
        path.tlen = al.path.tlen;
        
        deleteTrace = true;
        path.trace = new int[path.tlen];
        memcpy(path.trace, al.path.trace, path.tlen * sizeof(int));
    }
}



Alignment::~Alignment() {
    if (deleteTrace) {
        delete [] (int*)path.trace;
    }
}


void Alignment::Prepare(Sequence& _A, Sequence& _B, dalign::Work_Data* _workData, int _traceSpacing) {
    status = AS_PREPARED;
    workData = _workData;
    traceSpacing = _traceSpacing;
    alignment.path = &path;
    alignment.aseq = _A.ToDalignFromat();
    alignment.alen = _A.GetData().length();
    alignment.bseq = _B.ToDalignFromat();
    alignment.blen = _B.GetData().length();
    alignment.flags = 0;
}

int Alignment::GetLengthOnA() const {
    if (status < AS_ALIGNMENT) {
        return -1;
    }

    return path.aepos - path.abpos;
}

int Alignment::GetLengthOnB() const {
    if (status < AS_ALIGNMENT) {
        return -1;
    }

    return path.bepos - path.bbpos;
}

bool Alignment::ComputeTrace() {
    if (status < AS_ALIGNMENT) {
        return false;
    }

    dalign::Compute_Trace_MID(&alignment, workData, traceSpacing);
    status = AS_TRACE;
    return true;
}

bool Alignment::PrintAlignment(const string& filename) {
    if (status < AS_TRACE) {
        return false;
    }

    FILE *alignmentFile = fopen(filename.c_str(), "w");
    Print_Alignment(alignmentFile, &alignment, workData, 10, 80, 5, 1, 10);
    return true;
}

bool Alignment::GetAlignedPairs(vector<pair<int, int> >& pairs) const {
    pairs.clear();
    if (status < AS_TRACE) {
        return false;
    }
    int p, c;
    int i = path.abpos + 1;
    int j = path.bbpos + 1;
    // copy-pasterino
    for (c = 0; c < path.tlen; c++)
        if ((p = ((int*) path.trace)[c]) < 0) {
            p = -p;
            while (i != p) {
                pairs.push_back(make_pair(i - 1, j - 1));
                i += 1;
                j += 1;
            }
            j += 1;
        } else {
            while (j != p) {
                pairs.push_back(make_pair(i - 1, j - 1));
                i += 1;
                j += 1;
            }
            i += 1;
        }
    p = path.aepos;
    while (i <= p) {
        pairs.push_back(make_pair(i - 1, j - 1));
        i += 1;
        j += 1;
    }

    return true;
}

bool Alignment::GetCigarString(string& cigar) const {
    cigar.clear();
    if (status < AS_TRACE) {
        return false;
    }

    int p, c;
    int i = path.abpos + 1;
    int j = path.bbpos + 1;
    int counter = 0;
    // copy-pasterino
    for (c = 0; c < path.tlen; c++) {
        if (c > 0 && ((int*) path.trace)[c] != ((int*) path.trace)[c - 1]) {
            cigar.append(to_string(counter) + (((int*) path.trace)[c - 1] < 0 ? "I" : "D"));
        }
        if ((p = ((int*) path.trace)[c]) < 0) {
            p = -p;
            int Ms = 0;
            while (i != p) {
                Ms++;
                i += 1;
                j += 1;
            }
            if (Ms) {
                cigar.append(to_string(Ms) + "M");
                counter = 0;
            }
            j += 1;
        } else {
            int Ms = 0;
            while (j != p) {
                Ms++;
                counter = 0;
                i += 1;
                j += 1;
            }
            if (Ms) {
                cigar.append(to_string(Ms) + "M");
                counter = 0;
            }
            i += 1;
        }
        counter++;
    }
    if (counter) {
        cigar.append(to_string(counter) + (((int*) path.trace)[path.tlen - 1] < 0 ? "I" : "D"));
    }

    int Ms = path.aepos - i + 1;
    if (Ms) {
        cigar.append(to_string(Ms) + "M");
    }
    return true;
}

pair<int, int> Alignment::GetPosOnA() const {
    if (status < AS_ALIGNMENT) {
        return make_pair(0, 0);
    }
    
    return make_pair(path.abpos, path.aepos);
}

pair<int, int> Alignment::GetPosOnB() const {
    if (status < AS_ALIGNMENT) {
        return make_pair(0, 0);
    }
    
    return make_pair(path.bbpos, path.bepos);
}

pair<int, int> Alignment::GetStartPos() const {
    if (status < AS_ALIGNMENT) {
        return make_pair(0, 0);
    }
    
    return make_pair(path.abpos, path.bbpos);
}

pair<int, int> Alignment::GetEndPos() const {
    if (status < AS_ALIGNMENT) {
        return make_pair(0, 0);
    }
    
    return make_pair(path.aepos, path.bepos);
}

double Alignment::GetSimilarity() {
    return 0;
}

int Alignment::GetNumDifferences() {
    if (status < AS_TRACE) {
        return 0;
    }

    return path.diffs;
}
