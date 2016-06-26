#include "Sequence.h"
#include <algorithm>

char ToUpperCase(char c) {
    switch (c) {
        case 'a': case 'A': return 'A';
        case 'c': case 'C': return 'C';
        case 't': case 'T': return 'T';
        case 'g': case 'G': return 'G';
        default: return 0;
    }
}

char Remap(char c) {
    switch (c) {
        case 'a': case 'A': return 'T';
        case 'c': case 'C': return 'G';
        case 't': case 'T': return 'A';
        case 'g': case 'G': return 'C';
        default: return 0;
    }
}

Sequence::Sequence()
: dalignFromat(NULL) {
}

Sequence::Sequence(const Sequence& orig, bool remapReverse)
: id(orig.id), dalignFromat(NULL)  {
    if (remapReverse) {
        for (char c : orig.data) {
            data.push_back(Remap(c));
        }
        reverse(data.begin(), data.end());
    } else {
        data = orig.data;
    }
}

Sequence::~Sequence() {
    if (dalignFromat != NULL) {
        delete dalignFromat;
    }
}

Sequence::Sequence(const string& _data, string _id)
: data(_data), id(_id), dalignFromat(NULL) {

}

char* Sequence::ToDalignFromat() {
    if (dalignFromat != NULL) {
        return dalignFromat + 1;
    }
    int l = data.length();
    dalignFromat = new char[l + 2];
    dalignFromat[0] = 4;
    for (int i = 0; i < l; i++) {
        if (data[i] == 'A') dalignFromat[i + 1] = 0;
        if (data[i] == 'C') dalignFromat[i + 1] = 1;
        if (data[i] == 'G') dalignFromat[i + 1] = 2;
        if (data[i] == 'T') dalignFromat[i + 1] = 3;
    }
    dalignFromat[l + 1] = 4;
    return dalignFromat + 1;
}


FASTA::FASTA(const string& _filename)
: filename(_filename) {
}

FASTA& FASTA::operator>>(Sequence& seq) {
    ifstream is(filename);
    string data;
    string id;
    string line;

    if (getline(is, line)) {
        if (line.length()) {
            id = line.substr(1, line.length() - 1);
        }

        while (getline(is, line)) {
            for (char c : line) {
                data.push_back(ToUpperCase(c));
            }
        }
    }

    seq = Sequence(data, id);
    return *this;
}

FASTQ::FASTQ(const string& filename)
: is(filename.c_str()), isOk(true) {
}

FASTQ& FASTQ::operator>>(Sequence& seq) {
    string data;
    string id;
    string line;
    
    if (isOk) {
        isOk &= (bool)getline(is, line);
        if (line.length()) {
            id = line.substr(1, line.length() - 1);
        }
    }
    if (isOk) {
        isOk &= (bool)getline(is, line);
        for (char c : line) {
            data.push_back(ToUpperCase(c));
        }
    }
    if (isOk) {
        isOk &= (bool)getline(is, line);
    }
    if (isOk) {
        isOk &= (bool)getline(is, line);
    }
    if (isOk) {
        seq = Sequence(data, id);
    }
    else {
        seq = Sequence();
    }
    return *this;
}
