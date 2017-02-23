#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#define BUFFER_SIZE	1024
#define SQ(a) (a)*(a)

class Atom {
 public:
	Atom() : x_(0), y_(0), z_(0) {}
	Atom(int resi, float x, float y, float z, char* t) :
			resi_(resi), x_(x), y_(y), z_(z) {
    atomType_[0] = t[0];
    atomType_[1] = t[1];
    atomType_[2] = t[2];
    atomType_[3] = t[3];
  }

  int getResi() const { return resi_; }
	float getX() const {return x_;}
	float getY() const {return y_;}
	float getZ() const {return z_;}
  const char* getType() const {return atomType_;}
  bool isBackbone() {
    string atype = atomType_;
    if (atype == "CA" || atype == "O" || atype == "C" || atype == "N") {
      return true;
    }
    return false;
  }

 private:
	float x_, y_, z_;
  int resi_;
  char atomType_[4];
};

vector<Atom*> readPDBMolecule(string fn) {
	vector<Atom*> molecule;

	// Load each file
	FILE* inf = fopen(fn.c_str(), "r");
	if (!inf) {
		fprintf(stderr, "Could not load PDB file %f\n", fn.c_str());
		perror("Error");
		return molecule;
	}

	char atomLine[BUFFER_SIZE+10];
	float px, py, pz;
  int resi;
  char atomType[10];
	while (fgets(atomLine, BUFFER_SIZE, inf) != NULL) {
    if(sscanf(atomLine, "ATOM %*d %s %*s %*s %d %f %f %f", 
              atomType, &resi, &px, &py, &pz) != 4) {
      continue;
    }

		Atom* a = new Atom(resi, px, py, pz, atomType);
		molecule.push_back(a);
	}

	fclose(inf);

	return molecule;
}

float square(float a) {
  return a*a;
}

float RMSD_noca(const vector<Atom*> &a, const vector<Atom*> &b) {
  int aidx = 0;
  int bidx = 0;

  float rmsd = 0;
  int nCA = 0;

  while(aidx < a.size() && bidx < b.size()) {
    if (a[aidx]->getResi() < b[bidx]->getResi()) {
      aidx++; continue;
    }
    if (b[bidx]->getResi() < a[aidx]->getResi()) {
      bidx++; continue;
    }

    if (!a[aidx]->isBackbone()) {
      aidx++;
      continue;
    }
    if (!b[bidx]->isBackbone()) {
      bidx++;
      continue;
    }
      
    rmsd += square(a[aidx]->getX() - b[bidx]->getX()) 
        + square(a[aidx]->getY() - b[bidx]->getY()) 
        + square(a[aidx]->getZ() - b[bidx]->getZ());
    nCA++;
    aidx++;
    bidx++;
  }

  return sqrt(rmsd / nCA);
}

float RMSD(const vector<Atom*> &a, const vector<Atom*> &b) {
	assert(a.size() == b.size());
	float rmsd = 0;

  int nCA = 0;

	for (int i = 0; i < a.size(); i++) {
      
//    if (strncmp(a[i]->getType(), "CA", 2) == 0) {
      rmsd += square(a[i]->getX() - b[i]->getX()) + square(a[i]->getY() - b[i]->getY()) +
							square(a[i]->getZ() - b[i]->getZ());
//    }
    nCA++;
	}

	rmsd /= nCA;
	return sqrt(rmsd);
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
    fprintf(stderr, "Program will print the RMSD between two PDB files, but\n");
    fprintf(stderr, "without any movement, so the lines need to be identical\n");
    fprintf(stderr, "between both PDB files.\n\n");
		fprintf(stderr, "  usage: RMSDBasic <template.pdb> <target.pdb> [Ca only]\n");
		return -1;
	}

	vector<Atom*> tempMol = readPDBMolecule(argv[1]);
  vector<Atom*> targetMol = readPDBMolecule(argv[2]);

  float rmsd;
  if (argc == 4) {
    rmsd = RMSD_noca(tempMol, targetMol);
  } else {
    rmsd = RMSD(tempMol, targetMol);
  }
  printf("%f\n", rmsd);

  for (Atom* a : tempMol) {
    delete a;
  }
  for (Atom* a : targetMol) {
    delete a;
  }
}
