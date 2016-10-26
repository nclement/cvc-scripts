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
	Atom(float x, float y, float z) :
			x_(x), y_(y), z_(z) {}

	float getX() const {return x_;}
	float getY() const {return y_;}
	float getZ() const {return z_;}

 private:
	float x_, y_, z_;
};

vector<vector<Atom*>> readPDBQTMolecules(string fn, bool isPDB) {
	vector<vector<Atom*>> molecules;

	// Load each file
	FILE* inf = fopen(fn.c_str(), "r");
	if (!inf) {
		fprintf(stderr, "Could not load PDB file %f\n", fn.c_str());
		perror("Error");
		return molecules;
	}

	vector<Atom*> molecule;

	char atomLine[BUFFER_SIZE+10];
	float px, py, pz;
	while (fgets(atomLine, BUFFER_SIZE, inf) != NULL) {
		if (strncmp(atomLine, "ENDMDL", 6) == 0) {
			//printf("Found new molecule!\n");
			molecules.push_back(molecule);
			molecule.clear();
		}

    if (isPDB) {
      if(sscanf(atomLine, "ATOM %*d %*s %*s %*s %*d %f %f %f", &px, &py, &pz) != 3) {
        continue;
      }
    } else {
      if(sscanf(atomLine, "ATOM %*d %*s %*s %*d %f %f %f", &px, &py, &pz) != 3) {
        continue;
      }
    }

		Atom* a = new Atom(px, py, pz);
		molecule.push_back(a);
	}
	// Add it if there are any left
	if (molecule.size())
		molecules.push_back(molecule);

	fclose(inf);

	return molecules;
}

float RMSD(const vector<Atom*> &a, const vector<Atom*> &b) {
	assert(a.size() == b.size());
	float rmsd = 0;

	for (int i = 0; i < a.size(); i++) {
		rmsd += SQ(a[i]->getX() - b[i]->getX()) + SQ(a[i]->getY() - b[i]->getY()) +
							SQ(a[i]->getZ() - b[i]->getZ());
	}

	rmsd /= a.size();
	return sqrt(rmsd);
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		fprintf(stderr, "usage: RMSDMulti <template.pdbqt> <target_list of pdbqt> [pdb=1]\n");
		return -1;
	}

  bool isPDB = false;
  if (argc == 4 && atoi(argv[3]) == 1) {
    isPDB = true;
  }

	vector<Atom*> tempMol = readPDBQTMolecules(argv[1], isPDB)[0];

	// Open the PDBs file
  ifstream pdbs;
  pdbs.open(argv[2]);
  if (pdbs.is_open()) {
    string line;
    while(getline (pdbs, line)) {
			vector<vector<Atom*>> targetMols = readPDBQTMolecules(line, isPDB);

			for (const vector<Atom*> &m : targetMols) {
				float rmsd = RMSD(tempMol, m);
				fprintf(stdout, "%f\t", rmsd);
			}
			fprintf(stdout, "\n");

			// Clean up 
			for (vector<Atom*> &m : targetMols) {
				for(Atom* a : m) {
					delete a;
				}
			}
		}
	} else {
		perror("Could not open target_list file!");
	}
}
