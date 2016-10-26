
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<vector>

#define BUFFER_SIZE 300

using namespace std;

/*
 * This program takes an input xformfile (12 numbers from a 3x4 matrix) and 
 * applies them all successively to a single PDB file.
 */

void printUsage()
{
	cout<<"Usage:"<<endl;
	cout<<"conformationGenerator_single <PDB or PQR file> <xformfile> <output> [VERBOSE]"<<endl;
}

struct Atom {
  Atom(float x, float y, float z, string l) : _x(x), _y(y), _z(z), _line(l) {}
	float _x, _y, _z;
  string _line;
};

vector<Atom*> readPDBMolecule(string fn) {
	vector<Atom*> molecule;

	// Load each file
	FILE* inf = fopen(fn.c_str(), "r");
	if (!inf) {
		fprintf(stderr, "Could not load PDB file %s\n", fn.c_str());
		perror("Error");
		return molecule;
	}

	char atomLine[BUFFER_SIZE+10];
	float px, py, pz;
	while (fgets(atomLine, BUFFER_SIZE, inf) != NULL) {
    if(sscanf(atomLine, "ATOM %*d %*s %*s %*s %*d %f %f %f", 
              &px, &py, &pz) != 3) {
      continue;
    }

		Atom* a = new Atom(px, py, pz, atomLine);
		molecule.push_back(a);
	}

	fclose(inf);

	return molecule;
}

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		printUsage();
		return -1;		
	}

  bool verbose = false;
  if (argc == 5) {
    verbose = true;
  }


	//reading pdb
  vector<Atom*> pdb = readPDBMolecule(argv[1]);
	
	if(pdb.size() == 0)
	{
		cout << "Could not open pdb file" << endl;
		return -2;
	}

	//reading xformations and writing output conformations

	FILE* xforms = fopen(argv[2],"rt");
	if(!xforms)
	{
		cout<<"Could not open transformations file"<<endl;
		return -3;
	}


	while (!feof(xforms)) {
		double t00, t01, t02, t03;
		double t10, t11, t12, t13;
		double t20, t21, t22, t23;
		fscanf(xforms, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t00, &t01, &t02, &t03, &t10, &t11, &t12, &t13, &t20, &t21, &t22, &t23);
		if (verbose) printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23);

//		fscanf(xforms, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t00, &t01, &t02, &t03, &t10, &t11, &t12, &t13, &t20, &t21, &t22, &t23);
//		printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23);		


    // Multiply all the atoms by this transformation
		for(unsigned int j=0; j<pdb.size(); j++)
		{
      Atom* a = pdb[j];
			double nx = t00 * a->_x + t01 * a->_y + t02 * a->_z + t03;
			double ny = t10 * a->_x + t11 * a->_y + t12 * a->_z + t13;
			double nz = t20 * a->_x + t21 * a->_y + t22 * a->_z + t23;
      a->_x = nx;
      a->_y = ny;
      a->_z = nz;
		}
	
	}

	fclose(xforms);

  FILE* ofs = fopen(argv[3], "w");
  // Now write out the PDB
  for (unsigned int i = 0; i < pdb.size(); i++) {
    string currentLine = pdb[i]->_line;
    currentLine[ 30 ] = 0;
    fprintf(ofs, "%s%8.2f%8.2f%8.2f%s", 
            currentLine.c_str(), pdb[i]->_x, pdb[i]->_y, pdb[i]->_z, currentLine.c_str() + 54);
  }
  fclose(ofs);
}
