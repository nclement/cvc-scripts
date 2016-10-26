#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<vector>

#define BUFFER_SIZE 300

using namespace std;

class atom {
 private:
  double _x,_y,_z;
  unsigned int _id;
 public:
  atom(double x, double y, double z, unsigned int id) :
      _x(x), _y(y), _z(z), _id(id) {}

  unsigned int getID() {
    return _id;
  }

  double dist(atom& other) {
    double dist = pow( pow(_x-other._x, 2) + 
                       pow(_y-other._y, 2) +
                       pow(_z-other._z, 2), 0.5);

    return dist;
  }

  bool is_clash(atom& other, double clash_dist) {
    double dist = pow( pow(_x-other._x, 2) + 
                       pow(_y-other._y, 2) +
                       pow(_z-other._z, 2), 0.5);

    if (dist < clash_dist)
      return true;

    return false;
  }
};

void printUsage()
{
	cout<<"Usage:"<<endl;
	cout<<"clashDetector <PDB or PQR file>"<<endl;
}

int main(int argc, char *argv[])
{
	if(argc < 2)
	{
		printUsage();
		return -1;		
	}

	//reading pdb
  vector<atom> atoms;
	
	FILE* pdb = fopen(argv[1],"rt");
	if(!pdb)
	{
		cout<<"Could not open pdb file"<<endl;
		return -2;
	}

	vector<string> lines;
	char line[ BUFFER_SIZE+1 ];

  int j = 0;
	while (fgets(line, BUFFER_SIZE, pdb) != NULL)
	{
    if (strncmp(line, "ATOM", 4) && strncmp(line, "HETATM", 6))
    {
      continue;
    }
    double x, y, z;
    if (sscanf(line + 30, "%lf %lf %lf", &x, &y, &z) != 3)
    {
      fprintf(stderr, "\n\nError: Failed to read coordinates from line %d of input PDB file!\n\n", j);
      return -1;
    }
    unsigned int id;
    if (sscanf(line + 5, "%u", &id) != 1) {
      fprintf(stderr, "\n\nError: Failed to read coordinates from line %d of input PDB file!\n\n", j);
      return -1;
    }

    atoms.push_back(atom(x,y,z, id));
    ++j;
	}
	fclose(pdb);

  // Check for clashes
  double clash_dist = 1.0;
  int num_clashes = 0;
  double min_dist = 1.0;
  int min_i, min_j;
  for (int i = 0; i < atoms.size(); ++i) {
    for (int j = i+1; j < atoms.size(); ++j) {
      if (atoms[i].is_clash(atoms[j], clash_dist)) {
        num_clashes++;
        double cd = atoms[i].dist(atoms[j]);
        fprintf(stderr, "Clash between %d and %d: %f\n", 
                atoms[i].getID(), atoms[j].getID(), cd);
        if (cd < min_dist) {
          min_dist = cd;
          min_i = atoms[i].getID();
          min_j = atoms[j].getID();
        }
      }
    }
  }

  // Print stats
  fprintf(stdout, "Number of clashes: %d\n", num_clashes);
  fprintf(stdout, "Worst clash between %d and %d: %f\n", min_i, min_j, min_dist);
}
