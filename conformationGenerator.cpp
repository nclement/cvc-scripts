
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<vector>

#define BUFFER_SIZE 300

using namespace std;


void printUsage()
{
	cout<<"Usage:"<<endl;
	cout<<"conformationGenerator <PDB or PQR file> <xformfile> <startConf> <endConf>"<<endl;
}

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		printUsage();
		return -1;		
	}


	//reading pdb
	
	FILE* pdb = fopen(argv[1],"rt");
	if(!pdb)
	{
		cout<<"Could not open pdb file"<<endl;
		return -2;
	}

	vector<string> lines;
	char line[ BUFFER_SIZE+1 ];
	int numLines = 0;

	while (fgets(line, BUFFER_SIZE, pdb) != NULL)
	{
		string l(line);
		lines.push_back(l);
		numLines++;
	}
	fclose(pdb);


	//reading xformations and writing output conformations

	FILE* xforms = fopen(argv[2],"rt");
	if(!xforms)
	{
		cout<<"Could not open transformations file"<<endl;
		return -3;
	}


	int start = atoi(argv[3]);
	int end = atoi(argv[4]);
	int numXform;
	fscanf(xforms, "%d", &numXform);
	if(end > numXform) end = numXform;

	for(int i=0; i<=end; i++)
	{
		double t00, t01, t02, t03;
		double t10, t11, t12, t13;
		double t20, t21, t22, t23;
		double t30, t31, t32, t33;
		double r, s;
		fscanf(xforms, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t00, &t01, &t02, &t03, &t10, &t11, &t12, &t13, &t20, &t21, &t22, &t23, &t30, &t31, &t32, &t33);
		printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23, t30, t31, t32, t33);

//		fscanf(xforms, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t00, &t01, &t02, &t03, &t10, &t11, &t12, &t13, &t20, &t21, &t22, &t23);
//		printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23);		


		if(i<start) continue;
	
		string prefix(argv[1]);
		char num[6];
		sprintf(num, "%d", i);
		string newPDB = prefix.substr(0, prefix.length() - 4) + "_conf_" + num + ".pdb";
		FILE* newpdb = fopen(newPDB.c_str(),"wt");

		if(!newpdb)
		{
			cout<<"Failed to create new pdb "<<newPDB<<endl;
			return -4;
		}

		for(int j=0; j<numLines; j++)
		{
			char currentLine[500];
			strcpy(currentLine, lines[j].c_str());

			if (strncmp(currentLine, "ATOM", 4) && strncmp(currentLine, "HETATM", 6))
			{
				fprintf(newpdb, "%s", currentLine);
				continue;
			}
			double x, y, z;
			if (sscanf(currentLine + 30, "%lf %lf %lf", &x, &y, &z) != 3)
			{
				fprintf(stderr, "\n\nError: Failed to read coordinates from line %d of input PDB file!\n\n", j);
				return false;
			}

			double nx = t00 * x + t01 * y + t02 * z + t03;
			double ny = t10 * x + t11 * y + t12 * z + t13;
			double nz = t20 * x + t21 * y + t22 * z + t23;

			currentLine[ 30 ] = 0;

			fprintf(newpdb, "%s%7.3lf %7.3lf %7.3lf %s", currentLine, nx, ny, nz, currentLine + 54);
		}
	
		fclose(newpdb);
	}

	fclose(xforms);
}
