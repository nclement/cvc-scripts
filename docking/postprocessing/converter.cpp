#include <algorithm>    // reverse
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

vector<float> splitLine(string &l) {
  vector<float> ret;

  stringstream iss(l);
  do {
    float t;
    iss >> t;
    ret.push_back(t);
  } while(iss);

  return ret;
}

vector<float> getMatrix(string &l) {
  vector<float> matrix;

  // Get the vector we need
  vector<float> cells = splitLine(l);
  for (int i = 12; i > 0; --i) {
    matrix.push_back(cells[cells.size() - i - 2 - 1]);
  }

  return matrix;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    fprintf(stderr, "convert <F2Dock output file> <transformations file>\n");
    return -1;
  }

  ifstream f2if(argv[1]);

  if (!f2if.good()) {
    fprintf(stderr, "Error: Could not open input file, %s\n", argv[1]);
    return -1;
  }

  vector<vector<float> > lines;
  string line;
  bool reading = false;
  while(getline(f2if, line)) {
    if (line == "# START PEAKS") {
      reading = true;
      continue;
    } else if (line == "# END PEAKS") {
      break;
    }
    if (!reading) {
      continue;
    }

    lines.push_back(getMatrix(line));
  }
  f2if.close();

  // Need to reverse order this
  reverse(lines.begin(), lines.end());

  // Open the new file
  ofstream xof(argv[2]);
  if (!xof.good()) {
    fprintf(stderr, "Output file no good: %s\n", argv[2]);
    return -1;
  }

  // First write the number of transforms
  xof << lines.size() << "\n";
  // Then write the reversed transforms
  for (int i = 0; i < lines.size(); ++i) {
    for (int j = 0; j < lines[i].size() - 1; ++j) {
      xof << lines[i][j] << " ";
    }
    // write the last one with a newline instead of a space
    xof << lines[i][lines[i].size() - 1] << "\n";
  }

  xof.close();
  return 0;
}
