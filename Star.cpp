// Student: Alexander Gómez Del Carpio
// Subject: Computational Molecular Biology
// Description: In this piece of code we can see a use of dynamic programming in
// order to get the alignment of sequences in DNA or RNA with Needleman-Wunsch
// algorithm, whether you need the best result or all possible results. Date:
// 04/04/2024

#include <fstream>
#include <iostream>
#include <limits.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm> // To use remove_if

using namespace std;

// Little class to define a cell
class Cell {
public:
  int score;      // Score depending of previous
  vector<int> xy; // Store indexes of previous cells with the greater score
  Cell() : score(0), xy(0, 0) {}
};

void printMatrix(vector<vector<Cell>> &M) { // Prints the matrix of DP
  for (auto i : M) {
    for (auto j : i) {
      if (j.score >= 0)
        cout << " ";
      cout << j.score << " ";
    }
    cout << endl;
  }
}

void initialFill(
    vector<vector<Cell>> &M) { // Initialize the matrix first column and row
  for (int i = 1; i < M.size(); i++) {
    M[i][0].score = i * -2;
    M[i][0].xy = {i - 1, 0};
  }
  for (int i = 1; i < M[0].size(); i++) {
    M[0][i].score = i * -2;
    M[0][i].xy = {0, i - 1};
  }
}

// Calculate score and previous cell coordinates with best result
pair<int, vector<int>> calculateScore(const vector<vector<Cell>> &M, int i,
                                      int j, int diff_score) {
  int scores[3];
  scores[0] = M[i - 1][j - 1].score + diff_score; // 0 -> Diagonal
  scores[1] = M[i - 1][j].score - 2;              // 1 -> Up
  scores[2] = M[i][j - 1].score - 2;              // 2 -> Left

  vector<int> paths;     // Stores 0,1,2 according best score
  vector<int> paths_idx; // Stores coordinates previous cells with best score
  int max = INT_MIN;

  for (int k = 0; k < 3; k++) { // Get cell with biggest score
    if (scores[k] > max) {
      max = scores[k];
      paths.clear();
      paths.push_back(k);
    } else if (scores[k] == max) {
      paths.push_back(k);
    }
  }

  for (auto k : paths) { // Give the coordinates of the best score cell
    if (k == 0) {
      paths_idx.push_back(i - 1);
      paths_idx.push_back(j - 1);
    } else if (k == 1) {
      paths_idx.push_back(i - 1);
      paths_idx.push_back(j);
    } else {
      paths_idx.push_back(i);
      paths_idx.push_back(j - 1);
    }
  }

  return make_pair(max, paths_idx);
}

// Global variables in order to not use too many arguments
int minScore = INT_MAX; // Counts how many breaks for gaps occur, less better
pair<string, string> optimas = {"",""}; // Store optimal solution of two sequences
unordered_map<string, int> visited;  // In case the path is already visited with its previous score
string state; // Which coordinate we visit to store in visited
pair<bool, bool> gapAtTheEnd = {true, true}; // In case the sequences have gaps at the end in order to not
                 // count them as ruptures

// Function that finds the optimal solution
void rebuiltPath(const vector<vector<Cell>> &M, const string &a, const string &b, int i, int j, char ins1, char ins2, vector<string> ans, int scoreGaps) {
  ans[0].insert(ans[0].begin(), ins1); // Insert solution at the beginning
  ans[1].insert(ans[1].begin(), ins2);

  if (ins1 == '-' && !gapAtTheEnd.first) scoreGaps++; // In case if there's not gap at the end anymore
  else if (ins1 != '-'&&ins1!=' ') gapAtTheEnd.first = false;
  if (ins2 == '-' && !gapAtTheEnd.second) scoreGaps++;
  else if (ins2 != '-'&&ins2!=' ') gapAtTheEnd.second = false;
  // Check if is visited and store the new path in visited
  string state = to_string(i) + "," + to_string(j);

  if (visited.find(state) != visited.end()) { // We ignore the visited path only if has a lower or equal score in the same point
    if (scoreGaps > visited[state]) return;
  }

  visited[state] = scoreGaps;

  if (i == 0 && j == 0) { // In case the path arrives to the first cell (the
                          // end), verify if it is the best
    if (ans[0][0] == '-') scoreGaps--;
    if (ans[1][0] == '-') scoreGaps--;
    if (scoreGaps <= minScore) {
      minScore = scoreGaps;
      optimas.first = ans[0];
      optimas.second = ans[1];
    }

    return;
  }

  for (int k = 0; k < M[i][j].xy.size();k += 2) { // See what continues in the path
    ins1 = a[i]; // Ins1 and 2 next insertions
    ins2 = b[j];
    if (M[i][j].xy[k] == i) ins1 = '-'; // In case is vertical or horizontal there's a gap
    else if (M[i][j].xy[k + 1] == j) ins2 = '-';
    rebuiltPath(M, a, b, M[i][j].xy[k], M[i][j].xy[k + 1], ins1, ins2, ans, scoreGaps);
  }
}


pair<int, vector<string>> NW_Algorithm(string a, string b) {
	a="-"+a;
	b="-"+b;
  // Matrix init
  vector<vector<Cell>> matrix(a.size(), vector<Cell>(b.size(), Cell()));
  initialFill(matrix);

  int diff_score; // If we come from diagonal cell, how much is the penalty for
                  // being different
  pair<int, vector<int>> receiveScore; // Variable receives the result of every cell.

  // We go from second column and row and see what is the previous cell with the
  // best score
  for (int i = 1; i < a.size(); i++) {
    for (int j = 1; j < b.size(); j++) {
      a[i] == b[j] ? diff_score = 1 : diff_score = -1;
      receiveScore = calculateScore(matrix, i, j, diff_score);
      matrix[i][j].score = receiveScore.first;

      matrix[i][j].xy = receiveScore.second;
    }
  }

  vector<string> initial_ans(2, ""); // Store results during recursion

  rebuiltPath(matrix, a, b, matrix.size() - 1, matrix[0].size() - 1, ' ', ' ',
              initial_ans, 0);

  vector<string> result_alignment{optimas.first,optimas.second};

  //Restore to initial value
  minScore = INT_MAX;
  optimas = {"",""};
  visited.clear();
  gapAtTheEnd = {true, true};

  return make_pair(matrix[a.size() - 1][b.size() - 1].score,
                   result_alignment);
}

void star_alignment(vector<string> sequences){
  int max_score = INT_MIN; // Store max score of a sequence aligned with all the rest
  int temp_score = 0;
  int best_sequence_idx;
  pair<int,vector<string>> temp_NW;
  
  // Find the center with the best score
  for (int i = 0; i < sequences.size(); i++) {
    for (int j = 0; j < sequences.size(); j++) {
      if (i != j) {
        temp_NW = NW_Algorithm(sequences[i], sequences[j]);
        temp_score += temp_NW.first;
      }
    }
    if (temp_score > max_score) {
      max_score = temp_score;
      best_sequence_idx = i;
    }
    temp_score = 0;
  }

  // Find aligned sequence with the best one
  vector<string> new_aligned_sequences;
  string newSeq;

  for (int i = 0; i < sequences.size(); i++) {
    if (i == best_sequence_idx)
      new_aligned_sequences.emplace_back(sequences[best_sequence_idx]); //Insert the best sequence bc is the same
    else{
      temp_NW=NW_Algorithm(sequences[best_sequence_idx], sequences[i]);
      newSeq=temp_NW.second[1];
      newSeq.pop_back();      //Erasing the blank space at the end of the string because of NW_Alogorithm
      new_aligned_sequences.emplace_back(newSeq);
    }
  }

  int max_size = 0;   //Find out the maximum size among the sequences
  for (auto i : new_aligned_sequences) {
    if (i.size() > max_size) max_size = i.size();
  }

  for(auto i:new_aligned_sequences){      //Add - at the end until it reaches the same lenght of the largest one
		while(i.size()<max_size) i=i+"-";
		cout<<i<<endl;
  }
  return;
}


vector<string> readADNtxt(const string& type = "") {
    ifstream file("BRCA1.txt");
    vector<string> sequences;

    if (!file.is_open()) {
        cerr << "Failed to open file." << endl;
        return sequences; // Return empty vector if couldn't access
    }

    string line;
    bool includeAll = type.empty();
    bool includeF = (type == "F");
    bool includeR = (type == "R");

    while (getline(file, line)) {
        if (line.empty() || line.size() < 4) continue; 

        // See if the line has F forward or R reverse
        if (isalpha(line[9])) {
            // Get type (F o R)
            char seqType = line[9];
            if ((includeAll || (seqType == 'F' && includeF) || (seqType == 'R' && includeR)) && line.size() > 6) {
                // Find "-"'s positions
                size_t startPos = line.find('-');
                size_t endPos = line.find('-', startPos + 1);
                
                // Get sequence between "-"
                if (startPos != string::npos && endPos != string::npos && endPos > startPos) {
                    string sequence = line.substr(startPos + 1, endPos - startPos - 1);
                    
                    sequences.push_back(sequence);
                }
            }
        }
    }

    file.close();
    return sequences;
}

int main() {
  //Insert sequences manually
  vector<string> sequences{"ATTGCCATT", "ATGGCCATT", "ATCCAATTTT", "ATCTTCTT",
                           "ACTGACC"};
  //From the txt (parameter empty is forward an reverse,"F" forward, "R" reverse)
  sequences=readADNtxt();
  star_alignment(sequences);

  return 0;
}