#include "io.hpp"
#include <Eigen/Core>
#include <fstream>   // std::fileStream()
#include <iomanip>   // std::setw(), std::setprecision()
#include <iostream>  // std::cout
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;
typedef vector<vector<string>> dataFrame;

// load .csv file to dataFrame, implicitly in RowMajor
dataFrame readCSV_DF(string path) {
  /* Inspired by:
  https://stackoverflow.com/questions/34218040/how-to-read-a-csv-file-data-into-an-array
  http://forums.codeguru.com/showthread.php?462383-2D-Vector-of-Strings!-How-It-Can-Be-Done-!!!
*/

  dataFrame table;
  ifstream fileStream(path);  // open input fileStream

  if (!fileStream.is_open()) {
    cout << "Could not open \"" << path << "\"" << endl;
    return table;
  }

  string line;
  while (getline(fileStream, line)) {
    vector<string> row;
    stringstream lineStream(line);

    string token;
    while (getline(lineStream, token, ',')) row.push_back(token);

    table.push_back(row);
  }

  return table;
}

int inline min(int a, int b) { return (a < b) ? a : b; }
int inline max(int a, int b) { return (a > b) ? a : b; }

// print subtable [start_row:row_start + rows),[start_col:col_start + cols) of
// dataFrame
void printDF(dataFrame df, int row_start, int rows, int col_start, int cols) {
  int row_end = row_start + rows;
  int col_end = col_start + cols;

  // round range to avoid index out of bounds
  row_start = max(row_start, 0);
  col_start = max(col_start, 0);
  row_end = min(df.size(), row_end);
  col_end = min(df[0].size(), col_end);

  int width = 12;  // column width, RowMajor
  for (int i = row_start; i < row_end; i++) {
    for (int j = col_start; j < col_end; j++) {
      cout << setw(width) << df[i][j];
    }
    cout << "\n";
  }

  cout << endl;
}

// load subtable of dataFrame into double[][] in ColMajor
double *df2arr(dataFrame df, int row_start, int rows, int col_start, int cols) {
  int row_end = row_start + rows;
  int col_end = col_start + cols;

  // round range to avoid index out of bounds
  row_start = max(row_start, 0);
  col_start = max(col_start, 0);
  row_end = min(df.size(), row_end);
  col_end = min(df[0].size(), col_end);

  int m = row_end - row_start;
  int n = col_end - col_start;

  double *arr = new double[m * n];  // in ColMajor
  for (int j = col_start; j < col_end; j++) {
    for (int i = row_start; i < row_end; i++) {
      arr[rows * j + i] = stof(df[i][j]);  // string to double
    }
  }

  return arr;
}

// write double[][] to .csv file with precision
void writeCSV(string path, double arr[], int rows, int cols, int pres) {
  ofstream fileStream(path);  // open output fileStream

  if (!fileStream.is_open()) {
    cout << "Could not open \"" << path << "\"" << endl;
    return;
  }

  // print implicitly in RowMajor, access in ColMajor
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      fileStream << fixed << setprecision(pres) << arr[rows * j + i] << ", ";
    }
    fileStream << "\r";
  }

  return;
}

MatrixXd readCSV(string path) {
  dataFrame df = readCSV_DF(path);
  int rows = df.size();
  int cols = df[0].size();
  return Map<MatrixXd>(df2arr(df, 0, rows, 0, cols), rows, cols);
}

void writeCSV(MatrixXd mat, string path) {
  writeCSV(path, mat.data(), mat.rows(), mat.cols());
}