#ifndef IO_HPP
#define IO_HPP

#include <Eigen/Core>
#include <string>
#include <vector>

typedef std::vector<std::vector<std::string>> dataFrame;

/* load .csv file to dataFrame, implicitly in RowMajor */
dataFrame readCSV_DF(std::string path);

/* print subtable [start_row:row_start + rows),[start_col:col_start + cols) of
 * dataFrame */
void printDF(dataFrame df, int row_start, int rows, int col_start, int cols);

/* load subtable of dataFrame into double[][] in ColMajor */
double *df2arr(dataFrame df, int row_start, int rows, int col_start, int cols);

/* write ColMajor double[][] to .csv file with precision */
void writeCSV(std::string path, double arr[], int rows, int cols, int pres = 3);

/* load .csv file to Eigen::MatrixXd */
Eigen::MatrixXd readCSV(std::string path);

/* write Eigen::MatrixXd to .csv file */
void writeCSV(Eigen::MatrixXd mat, std::string path);

#endif  // IO_HPP
