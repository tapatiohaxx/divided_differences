/*****************************************************************************************************************************************************************************
 * Name: Matthew Plascencia                                                                                                                                                  *
 * ID #: 012600809                                                                                                                                                           *
 * CS 3010                                                                                                                                                                   *
 * Section #: 2                                                                                                                                                              *
 * Project #: 4                                                                                                                                                              *
 * Description: This project will implement Newton's Divided Differences method and use it to find interpolating                                                             *
 *              polynomials in lagrange and simplified forms.                                                                                                                *
 *****************************************************************************************************************************************************************************/
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <fstream>
#include <ios>

int elements_size; //global variable to store the number of rows
std::vector<float> coefficients; //The coefficients that we use for the 
float difference_matrix[50][51]; //matrix used to hold the difference table

void inputFromConsole(){	
	std::cout << "Enter the total number of set of points : ";
	std::cin >> elements_size;
	for (int i = 1; i <= elements_size; ++i)
	{
		std::cout << "\nX" << i << " = ";
		std::cin >> difference_matrix[i][0];
		std::cout << "Y" << i << " = ";
		std::cin >> difference_matrix[i][1];
	}
}
void inputFromFile(){
	std::vector<float> file_values;
	float input;
	std::string file_name;
	std::cout << "Enter a file name to fetch values from: \n";
	std::cin.get();
	getline(std::cin, file_name);
	std::ifstream table_input(file_name);
	while (table_input >> input) 
		file_values.push_back(input);
	table_input.close();
	if (file_values.size() > 100) {
		std::cout << "TOO MANY VALUES!!! EXITING.....";
		exit(0);
	}
	if (file_values.size() < 6){
		std::cout << "TOO FEW VALUES!!! EXITING.....";
		exit(0);
	}
	if (file_values.size() % 2 != 0) {
		std::cout << "UNEVEN AMOUNT OF VALUES IN FILE. EXITING....";
		exit(0);
	}
	for (int i = 0; i < file_values.size()/2; ++i)
		difference_matrix[i + 1][0] = file_values[i];
	for(int i = file_values.size()/2; i < file_values.size(); ++i)
		difference_matrix[i - file_values.size() / 2 + 1][1] = file_values[i];
	elements_size = file_values.size() / 2;
}

void menuandInputs() {
	int input;
	std::cout << "Would you like to interpolate a polynomial from the console or from a file? \n"
			<< "1: FROM FILE \n" <<  "2: FROM CONSOLE \n\n";
	std::cin >> input;
	if (input == 1) 
		inputFromFile();
	else if (input == 2)
		inputFromConsole();	
	else
		exit(0);
}
//create and display the difference matrix and display it
void makeDifferenceTable(){
	for (int j = 1; j < elements_size; ++j){
		for (int i = 1; i <= elements_size - j; ++i){
			if (j == 1)
				difference_matrix[i][j + 1] = (difference_matrix[i + 1][j] - difference_matrix[i][j]) / (difference_matrix[i + 1][0] - difference_matrix[i][0]);
			else
				difference_matrix[i][j + 1] = (difference_matrix[i + 1][j] - difference_matrix[i][j]) / (difference_matrix[i + j][0] - difference_matrix[i][0]);
		}
	}
	std::cout << "\nDIVIDED DIFFERENCE TABLE \n\n";
	std::cout << "A number and a p (i p) denotes a number of pairs past three pairs.\n";
	std::cout << std::setw(10) << "\t  x\t  ";
	std::cout << std::setw(10) << " f[]\t  ";	
	for (int i = 1; i < elements_size; ++i){
		if (i == 1)
			std::cout << std::setw(10) << " f[,]\t ";
		else if (i == 2)
			std::cout << std::setw(10) << " f[,,]\t ";
		else if (i == 3)
			std::cout << std::setw(10) << " f[,,,]\t ";
		else
			std::cout << std::setw(9) << "f[" << i << "]\t";
	}
	std::cout << "\n\n";
	int cols = 1;
	for (int i = 1; i <= elements_size; ++i){
		std::cout << std::setw(10) << difference_matrix[i][0] << std::setprecision(5) << "\t";
		std::cout << std::setw(10) << difference_matrix[i][1] << std::setprecision(5) << "\t";
		for (int j = 2; j <= elements_size -i  + 1; ++j){
			/*if (j > 1) {
				for (int k = 1; k < j; ++k)
					std::cout << "\n";
				std::cout << std::setw(10) << difference_matrix[i][j] << std::fixed << std::setprecision(5) << "\t";
			}
			else*/
				std::cout << std::setw(10) << difference_matrix[i][j] << std::setprecision(5) << "\t";			
		}
		std::cout << "\n\n";
	}
}
//Constructs the large polynomial by looking at the first row of the difference table and the first column and using them accordingly.
void makeNewton(){
	std::cout << "\n\nUnNSIMPLIFIED POLYNOMIAL:\n\n";
	coefficients.push_back(difference_matrix[1][1]);
	for (int f = 2; f <= elements_size; ++f)
		coefficients.push_back(difference_matrix[1][f]);
	for (int i = 0; i < coefficients.size(); ++i)
	{
		if (i == 0)
			std::cout << coefficients[i];
		else{
			//print out the coefficient before the factored part of the polynomial.
			if (coefficients[i] > 0)
				std::cout << " + " << coefficients[i];
			else if (coefficients[i] < 0)
				std::cout << " - " << fabs(coefficients[i]);
			else {
				++i;
				continue;
			}
			for (int j = 0; j < i; ++j) {
				//loop for displaying the (x-an) parts of the polynomial
				if (difference_matrix[j + 1][0] == 0)
					std::cout << "x ";
				else{
					if (difference_matrix[j + 1][0] > 0)
						std::cout << "(x - " << difference_matrix[j + 1][0] << ")";
					else
						std::cout << "(x + " << fabs(difference_matrix[j + 1][0]) << ")";
				}
			}
		}
	}
	std::cout << "\n";
}
std::vector<float> multiplyPolynomial(std::vector<float> poly_1, std::vector<float> poly_2){
	std::vector<float> product(poly_1.size() + poly_2.size() - 1, 0);
	// Multiply two polynomials term by term and take every term of first polynomial
	for (int i = 0; i < poly_1.size(); i++)	
		// Multiply the current term of first polynomial with every term of second polynomial.
		for (int j = 0; j < poly_2.size(); j++)
			product[i + j] += poly_1[i] * poly_2[j];
	return product;
}
void makeLagrange(){
	std::vector<float> x_coefficients, y_coefficients, lagrange_coefficients;
	float coeff_over_y;
	std::cout << "The Lagrange Polynomial is:\n";
	for (int i = 0; i < elements_size; ++i) {
		x_coefficients.push_back(difference_matrix[i + 1][0]);
		y_coefficients.push_back(difference_matrix[i + 1][1]);
	}
	for (int i = 0; i < elements_size; ++i) {
		float denominator = 1;
		for (int j = 0; j < elements_size; ++j)
			if (j != i)
				denominator *= x_coefficients[i] - x_coefficients[j];
		lagrange_coefficients.push_back((1/denominator) * y_coefficients[i]);
		std::wcout << lagrange_coefficients[i] << " ";
	}
	std::wcout << "\n";
	for (int i = 0; i < elements_size; ++i){
		if (i == 0) {
			std::cout << lagrange_coefficients[i];
			for (int j = 0; j < elements_size; ++j) {
				//loop for displaying the (x-an) parts of the polynomial
				if (j != i) {
					if (x_coefficients[j] == 0)
						std::cout << "x ";
					else {
						if (x_coefficients[j] > 0)
							std::cout << "(x - " << x_coefficients[j] << ")";
						else
							std::cout << "(x + " << fabs(x_coefficients[j]) << ")";
					}
				}

			}
		}
		else {
			if (lagrange_coefficients[i] > 0)
				std::cout << " + " << lagrange_coefficients[i];
			else if (lagrange_coefficients[i] < 0)
				std::cout << " - " << fabs(lagrange_coefficients[i]);
			else {
				++i;
				continue;
			}
			for (int j = 0; j < elements_size; ++j) {
				//loop for displaying the (x-an) parts of the polynomial
				if (j != i) {
					if (x_coefficients[j] == 0)
						std::cout << "x ";
					else {
						if (x_coefficients[j] > 0)
							std::cout << "(x - " << x_coefficients[j] << ")";
						else
							std::cout << "(x + " << fabs(x_coefficients[j]) << ")";
					}
				}	
			}
		}

	}
	std::cout << "\n";
}
void printSimplePoly(std::vector<float> in_vec){
	for (int i = in_vec.size(); i >= 0; --i){
		if (i > 1) {
			if (in_vec[i] > 0) {
				if (in_vec[i - 1] >= 0) 
					std::cout << in_vec[i] << "x^" << i << " + ";
				else if (in_vec[i - 1] < 0)
					std::cout << in_vec[i] << "x^" << i << " - ";
			}
			else if (in_vec[i] < 0) {
				if (in_vec[i - 1] >= 0)
					std::cout << fabs(in_vec[i]) << "x^" << i << " + ";
				else if (in_vec[i - 1] < 0) 
					std::cout << fabs(in_vec[i]) << "x^" << i << " - ";
			}
			else {
				if (in_vec[i - 1] >= 0)
					std::cout << " + ";
				else if (in_vec[i - 1] < 0) 
					std::cout << " - ";
			}
		}
		else if (i == 1) {
			if (in_vec[i] > 0) {
				if (in_vec[i - 1] >= 0)
					std::cout << in_vec[i] << "x" << " + ";
				else if (in_vec[i - 1] < 0)
					std::cout << in_vec[i] << "x" << " - ";
			}
			else if (in_vec[i] < 0) {
				if (in_vec[i - 1] >= 0)
					std::cout << fabs(in_vec[i]) << "x" << " + ";
				else if (in_vec[i - 1] < 0)
					std::cout << fabs(in_vec[i]) << "x" << " - ";
			}
			else {
				if (in_vec[i - 1] >= 0)
					std::cout << " + ";
				else if (in_vec[i - 1] < 0)
					std::cout << " - ";
			}
		}
		else {
			if (in_vec[i] > 0)
				std::cout << in_vec[i] << "\n";
			else if (in_vec[i] < 0)
				std::cout << fabs(in_vec[i]) << "\n";
			else 
				std::cout << "\n";
		}
	}
}
void makeSimplePoly(){
	std::vector<float> binomial_1(2, 1); std::vector<float> binomial_2(2, 1); //binomials to store 
	std::vector<float> matrix_numbers;
	std::vector<float> new_poly; //intermediate polynomial
	int size = 2;
	std::vector<float> final_poly; //what we print
	std::cout << "\nSIMPLIFIED POLYNOMIAL:\n\n";
	for (int i = 0; i < coefficients.size(); ++i)
		final_poly.push_back(0);
	for (int i = 1; i < elements_size - 1; ++i)//push all of the elements in first column into a std::vector to be used in multiplying the (x - an) terms
		matrix_numbers.push_back(difference_matrix[i][0]);
	final_poly[0] = coefficients[0];//pushing the first term to the first element of the array
	for (int i = 1; i < coefficients.size(); ++i){
		//multiply the x term coefficient by x term; add result to final polynomual.
		if (i == 1){
			binomial_1[0] = -matrix_numbers[0]; //c(1+x)
			for (int j = 0; j < 2; ++j)
				final_poly[j] += coefficients[1] * binomial_1[j];
		}
		else{
			//multiply the x^2 term coefficient by x^2 term; add result to final polynomual.
			if (i == 2){				
				binomial_2[0] = -matrix_numbers[i - 1];
				binomial_1 = multiplyPolynomial(binomial_1, binomial_2);
				for (int k = 0; k < binomial_1.size(); ++k)
					final_poly[k] += coefficients[i] * binomial_1[k];
			}
			else {
				for (int j = 3; j < i; ++j){
					//multiply the x^n term coefficient by x^n term; add result to final polynomual.
					binomial_2[0] = -matrix_numbers[j - 1];
					if (binomial_1.size() <= coefficients.size())
						binomial_1 = multiplyPolynomial(binomial_1, binomial_2);
					
					for (int k = 0; k < binomial_1.size(); ++k)
						final_poly[k] += coefficients[i] * binomial_1[k];
				}
			}
		}
	}
	printSimplePoly(final_poly);
}
int main() {
	menuandInputs();
	makeDifferenceTable();
	makeNewton();
	makeLagrange();
	makeSimplePoly();
}