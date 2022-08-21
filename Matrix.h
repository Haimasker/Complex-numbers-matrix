#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>

using namespace std::complex_literals;

class Matrix {
private:
	std::vector<std::vector<std::complex<double>>> mat;
	unsigned rows;
	unsigned cols;

public:
	Matrix(unsigned = 0, unsigned = 0, std::complex<double> = 0.0 + 0.0i) noexcept;
	Matrix(std::vector<std::vector<std::complex<double>>>*) noexcept;
	Matrix(unsigned, unsigned, std::vector<std::vector<std::complex<double>>>*, std::complex<double> = 0.0 + 0.0i) noexcept;
	Matrix(unsigned, unsigned, std::vector<std::complex<double>>*, std::complex<double> = 0.0 + 0.0i) noexcept;
	Matrix(const Matrix&) noexcept;
	Matrix&	operator = (const Matrix&) noexcept;
	Matrix(const Matrix&&) noexcept;
	Matrix&	operator = (const Matrix&&) noexcept;

	virtual	~Matrix();

	std::complex<double>&	operator () (const unsigned&, const unsigned&);
	const std::complex<double>&	operator () (const unsigned&, const unsigned&) const;

	Matrix&	operator = (const Matrix&);
	Matrix	operator + (const Matrix&);
	Matrix&	operator += (const Matrix&);
	Matrix	operator - (const Matrix&);
	Matrix&	operator -= (const Matrix&);
	Matrix	operator * (const Matrix&);
	Matrix&	operator *= (const Matrix&);

	Matrix&	operator = (const std::complex<double>&);
	Matrix	operator + (const std::complex<double>&);
	Matrix&	operator += (const std::complex<double>&);
	Matrix	operator - (const std::complex<double>&);
	Matrix&	operator -= (const std::complex<double>&);
	Matrix	operator * (const std::complex<double>&);
	Matrix&	operator *= (const std::complex<double>&);
	Matrix	operator / (const std::complex<double>&);
	Matrix&	operator /= (const std::complex<double>&);

	bool	operator == (const Matrix&);
	bool	operator != (const Matrix&);

	void		elementwiseAddition(const Matrix&);
	static Matrix	elementwiseAddition(const Matrix&, const Matrix&);
	void		elementwiseSubtraction(const Matrix&);
	static Matrix	elementwiseSubtraction(const Matrix&, const Matrix&);
	void		elementwiseProduct(const Matrix&);
	static Matrix	elementwiseProduct(const Matrix&, const Matrix&);
	void		elementwiseDivision(const Matrix&);
	static Matrix	elementwiseDivision(const Matrix&, const Matrix&);

	void		conjugate();
	static Matrix	conjugate(const Matrix&);

	void		transpose();
	static Matrix	transpose(const Matrix&);

	void		hermitianConjugate();
	static Matrix	hermitianConjugate(const Matrix&);

	void		tensorProduct(const Matrix&);
	static Matrix	tensorProduct(const Matrix&, const Matrix&);

	bool	cofactor();
	bool	adjugate();
	bool	inverse();

	void	resize(const unsigned&, const unsigned&, const std::complex<double> & = 0.0 + 0.0i);
	bool	reshape(const unsigned&, const unsigned&);

	std::complex<double>	trace();
	std::complex<double>	determinant();

	static Matrix	identity(const unsigned& = 1);

	std::vector<std::vector<std::complex<double>>>	getMat();

	unsigned	getRows();
	unsigned	getCols();

	friend std::ostream&	operator << (std::ostream&, Matrix&) noexcept;
};
