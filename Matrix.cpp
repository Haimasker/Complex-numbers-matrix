#include "Matrix.h"

Matrix::Matrix(unsigned row, unsigned col, std::complex<double> init) noexcept {
	mat.resize(row);
	for (unsigned i = 0; i < row; i++)
		mat[i] = std::vector<std::complex<double>>(col, init);
	rows = row;
	cols = (row == 0) ? 0 : col;
}

Matrix::Matrix(std::vector<std::vector<std::complex<double>>>* m) noexcept {
	bool square = true;
	for (unsigned i = 1; i < (*m).size() && square; i++)
		square = (*m)[i].size() == (*m)[0].size();
	if (!square) {
		this->cols = this->rows = 0;
		this->mat.resize(0);
		return;
	}
	this->mat = *m;
	this->rows = (unsigned)m->size();
	this->cols = (unsigned)(*m)[0].size();
}

Matrix::Matrix(unsigned row, unsigned col, std::vector<std::vector<std::complex<double>>>* m, std::complex<double> init) noexcept {
	std::vector<std::vector<std::complex<double>>> res(row);
	for (size_t i = 0; i < res.size(); i++)
		res[i] = std::vector<std::complex<double>>(col);

	unsigned mSize = 0;
	for (unsigned i = 0; i < (*m).size(); i++)
		mSize += (unsigned)(*m)[i].size();

	for (unsigned i = 0; i < std::min(row * col, mSize); i++)
		res[i / col][i % col] = (*m)[i / (*m)[0].size()][i % (*m)[0].size()];
	if (mSize < (row * col))
		for (; mSize < row * col; mSize++)
			res[mSize / col][mSize % col] = init;

	this->mat = res;
	this->rows = row;
	this->cols = (row == 0) ? 0 : col;
}

Matrix::Matrix(unsigned row, unsigned col, std::vector<std::complex<double>>* m, std::complex<double> init) noexcept {
	std::vector<std::vector<std::complex<double>>> res(row);
	for (unsigned i = 0; i < res.size(); i++)
		res[i] = std::vector<std::complex<double>>(col);

	for (unsigned i = 0; i < std::min((size_t)row * col, (*m).size()); i++)
		res[i / col][i % col] = (*m)[i];
	if ((*m).size() < (size_t)row * col)
		for (unsigned i = (unsigned)(*m).size(); i < row * col; i++)
			res[i / col][i % col] = init;

	this->mat = res;
	this->rows = row;
	this->cols = (row == 0) ? 0 : col;
}

Matrix::Matrix(const Matrix& m) noexcept {
	this->mat = m.mat;
	this->rows = m.rows;
	this->cols = m.cols;
}

Matrix::Matrix(const Matrix&& m) noexcept :
	rows(std::exchange(m.rows, 0)),
	cols(std::exchange(m.cols, 0)),
	mat(std::move(m.mat)) {
}

Matrix::~Matrix() {}

std::complex<double>& Matrix::operator () (const unsigned& row, const unsigned& col) {
	return this->mat[row % this->rows][col % this->cols];
}

const std::complex<double>& Matrix::operator () (const unsigned& row, const unsigned& col) const {
	return this->mat[row % this->rows][col % this->cols];
}

Matrix& Matrix::operator = (const Matrix& m) {
	if (this == &m || this->mat == m.mat)
		return *this;
	this->mat = m.mat;
	this->rows = m.rows;
	this->cols = m.cols;
	return *this;
}

Matrix Matrix::operator + (const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return Matrix();
	Matrix res(this->rows, this->cols, 0.0 + 0.0i);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			res(i, j) = this->mat[i][j] + m(i, j);
	return res;
}

Matrix& Matrix::operator += (const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return *this;
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] += m(i, j);
	return *this;
}

Matrix Matrix::operator - (const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return Matrix();
	Matrix res(this->rows, this->cols, 0.0 + 0.0i);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			res(i, j) = this->mat[i][j] - m(i, j);
	return res;
}

Matrix& Matrix::operator -= (const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return *this;
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] -= m(i, j);
	return *this;
}

Matrix Matrix::operator * (const Matrix& m) {
	if (this->cols != m.rows)
		return Matrix();
	Matrix res(this->rows, m.cols, 0.0 + 0.0i);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < m.cols; j++)
			for (unsigned k = 0; k < m.rows; k++)
				res(i, j) += this->mat[i][k] * m(k, j);
	return res;
}

Matrix& Matrix::operator *= (const Matrix& m) {
	Matrix res = (*this) * m;
	(*this) = res;
	return *this;
}

Matrix& Matrix::operator = (const std::complex<double>& num) {
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] = num;
	return *this;
}

Matrix Matrix::operator + (const std::complex<double>& num) {
	Matrix res(this->rows, this->cols);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			res(i, j) = this->mat[i][j] + num;
	return res;
}

Matrix& Matrix::operator += (const std::complex<double>& num) {
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] += num;
	return *this;
}

Matrix Matrix::operator - (const std::complex<double>& num) {
	Matrix res(this->rows, this->cols);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			res(i, j) = this->mat[i][j] - num;
	return res;
}

Matrix& Matrix::operator -= (const std::complex<double>& num) {
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] -= num;
	return *this;
}

Matrix Matrix::operator * (const std::complex<double>& num) {
	Matrix res(this->rows, this->cols);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			res(i, j) = this->mat[i][j] * num;
	return res;
}

Matrix& Matrix::operator *= (const std::complex<double>& num) {
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] *= num;
	return *this;
}

Matrix Matrix::operator / (const std::complex<double>& num) {
	Matrix res(this->rows, this->cols);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			res(i, j) = this->mat[i][j] / num;
	return res;
}

Matrix& Matrix::operator /= (const std::complex<double>& num) {
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] /= num;
	return *this;
}

bool Matrix::operator == (const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return false;
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			if (this->mat[i][j] != m(i, j))
				return false;
	return true;
}

bool Matrix::operator != (const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return true;
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			if (this->mat[i][j] != m(i, j))
				return true;
	return false;
}

void Matrix::elementwiseAddition(const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return;
	for (unsigned i = 0; i < this->rows * this->cols; i++)
		this->mat[i / this->cols][i % this->cols] += m(i / m.cols, i % m.cols);
	return;
}

Matrix Matrix::elementwiseAddition(const Matrix& m1, const Matrix& m2) {
	if (m1.rows != m2.rows || m1.cols != m2.cols)
		return Matrix();
	Matrix res(m1);
	for (unsigned i = 0; i < res.rows * res.cols; i++)
		res(i / res.cols, i % res.cols) += m2(i / m2.cols, i % m2.cols);
	return res;
}

void Matrix::elementwiseSubtraction(const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return;
	for (unsigned i = 0; i < this->rows * this->cols; i++)
		this->mat[i / this->cols][i % this->cols] -= m(i / m.cols, i % m.cols);
	return;
}

Matrix Matrix::elementwiseSubtraction(const Matrix& m1, const Matrix& m2) {
	if (m1.rows != m2.rows || m1.cols != m2.cols)
		return Matrix();
	Matrix res(m1);
	for (unsigned i = 0; i < res.rows * res.cols; i++)
		res(i / res.cols, i % res.cols) -= m2(i / m2.cols, i % m2.cols);
	return res;
}

void Matrix::elementwiseProduct(const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return;
	for (unsigned i = 0; i < this->rows * this->cols; i++)
		this->mat[i / this->cols][i % this->cols] *= m(i / m.cols, i % m.cols);
	return;
}

Matrix Matrix::elementwiseProduct(const Matrix& m1, const Matrix& m2) {
	if (m1.rows != m2.rows || m1.cols != m2.cols)
		return Matrix();
	Matrix res(m1);
	for (unsigned i = 0; i < res.rows * res.cols; i++)
		res(i / res.cols, i % res.cols) *= m2(i / m2.cols, i % m2.cols);
	return res;
}

void Matrix::elementwiseDivision(const Matrix& m) {
	if (this->rows != m.rows || this->cols != m.cols)
		return;
	for (unsigned i = 0; i < this->rows * this->cols; i++)
		this->mat[i / this->cols][i % this->cols] /= m(i / m.cols, i % m.cols);
	return;
}

Matrix Matrix::elementwiseDivision(const Matrix& m1, const Matrix& m2) {
	if (m1.rows != m2.rows || m1.cols != m2.cols)
		return Matrix();
	Matrix res(m1);
	for (unsigned i = 0; i < res.rows * res.cols; i++)
		res(i / res.cols, i % res.cols) /= m2(i / m2.cols, i % m2.cols);
	return res;
}

void Matrix::conjugate() {
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			this->mat[i][j] = std::conj(this->mat[i][j]);
	return;
}

Matrix Matrix::conjugate(const Matrix& m) {
	Matrix res(m);
	for (unsigned i = 0; i < res.rows; i++)
		for (unsigned j = 0; j < res.cols; j++)
			res.mat[i][j] = std::conj(res.mat[i][j]);
	return res;
}

void Matrix::transpose() {
	std::vector<std::vector<std::complex<double>>> m(this->cols);
	for (unsigned i = 0; i < this->cols; i++)
		m[i] = std::vector<std::complex<double>>(this->rows);
	for (unsigned i = 0; i < this->rows; i++)
		for (unsigned j = 0; j < this->cols; j++)
			m[j][i] = this->mat[i][j];
	this->mat = m;
	this->rows ^= this->cols ^= this->rows ^= this->cols;
	return;
}

Matrix Matrix::transpose(const Matrix& m) {
	std::vector<std::vector<std::complex<double>>> res(m.cols);
	for (unsigned i = 0; i < m.cols; i++)
		res[i] = std::vector<std::complex<double>>(m.rows);
	for (unsigned i = 0; i < m.rows; i++)
		for (unsigned j = 0; j < m.cols; j++)
			res[j][i] = m.mat[i][j];
	return Matrix(&res);
}

void Matrix::hermitianConjugate() {
	this->transpose();
	this->conjugate();
	return;
}

Matrix Matrix::hermitianConjugate(const Matrix& m) {
	Matrix res(m);
	res.transpose();
	res.conjugate();
	return res;
}

void Matrix::tensorProduct(const Matrix& m) {
	std::vector<std::vector<std::complex<double>>> res;
	res.resize((size_t)this->rows * m.rows);
	for (unsigned i = 0; i < res.size(); i++)
		res[i] = std::vector<std::complex<double>>(this->cols * m.cols);
	for (unsigned i = 0; i < this->rows * m.rows; i++)
		for (unsigned j = 0; j < this->cols * m.cols; j++)
			res[i][j] = this->mat[i / m.rows][j / m.cols] * m(i % m.rows, j % m.cols);
	this->mat = res;
	this->rows *= m.rows;
	this->cols *= m.cols;
	return;
}

Matrix Matrix::tensorProduct(const Matrix& m1, const Matrix& m2) {
	std::vector<std::vector<std::complex<double>>> res;
	res.resize((size_t)m1.rows * m2.rows);
	for (unsigned i = 0; i < res.size(); i++)
		res[i] = std::vector<std::complex<double>>(m1.cols * m2.cols);
	for (unsigned i = 0; i < m1.rows * m2.rows; i++)
		for (unsigned j = 0; j < m1.cols * m2.cols; j++)
			res[i][j] = m1(i / m2.rows, j / m2.cols) * m2(i % m2.rows, j % m2.cols);
	return Matrix(&res);
}

bool Matrix::cofactor() {
	if (this->rows != this->cols)
		return false;
	std::vector<std::vector<std::complex<double>>> res(this->rows);
	for (unsigned i = 0; i < res.size(); i++)
		res[i] = std::vector<std::complex<double>>(this->cols);
	for (unsigned i = 0; i < this->rows * this->cols; i++) {
		Matrix m(*this);
		m.mat.erase(m.mat.begin() + i % this->cols);
		for (unsigned j = 0; j < m.mat.size(); j++)
			m.mat[j].erase(m.mat[j].begin() + i % this->cols);
		res[i / this->cols][i % this->cols] = m.determinant();
	}
	this->mat = res;
	return true;
}

bool Matrix::adjugate() {
	if (!(this->cofactor()))
		return false;
	this->transpose();
	return true;
}

bool Matrix::inverse() {
	if (this->determinant() == 0.0 + 0.0i)
		return false;
	this->adjugate();
	*this /= this->determinant();
	return true;
}

void Matrix::resize(const unsigned& row, const unsigned& col, const std::complex<double>& init) {
	if (this->rows == row && this->cols == col)
		return;
	if (row < this->rows)
		this->mat.resize(row);
	else if (row > this->rows)
		this->mat.resize(row, std::vector<std::complex<double>>(this->cols, init));
	if (col < this->cols)
		for (unsigned i = 0; i < this->mat.size(); i++)
			this->mat[i].resize(col);
	else if (col > this->cols)
		for (unsigned i = 0; i < this->mat.size(); i++)
			this->mat[i].resize(col, init);
	this->rows = row;
	this->cols = col;
	return;
}

bool Matrix::reshape(const unsigned& row, const unsigned& col) {
	if ((this->rows * this->cols != row * col) || (this->rows == row && this->cols == col))
		return false;
	std::vector<std::vector<std::complex<double>>> m(row);
	for (unsigned i = 0; i < m.size(); i++)
		m[i] = std::vector<std::complex<double>>(col);
	for (unsigned i = 0; i < row * col; i++)
		m[i / col][i % col] = this->mat[i / this->cols][i % this->cols];
	this->mat = m;
	this->rows = row;
	this->cols = col;
	return true;
}

std::complex<double> Matrix::trace() {
	std::complex<double> res = 0;
	for (unsigned i = 0; i < std::min(this->rows, this->cols); i++)
		res += this->mat[i][i];
	return res;
}

std::complex<double> Matrix::determinant() {
	std::complex<double> res = 0.0 + 0.0i;
	if (this->rows != this->cols)
		return res;
	if (this->rows == 1)
		return this->mat[0][0];
	if (this->rows == 2)
		return this->mat[0][0] * this->mat[1][1] - this->mat[0][1] * this->mat[1][0];

	std::complex<double> degree = 1.0 + 0.0i;
	Matrix m(*this);

	for (unsigned i = 0; i < this->rows; i++) {
		Matrix tmp(m);
		tmp.mat.erase(tmp.mat.begin());
		for (unsigned j = 0; j < tmp.mat.size(); j++)
			tmp.mat[j].erase(tmp.mat[j].begin() + i);
		tmp.rows--;
		tmp.cols--;
		res += degree * m.mat[0][i] * tmp.determinant();
		degree *= -1;
	}
	return res;
}

Matrix Matrix::identity(const unsigned& size) {
	Matrix res(size, size);
	for (unsigned i = 0; i < size; i++)
		res(i, i) = 1.0 + 0.0i;
	return res;
}

std::vector<std::vector<std::complex<double>>> Matrix::getMat() {
	return this->mat;
}

unsigned Matrix::getRows() {
	return this->rows;
}

unsigned Matrix::getCols() {
	return this->cols;
}

std::ostream& operator << (std::ostream& out, Matrix& m) noexcept {
	if (&m != nullptr && m.rows)
		for (std::vector<std::complex<double>> i : m.mat) {
			for (std::complex<double> j : i)
				out << std::setw(5) << j << " ";
			out << "\n";
		}
	else
		out << "Empty matrix\n";
	out << "\n";
	return out;
}
