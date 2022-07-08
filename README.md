<h1 align="center">
	Complex numbers matrix
</h1>

<p align="center">
	Simple C++ implementation of complex numbers matrix using vector of vectors from std.
</p>

<p align="center">
	<img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/haimasker/Complex-numbers-matrix?color=blue" />
	<img alt="Number of lines of code" src="https://img.shields.io/tokei/lines/github/haimasker/Complex-numbers-matrix?color=blue" />
	<img alt="Code language count" src="https://img.shields.io/github/languages/count/haimasker/Complex-numbers-matrix?color=blue" />
	<img alt="GitHub top language" src="https://img.shields.io/github/languages/top/haimasker/Complex-numbers-matrix?color=blue" />
</p>

<h3 align="center">
	<a href="#preamble">Preamble</a>
	<span> · </span>
  <a href="#fields">Class fields</a>
	<span> · </span>
	<a href="#constructors">Constructors and destructor</a>
	<span> · </span>
  <a href="#operators">Operators</a>
  <span> · </span>
  <a href="#methods">Methods</a>
</h3>

---

<a name="preamble"></a>
## Preamble

The purpose of this project is to create `Matrix` class that will be helpful in future projects. <br>
Matrix is presented via `std::vector<>` of `std::vector<>` of `std::complex<>` of `double`

---

<a name="fields"></a>
## Class fields

All of `Matrix` class fields are private <br>

1. `rows` - number of rows in matrix
```cpp
unsigned rows;
```

2. `cols` - number of columns in matrix
```cpp
unsigned cols;
```

3. `mat` - matrix of elements (`rows` x `cols` sized)
```cpp
std::vector<std::vector<std::complex<double>>> mat;
```

---

<a name="constructors"></a>
## Constructors and destructor

1. The Constructor receiving the number of `rows`, `cols` and initial value for elements. <br>
By default have values `0`, `0` and `0+0i` respectively.
```cpp
Matrix(unsigned = 0, unsigned = 0, std::complex<double> = 0.0 + 0.0i);
```

<br>

2. The constructor receiving a vector that represents a matrix. <br>
Vector should be square.
```cpp
Matrix(std::vector<std::vector<std::complex<double>>>*);
```

<br>

3. The constructor receiving the number of `rows`, `cols`, vector that represents a matrix and initial value for elements. <br>
Makes `rows` x `cols` sized matrix, fills it elementwise from given vector. <br>
If there are not enough elements in the given vector, then fills remaining elements with initial value.
```cpp
Matrix(unsigned, unsigned, std::vector<std::vector<std::complex<double>>>*, std::complex<double> = 0.0 + 0.0i);
```

<br>

4. The constructor receiving the number of `rows`, `cols`, vector of elements and initial value for elements. <br>
Makes `rows` x `cols` sized matrix, fills it elementwise from given vector. <br>
If there are not enough elements in the given vector, then fills remaining elements with initial value.
```cpp
Matrix(unsigned, unsigned, std::vector<std::complex<double>>*, std::complex<double> = 0.0 + 0.0i);
```

<br>

5. Copy constructor
```cpp
Matrix(const Matrix&);
```

<br>

6. Destructor
```cpp
virtual ~Matrix();
```

---

<a name="operators"></a>
## Operators

1. Call operator for non-const objects receiving `row` and `col`. <br>
Returns `mat[row][col]`. <br>
Indices are cycled, so there is no out-of-reange error
```cpp
std::complex<double>& operator () (const unsigned&, const unsigned&);
```

<br>

2. Call operator for const objects receiving `row` and `col`. <br>
Returns `mat[row][col]`. <br>
Indices are cycled, so there is no out-of-reange error
```cpp
const std::complex<double>& operator () (const unsigned&, const unsigned&) const;
```

<br>

3. Copy assignment operator.
```cpp
Matrix& operator = (const Matrix&);
```

<br>

4. Addition operator. <br>
Sums matrices elementwise.
```cpp
Matrix operator + (const Matrix&);
```

<br>

5. Addition assignment operator. <br>
Sums matrices elementwise.
```cpp
Matrix& operator += (const Matrix&);
```

<br>

6. Subtraction operator. <br>
Subtracts matrices elementwise.
```cpp
Matrix operator - (const Matrix&);
```

<br>

7. Subtraction assignment operator. <br>
Subtracts matrices elementwise.
```cpp
Matrix& operator -= (const Matrix&);
```

<br>

8. Multiplication operator. <br>
Multiplies matrices elementwise.
```cpp
Matrix operator * (const Matrix&);
```

<br>

9. Multiplication assignment operator. <br>
Multiplies matrices elementwise.
```cpp
Matrix& operator *= (const Matrix&);
```

<br>

10. Assignment operator. <br>
Assigns matrix elements to given value.
```cpp
Matrix& operator = (const std::complex<double>&);
```

<br>


11. Addition operator. <br>
Sums matrix elements with given value.
```cpp
Matrix operator + (const std::complex<double>&);
```

<br>

12. Addition assignment operator. <br>
Sums matrix elements with given value.
```cpp
Matrix& operator += (const std::complex<double>&);
```

<br>

13. Subtraction operator. <br>
Subtracts given value from matrix elements.
```cpp
Matrix operator - (const std::complex<double>&);
```

<br>

14. Subtraction assignment operator. <br>
Subtracts given value from matrix elements.
```cpp
Matrix& operator -= (const std::complex<double>&);
```

<br>

15. Multiplication operator. <br>
Multiplies matrix elements by given value.
```cpp
Matrix operator * (const std::complex<double>&);
```

<br>

16. Multiplication assignment operator. <br>
Multiplies matrix elements by given value.
```cpp
Matrix& operator *= (const std::complex<double>&);
```

<br>

17. Division operator. <br>
Divides matrix elements by given value.
```cpp
Matrix operator / (const std::complex<double>&);
```

<br>

18. Division assignment operator. <br>
Divides matrix elements by given value.
```cpp
Matrix& operator /= (const std::complex<double>&);
```

<br>

19. Equal operator. <br>
Determines if matrices are equal.
```cpp
bool operator == (const Matrix&);
```

<br>

20. Non-equal operator. <br>
Determines if matrices are not-equal.
```cpp
bool operator != (const Matrix&);
```

<br>

21. Output operator. <br>
Prints matrix to output.
```cpp
friend std::ostream& operator << (std::ostream&, const Matrix&);
```

---

<a name="methods"></a>
## Methods

1. `elementwiseProduct` <br>
Multiplies `this->mat` by given matrix elementwise.
```cpp
void elementwiseProduct(const Matrix&);
```

<br>

2. `elementwiseProduct` <br>
Returns elementwise product of given matrices.
```cpp
static Matrix elementwiseProduct(const Matrix&, const Matrix&);
```

<br>

3. `elementwiseDivision` <br>
Divides `this->mat` by given matrix elementwise.
```cpp
void elementwiseDivision(const Matrix&);
```

<br>

4. `elementwiseDivision` <br>
Return elemetwise division of given matrices.
```cpp
static Matrix elementwiseDivision(const Matrix&, const Matrix&);
```

<br>

5. `conjugate` <br>
Transforms elements of `this->mat` to it's conjugate.
```cpp
void conjugate();
```

<br>

6.  `conjugate` <br>
Returns conjugate of given matrix.
```cpp
static Matrix conjugate(const Matrix&);
```

<br>

7. `transpose` <br>
Transposes `this->mat`.
```cpp
void transpose();
```

<br>

8. `transpose` <br>
Returns transposition of given matrix.
```cpp
static Matrix transpose(const Matrix&);
```

<br>

9. `hermitianConjugate` <br>
Transform `this->mat` to it's Hermitian-conjugate.
```cpp
void hermitianConjugate();
```

<br>

10. `hermitianConjugate` <br>
Returns Hermitian-conjugate of given matrix.
```cpp
static Matrix hermitianConjugate(const Matrix&);
``` 

<br>

11. `tensorProduct` <br>
Produces a tensor product with the given matrix.
```cpp
void tensorProduct(const Matrix&);
```

<br>

12. `tensorProduct` <br>
Returns tensor product of given matrices.
```cpp
static Matrix tensorProduct(const Matrix&, const Matrix&);
```

<br>

13. `cofactor` <br>
Returns `false` if matrix is not square. <br>
Returns `true` and transforms `this->mat` to it's co-factor matrix if matrix is square.
```cpp
bool cofactor();
```

<br>

14. `adjugate` <br>
Returns `false` if matrix is not square. <br>
Returns `true` and transforms `this->mat` to it's adjugate matrix if matrix is square.
```cpp
bool adjugate();
```

<br>

15. `inverse` <br>
Returns `false` if matrix is not square. <br>
Returns `true` and transforms `this->mat` to it's inverse matrix if matrix is square.
```cpp
bool inverse();
```

<br>

16. `resize` <br>
Receive `row`, `col` and initial value `init`. <br>
Resizes matrix to `row` x `col`. <br>
If matrix becomes greater (in any dimension), then new empty elements are filled with `init` value.
```cpp
void resize(const unsigned&, const unsigned&, const std::complex<double>& = 0.0 + 0.0i);
```

<br>

17. `rehape` <br>
Returns `false` if reshaping is impossible. <br>
Returns `true` and reshapes `this->mat` if it is possible.
```cpp
bool reshape(const unsigned&, const unsigned&);
```

<br>

18. `trace` <br>
Returns trace of matrix.
```cpp
std::complex<double> trace();
```

<br>

19. `determinant` <br>
Returns `0+0i` if matrix is not square. <br>
Returns determinant of matrix if it is square.
```cpp
std::complex<double> determinant();
```

<br>

20. `identity` <br>
Returns identity matrix with given size.
```cpp
static Matrix identity(const unsigned& = 1);
```

<br>

21. `getRows` <br>
Returns `rows` value.
```cpp
unsigned getRows();
```

<br>

22. `getCols` <br>
Returns `cols` value.
```cpp
unsigned getCols();
```
