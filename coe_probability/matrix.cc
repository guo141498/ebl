#include "matrix.h"
#include <cstdlib>

using namespace std;

bool matrix::check_data() const
{
  if (data.empty()) return false;

  int colsize = data[0].size();
  for (const auto& d : data)
    if (d.size() != colsize) return false;

  return true;
}

matrix::matrix() {}

matrix::matrix(const vector<vector<complex<double> > >& data_) : data(data_)
{
  if(!check_data()) {
    cout << "matrix::ERROR::The input data is ill, exit!" << endl;
    exit(1);
  }

  nrow = data.size();
  ncol = data[0].size();
}

matrix::matrix(unsigned nrow_, unsigned ncol_, const vector<complex<double> >& data_)
  : nrow(nrow_), ncol(ncol_)
{
  if (data_.size() < nrow * ncol) {
    cout << "matrix::ERROR::The input data is ill, exit!" << endl;
    exit(1);
  }

  unsigned iter = 0;

  data.resize(nrow);
  for (unsigned irow = 0; irow < nrow; irow++) {
    data[irow].resize(ncol);
    for (unsigned icol = 0; icol < ncol; icol++, iter++)
      data[irow][icol] = data_[iter];
  }
}

matrix& matrix::operator*=(const matrix& rhs) {
  if (ncol != rhs.nrow) {
    cout << "matrix::*=::ERROR::You are trying to do producting with two unmatch matrix, please check" << endl;
    cout << "M1 " << (*this) << "M2 " << rhs;
    exit(0);
  }

  vector<vector<complex<double> > > result;
  result.resize(nrow);
  for (unsigned irow = 0; irow < nrow; irow++) {
    result[irow].resize(rhs.ncol, polar(0.0, 0.0));
    for (unsigned icol = 0; icol < rhs.ncol; icol++)
      for (unsigned i = 0; i < ncol; i++)
        result[irow][icol] += data[irow][i] * rhs[i][icol];
  }
  data = result;

  return (*this);
}

matrix matrix::operator*(const matrix& rhs) const
{
  matrix result;
  result = (*this);
  result *= rhs;
  return result;
}

matrix& matrix::operator+=(const matrix& rhs) {
  if (ncol != rhs.ncol || nrow != rhs.nrow) {
    cout << "matrix::+=::ERROR::You are trying to do + with two unmatch matrix, please check" << endl;
    cout << "M1 " << (*this) << "M2 " << rhs;
    exit(0);
  }

  for (unsigned irow = 0; irow < nrow; irow++)
    for (unsigned icol = 0; icol < ncol; icol++)
        data[irow][icol] += rhs[irow][icol];

  return (*this);
}

matrix matrix::operator+(const matrix& rhs) const
{
  matrix result;
  result = (*this);
  result += rhs;
  return result;
}

ostream& operator<<(std::ostream& os, const matrix& m) {
  for (const auto& v : m.data) {
    os << endl << "| ";
    for (const auto& c : v)
      os << c << " ";
    os << "|";
  }

  os << endl;
  return os;
}
   

matrix matrix::conjugate() const {
	vector<vector<complex<double> > > conj_data;
	conj_data.resize(ncol);

	for (auto& row : conj_data) row.resize(nrow);
	for (unsigned irow = 0; irow < nrow; irow++)
		for (unsigned icol = 0; icol < ncol; icol++)
			conj_data[icol][irow] = conj(data[irow][icol]);

    matrix result(conj_data);
    return result;
}

//int main()
//{
//    double pi = 3.1415926;
//  vector<complex<double> > m1 = {
//    polar(1.0, 0.0),
//    polar(1.0, pi/6),
//    0,
//
//    0,
//    polar(1.0, 3.14),
//    0,
//
//    0,
//    0,
//    polar(1.0, 3.14)
//  };
//
//  matrix M1(3, 3, m1);
//  cout << M1;
//
//  vector<complex<double> > m2 = {
//    polar(1.0, 0.0),
//    0,
//    0,
//
//    0,
//    polar(1.0, 3.1415),
//    0,
//
//    0,
//    0,
//    polar(1.0, 0.0)
//  };
//
//  matrix M2(3, 3, m2);
//  cout <<"M1 * M2= "<< M1 * M2;
//  cout <<"M2 * M1= "<<  M2 + M1;
//  
//  vector<complex<double> > m2 = { polar(1.0, 0.0), 0,
//    polar(0.5, 1.68), polar(1.0, 3.14) };
//
//
//  matrix M2(2, 2, m2);
//  cout << M1 * M2;
//  cout << M2;
//
//  return 0;
//}
