#include <iostream>
#include <complex>
#include <vector>

class matrix {
  private:
    bool check_data() const;
    std::vector<std::vector<std::complex<double> > > data;
    unsigned nrow, ncol;

  public:
    matrix();
    matrix(const std::vector<std::vector<std::complex<double> > >& data_);
    matrix(unsigned nrow_, unsigned ncol_, const std::vector<std::complex<double> >& data_);

	matrix conjugate() const;

    matrix& operator*=(const matrix& rhs);
    matrix operator*(const matrix& rhs) const;

    matrix& operator+=(const matrix& rhs);
    matrix operator+(const matrix& rhs) const;
    std::vector<std::complex<double> >& operator[](unsigned i) { return data[i]; }
    const std::vector<std::complex<double> >& operator[](unsigned i) const { return data[i]; }

    friend std::ostream& operator<<(std::ostream& os, const matrix& m);
};

	

std::ostream& operator<<(std::ostream& os, const matrix& m);
