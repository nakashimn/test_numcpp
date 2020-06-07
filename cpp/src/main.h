#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <NumCpp.hpp>

using RMatrixXd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;


class bin{

public:
    bin();
    ~bin();
    bool read_bin(const std::string& filepath,
        std::vector<std::vector<double>>& vecvec);
    bool write_bin(const std::string& filepath,
        const std::vector<std::vector<double>>& vecvec);
};



class eigen {
public:
    eigen();
    ~eigen();
    bool vecvec_to_matrix(const std::vector<std::vector<double>>& vecvec,
	      RMatrixXd& matrix);
    bool eigen_dot(const RMatrixXd& src_l,
	      const RMatrixXd& src_r, RMatrixXd& dest);
    bool matrix_to_vecvec(const RMatrixXd& matrix,
        std::vector<std::vector<double>>& vecvec);
};

class blas {
public:
    blas();
    ~blas();
};

class numcpp {
public:
    numcpp();
    ~numcpp();
	bool vecvec_to_ndarray(const std::vector<std::vector<double>>& vecvec,
		nc::NdArray<double>& ndarray);
	bool ndarray_dot(const nc::NdArray<double>& src_l,
		const nc::NdArray<double>& src_r, nc::NdArray<double>& dest);
	bool ndarray_to_vecvec(const nc::NdArray<double>& ndarray,
		std::vector<std::vector<double>>& vecvec);
};
