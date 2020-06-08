#include "main.h"

int main(int, char**) {
    std::cout << "Test_NumCpp\n";
    bin bin;
    eigen eigen;
    blas blas;
    numcpp numcpp;

    std::vector<double> tmp_vec(10, 0.0);
    std::vector<std::vector<double>> input_x(10, tmp_vec);
    std::vector<std::vector<double>> input_y(10, tmp_vec);

    bin.read_bin("../../data/input_x.bin", input_x);
    bin.read_bin("../../data/input_y.bin", input_y);

    // test Eigen
    RMatrixXd eigen_input_x;
    RMatrixXd eigen_input_y;
    RMatrixXd eigen_output_z;

    eigen.vecvec_to_matrix(input_x, eigen_input_x);
    eigen.vecvec_to_matrix(input_y, eigen_input_y);
    eigen.eigen_dot(eigen_input_x, eigen_input_y, eigen_output_z);

    std::vector<std::vector<double>> eigen_to_vecvec_output_z;

    eigen.matrix_to_vecvec(eigen_output_z, eigen_to_vecvec_output_z);
    bin.write_bin("../../data/eigen_output_z.bin", eigen_to_vecvec_output_z);

    // test NumCpp
    nc::NdArray<double> numcpp_input_x;
    nc::NdArray<double> numcpp_input_y;
    nc::NdArray<double> numcpp_output_z;

    numcpp.vecvec_to_ndarray(input_x, numcpp_input_x);
    numcpp.vecvec_to_ndarray(input_y, numcpp_input_y);
    numcpp.ndarray_dot(numcpp_input_x, numcpp_input_y, numcpp_output_z);

    std::vector<std::vector<double>> numcpp_to_vecvec_output_z;

    numcpp.ndarray_to_vecvec(numcpp_output_z, numcpp_to_vecvec_output_z);
    bin.write_bin("../../data/numcpp_output_z.bin", numcpp_to_vecvec_output_z);

}


bin::bin(){}
bin::~bin(){}

bool bin::read_bin(
    const std::string& filepath,
    std::vector<std::vector<double>>& vecvec
){
    std::ifstream fin(filepath, std::ios::binary);
    for (unsigned int i = 0; i < vecvec.size(); i++) {
        fin.read((char *)&vecvec[i][0], sizeof(double) * vecvec[i].size());
    }
    return true;
}

bool bin::write_bin(
    const std::string& filepath,
    const std::vector<std::vector<double>>& vecvec
){
    std::ofstream fout(filepath, std::ios::binary);
    for (unsigned int i = 0; i < vecvec.size(); i++) {
        fout.write((char *)&vecvec[i][0], sizeof(double) * vecvec[i].size());
    }
    return true;
}


eigen::eigen() {}
eigen::~eigen() {}

bool eigen::vecvec_to_matrix(
    const std::vector<std::vector<double>>& vecvec,
    RMatrixXd& matrix
) {
    RMatrixXd tmp_matrix(vecvec.size(), vecvec[0].size());
    Eigen::VectorXd tmp_vector;
    for (int i = 0; i < tmp_matrix.rows(); i++) {
        tmp_matrix.row(i)= Eigen::VectorXd::Map(
            vecvec[i].data(), vecvec[i].size());
    }
    matrix = tmp_matrix;
    return true;
}

bool eigen::eigen_dot(
    const RMatrixXd& src_l,
    const RMatrixXd& src_r,
    RMatrixXd& dest
) {
    dest = src_l * src_r;
    return true;
}

bool eigen::matrix_to_vecvec(
    const RMatrixXd& matrix,
    std::vector<std::vector<double>>& vecvec
){
    std::vector<double> tmp_vec(matrix.cols(), 0.0);
    std::vector<std::vector<double>> tmp_vecvec(matrix.rows(), tmp_vec);
    for(unsigned int i = 0; i < tmp_vecvec.size(); i++) {
        tmp_vecvec[i] = std::vector<double>(matrix.row(i).data(),
            matrix.row(i).data()+matrix.row(i).cols());
    }
    vecvec = tmp_vecvec;
    return true;
}


blas::blas() {}
blas::~blas() {}


numcpp::numcpp() {}
numcpp::~numcpp() {}

bool numcpp::vecvec_to_ndarray(
    const std::vector<std::vector<double>>& vecvec,
    nc::NdArray<double>& ndarray
) {
    nc::NdArray<double> tmp_ndarray(vecvec.size(), vecvec[0].size());
    for (unsigned int i = 0; i < tmp_ndarray.numRows(); i++) {
        for(unsigned int j = 0; j < tmp_ndarray.numCols(); j++) {
            tmp_ndarray(i, j) = vecvec[i][j];
        }
    }
    ndarray = tmp_ndarray;
    return true;
}

bool numcpp::ndarray_dot(
    const nc::NdArray<double>& src_l,
    const nc::NdArray<double>& src_r,
    nc::NdArray<double>& dest
) {
    dest = nc::dot(src_l, src_r);
    return true;
}

bool numcpp::ndarray_to_vecvec(
    nc::NdArray<double>& ndarray,
    std::vector<std::vector<double>>& vecvec
) {
    std::vector<double> tmp_vec(ndarray.numCols(), 0.0);
    std::vector<std::vector<double>> tmp_vecvec(ndarray.numRows(), tmp_vec);
    for (unsigned int i = 0; i < tmp_vecvec.size(); i++) {
        for (unsigned int j = 0; j < tmp_vecvec[i].size(); j++) {
            tmp_vecvec[i][j] = ndarray(i, j);
        }
    }
    vecvec = tmp_vecvec;
    return true;
}
