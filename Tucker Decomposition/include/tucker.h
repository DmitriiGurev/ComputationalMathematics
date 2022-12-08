#pragma once

#include <eigen-3.4.0/unsupported/Eigen/CXX11/Tensor>
#include <eigen-3.4.0/Eigen/Dense>

class Tucker
{
public:
	// Default constructor
	Tucker();

	// Zero tensor with given dimensions and ranks 
	Tucker(int n0, int n1, int n2, int r0, int r1, int r2);

	// Compress a tensor with a given accuracy
	Tucker(const Eigen::Tensor<double, 3>& tensor, double eps = 1e-14, int rmax = 1e+6);

	// Tucker tensor's dimensions, ranks, factors, and core
	std::vector<long long> Dimensions() const;
	std::vector<long long> Ranks() const;
	std::vector<Eigen::MatrixXd> U() const;
	Eigen::Tensor<double, 3> Core() const;

	// Value at (i0, i1, i2)
	double At(int i0, int i1, int i2) const;

	// Full tensor
	Eigen::Tensor<double, 3> Reconstructed() const;

	// Sum of the elements
	double Sum() const;

	// Frobenius norm
	double Norm() const;

	// Recompress the tensor with a given accuracy
	void Recompress(double eps = 1e-14, int rmax = 1e+6);

	// Print the reconstructed tensor
	friend std::ostream& operator <<(std::ostream& out, const Tucker& tucker);

	// Element-wise operations
	friend Tucker operator +(const Tucker& t1, const Tucker& t2);
	friend Tucker operator -(const Tucker& t1, const Tucker& t2);
	friend Tucker operator *(const Tucker& t1, const Tucker& t2);
	friend Tucker operator *(const double alpha, const Tucker& t);
	friend Tucker operator *(const Tucker& t, const double alpha);
	friend Tucker operator -(const Tucker& t);

private:
	// Convert a tensor to a matrix
	template<typename Scalar, int rank, typename sizeType>
	Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
	TensorToMatrix(const Eigen::Tensor<Scalar, rank>& tensor, const sizeType rows, const sizeType cols) const;

	// k-mode unfolding of the tensor  
	Eigen::MatrixXd Unfolding(const Eigen::Tensor<double, 3>& tensor, int k) const;

	// Reconstruct a tensor from the k-mode unfolding
	Eigen::Tensor<double, 3> Folding(int I0, int I1, int I2, const Eigen::MatrixXd& unfolding, int k) const;

	// Compute the factors
	void ComputeU(const Eigen::Tensor<double, 3>& tensor, std::vector<Eigen::MatrixXd>& u, double eps = 1e-14, int rmax = 1e+6);

private:
	std::vector<long long> _n = std::vector<long long>(3);
	std::vector<long long> _r = std::vector<long long>(3);
	std::vector<Eigen::MatrixXd> _u = std::vector<Eigen::MatrixXd>(3);
	Eigen::Tensor<double, 3> _core;
};