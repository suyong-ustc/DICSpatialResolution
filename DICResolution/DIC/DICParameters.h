#pragma once
#include <opencv2/core.hpp>
#include <armadillo>


class DICParameters
{
public:
	DICParameters() {};
	~DICParameters() {};

	// 感兴趣区域
	cv::Rect roi{ cv::Rect{100, 100, 300, 300} };

	// 网格间距
	int grid_step { 1 };

	// 子区尺寸
	int subset_size { 19 };

	// ZNCC 阈值
	double zncc_threshold { 0.8 };

	// 最大迭代次数
	int max_iteration_times { 10 };

	// 迭代阈值
	double error_threshold { 2e-4 };

	// 插值阶数
	int bspline_interpolation_order { 3 };

	// 形函数阶数
	int shape_function_order { 0 };

	// 权重
	int weight { 0 };

	// 应变算法
	int strain_method { 0 };

	// 线程数目
	int thread_number { 8 };

};

