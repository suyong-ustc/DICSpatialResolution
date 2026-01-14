#include <iostream>
#include <format>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
#include "DIC/DICAlgorithm.h"
#include "DIC/DICParameters.h"
#include "Interpolation/Interpolator.h"
#include "Interpolation/BicubicMOMSInterpolator.h"
#include "Interpolation/BiquinticBSplineInterpolatror.h"
#include "Interpolation/BisepticBSplineInterpolator.h"
using namespace arma;


bool ReadImage(const std::string& image_path, mat& image)
{
	std::cout << "Import image with path " << image_path << std::endl;

	// 读取图像
	cv::Mat cvmat = cv::imread(image_path, cv::IMREAD_GRAYSCALE);

	if (cvmat.empty())
	{
		std::cerr << "I can not import the image!" << std::endl;
		return false;
	}

	// 将 opencv 的数据格式转化为 armadillo 矩阵
	image.zeros(cvmat.rows, cvmat.cols);

	for (int r = 0; r < cvmat.rows; ++r)
	{
		for (int c = 0; c < cvmat.cols; ++c)
		{
			image(r, c) = cvmat.at<uchar>(r, c);
		}

	}

	return true;
}



bool RegisterDisplacement(const DICParameters& dic_parameters, int matlab_first_step_position, int matlab_second_step_position)
{
	// 读取图像
	const std::string prefix_image{ "..\\images" };
	const std::string refer_image_path = std::format("{0}\\a_{1}.bmp", prefix_image, matlab_second_step_position);
	const std::string deform_image_path = std::format("{0}\\a_{1}.bmp", prefix_image, 0);

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);

	// 网格点
	mat grid_x;
	mat grid_y;
	ConstructGrid(dic_parameters.roi, dic_parameters.grid_step, grid_x, grid_y);

	// 真实位移
	mat real_u(grid_x.n_rows, grid_x.n_cols, fill::zeros);
	mat real_v(grid_x.n_rows, grid_x.n_cols, fill::zeros);

	int first_step_position = matlab_first_step_position - 2;
	int second_step_position = matlab_second_step_position - 2;

	for (int c = 0; c < grid_x.n_cols; ++c) 
	{
		int x = static_cast<int>(grid_x(0, c));

		if (x < first_step_position)
			real_u.col(c).fill(0);
		else if (x < second_step_position)
			real_u.col(c).fill(1);
		else
			real_u.col(c).fill(2);
	}

	// 初始化计算结果
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order));

	// 设置初值
	dic_output->set_u(real_u);

	// 相关计算
	RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);

	// 计算应变
	mat exx, eyy, exy;
	EvaluateStrain(dic_output->u(), dic_output->v(), dic_parameters.grid_step, dic_parameters.strain_method, exx, eyy, exy);

	// 结果输出
	const std::string prefix_result = std::format("..\\results\\Step{0}Weight{1}Subset{2}Shape{3}_", 
		matlab_second_step_position, dic_parameters.weight, dic_parameters.subset_size, dic_parameters.shape_function_order);

	const std::string x_path = prefix_result + "x.csv";
	const std::string y_path = prefix_result + "y.csv";
	dic_output->write_grid_coordinate(x_path, y_path, arma::csv_ascii);

	const std::string u_path = prefix_result + "u.csv";
	const std::string v_path = prefix_result + "v.csv";
	dic_output->write_displacement_field(u_path, v_path, arma::csv_ascii);

	const std::string iter_times_path = prefix_result + "iter_times.csv";
	dic_output->write_iteration_times(iter_times_path, arma::csv_ascii);

	const std::string zncc_path = prefix_result + "zncc.csv";
	dic_output->write_zncc(zncc_path, arma::csv_ascii);

	const std::string valid_sign_path = prefix_result + "valid_sign.csv";
	dic_output->write_valid_sign(valid_sign_path, arma::csv_ascii);

	const std::string exx_path = prefix_result + "exx.csv";
	const std::string eyy_path = prefix_result + "eyy.csv";
	const std::string exy_path = prefix_result + "exy.csv";
	exx.save(exx_path, csv_ascii);
	eyy.save(eyy_path, csv_ascii);
	exy.save(exy_path, csv_ascii);

	delete dic_output;

	std::cout << " Completed!" << std::endl;

	return true;
}



int main()
{
	// 相关计算参数
	DICParameters dic_parameters;
	dic_parameters.roi = cv::Rect(100, 100, 300, 300);	// 感兴趣区域	
	dic_parameters.grid_step = 1;						// 网格间距	
	dic_parameters.subset_size = 29;					// 子区尺寸	
	dic_parameters.zncc_threshold = 0.8;				// ZNCC阈值	
	dic_parameters.max_iteration_times = 10;			// 最大迭代次数
	dic_parameters.error_threshold = 1e-4;				// 迭代阈值
	dic_parameters.bspline_interpolation_order = 3;		// 插值阶数：3 三次B样条；5 五次B样条；7 七次B样条
	dic_parameters.shape_function_order = 2;			// 形函数阶数
	dic_parameters.weight = 0;							// 权重：0 均匀权重；1 高斯权重
	dic_parameters.strain_method = 5;					// 应变算法：0 Sobel算子；3 三次B样条；5 五次B样条；7 七次B样条
	dic_parameters.thread_number = 12;					// 多线程数目

	// 间断位置
	int matlab_first_step_position = 250;

	// 计算
	for (int matlab_second_step_position = 251; matlab_second_step_position <= 301; matlab_second_step_position++)	
	{
		for (int subset_size = 15; subset_size <= 45; subset_size += 2)
		{
			// 设置子区尺寸
			dic_parameters.subset_size = subset_size;

			// 输出当前参数
			std::cout << std::format("Process Step at {} with Subset Size {}", matlab_second_step_position, dic_parameters.subset_size) << std::endl;

			// 相关计算
			RegisterDisplacement(dic_parameters, matlab_first_step_position, matlab_second_step_position);
		}
	}


}