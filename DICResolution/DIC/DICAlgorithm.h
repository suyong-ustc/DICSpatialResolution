#pragma once
#include <armadillo>
#include "DICParameters.h"
#include "DICOutput.h"
#include "../Interpolation/Interpolator.h"


/*******************************************************************************************************************/
/**********************************               形函数                ********************************************/
/******************************************************************************************************************/

bool ShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, const int& order, arma::mat& deform_x, arma::mat& deform_y);

bool ZeroOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FirstOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool SecondOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool ThirdOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FourthOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);

bool FifthOrderShapeFunction(const double& x0, const double& y0, const arma::mat& x, const arma::mat& y, const arma::vec& p, arma::mat& deform_x, arma::mat& deform_y);



/*******************************************************************************************************************/
/**********************************               伪逆矩阵                ********************************************/
/******************************************************************************************************************/

bool PseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, const int& order, arma::mat& pseudo);

bool ZeroOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool FirstOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool SecondOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool ThirdOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool FourthOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);

bool FifthOrderPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, arma::mat& pseudo);



bool WeightedPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, const arma::vec& weight, const int& order, arma::mat& pseudo);

bool ZeroOrderWeightedPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, const arma::vec& weight, arma::mat& pseudo);

bool FirstOrderWeightedPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, const arma::vec& weight, arma::mat& pseudo);

bool SecondOrderWeightedPseudoInverseMatrix(const arma::mat& gx, const arma::mat& gy, const arma::mat& x, const arma::mat& y, const arma::vec& weight, arma::mat& pseudo);



/*******************************************************************************************************************/
/**********************************            计算全场变形             ********************************************/
/******************************************************************************************************************/


DICOutput* RigsterFullFieldDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const DICParameters& dic_parameters);

bool EstimateInitialDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const DICParameters& dic_parameters, DICOutput* dic_output);

bool RegisterSubpixelDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const DICParameters& dic_parameters, DICOutput* dic_output);


Interpolator* ConstructInterpolator(const arma::mat& image, const int& n);

void ConstructGrid(const cv::Rect& roi, const int& grid_step, arma::mat& x, arma::mat& y);

void SubsetPointsRelativeCoordinate(const int& m, arma::mat& dx, arma::mat& dy);

arma::mat SubsetIntensities(const arma::mat& image, const int& x0, const int& y0, const int& m);

bool isConverge(const arma::vec& dp, const double& threshold);



/*******************************************************************************************************************/
/**********************************            计算全场应变             ********************************************/
/******************************************************************************************************************/

void EvaluateStrain(const arma::mat& u, const arma::mat& v, int grid_step, int method, arma::mat& exx, arma::mat& eyy, arma::mat& exy);

void EvaluateStrainSobel(const arma::mat& u, const arma::mat& v, int grid_step, arma::mat& exx, arma::mat& eyy, arma::mat& exy);

void EvaluateStrainBSpline(const arma::mat& u, const arma::mat& v, int grid_step, int method, arma::mat& exx, arma::mat& eyy, arma::mat& exy);


/*******************************************************************************************************************/
/**********************************            辅助计算函数              ********************************************/
/******************************************************************************************************************/

int ParameterTotalNumber(const int& order);

double NormalizeVectorize(const arma::mat& m, arma::vec& v);

double ZNCC(const arma::mat& a, const arma::mat& b);