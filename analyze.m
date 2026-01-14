clc;
clear;
close all;

%% 读取图片
% refer_image = imread("..\images\a_280.bmp");
% target_image = imread("..\images\a_0.bmp");
% refer_image = rgb2gray(refer_image);
% target_image = rgb2gray(target_image);

%% 参数
matlab_first_step_position = 250;
matlab_second_step_position = 280;
first_step_position = matlab_first_step_position - 2;      % 第一个间断位置
second_step_position = matlab_second_step_position - 2;     % 第二个间断位置
subset_size = 45;   % 子区尺寸
N = 2;              % 形函数阶数
weight = 0;         % 权重

M = (subset_size - 1) / 2;
m = M + 0.5;
R = M / 3;

%% 读取数据
prefix = ['results\\Step', num2str(matlab_second_step_position), ...
    'Weight', num2str(weight), ...
    'Subset', num2str(subset_size), ...
    'Shape', num2str(N), '_'];
x = readmatrix([prefix, 'x.csv']);
y = readmatrix([prefix, 'y.csv']);
u = readmatrix([prefix, 'u.csv']);
v = readmatrix([prefix, 'v.csv']);

u_dic_avg = mean(u);
u_dic_std = std(u);

e_bspline = readmatrix([prefix, 'exx.csv']);
e_bspline_dic_avg = mean(e_bspline,'omitnan');
e_bspline_dic_std = std(e_bspline,'omitnan');

%% 计算应变
GX = [1, 0, -1; 
      2, 0, -2; 
      1, 0, -1] / 8;
e_sobel = conv2(u, GX, 'same');
e_sobel(:,end) = 0;

e_soble_dic_avg = mean(e_sobel);
e_soble_dic_std = std(e_sobel);

figure;
imagesc(e_sobel);
colorbar;

%% 比较计算值与真实值

xx = x(1,:);
xxx = linspace(xx(1), xx(end), 2000);

if weight == 0
    [u_theory1, e_theory1] = DICStepDisplacementStrain(xxx, first_step_position - 0.5, m, N, weight);
    [u_theory2, e_theory2] = DICStepDisplacementStrain(xxx, second_step_position - 0.5, m, N, weight);
elseif weight == 1
    [u_theory1, e_theory1] = DICStepDisplacementStrain(xxx, first_step_position - 0.5, R, N, weight);
    [u_theory2, e_theory2] = DICStepDisplacementStrain(xxx, second_step_position - 0.5, R, N, weight);
end

u_theory = u_theory1 + u_theory2;
e_theory = e_theory1 + e_theory2;

figure;
hold on;
errorbar(xx, u_dic_avg, u_dic_std);
plot(xxx, u_theory);

figure;
hold on;
% errorbar(xx, e_bspline_dic_avg, e_bspline_dic_std,'Color','blue');
errorbar(xx, e_soble_dic_avg, e_soble_dic_std,'Color','black');
plot(xxx, e_theory, 'r');
legend('Sobel','DIC');
