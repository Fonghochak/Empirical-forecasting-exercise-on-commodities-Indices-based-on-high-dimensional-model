% 读取数据
targetData = readtable('MonthlyTarget.xlsx');
financialData = readtable('Financial variables.xlsx');
transFredData = readtable('Trans FRED.xlsx');

% 移除日期列，保留数值列
targetData = targetData(:, 2:end);

% 清理列名
targetData.Properties.VariableNames = matlab.lang.makeValidName(targetData.Properties.VariableNames);

% 计算对数收益率并转换为百分比收益率
logReturns = varfun(@(x) [NaN; diff(log(x))] * 100, targetData);
logReturns.Properties.VariableNames = strcat('log_return_', targetData.Properties.VariableNames);

% 计算滚动波动率
windowSize = 12;
rollingVolatility = varfun(@(x) movstd(x, windowSize, 'omitnan'), logReturns);
rollingVolatility.Properties.VariableNames = strcat('volatility_', logReturns.Properties.VariableNames);

% 创建一个图形窗口
figure;

% 绘制价格曲线
subplot(3, 1, 1);
plot(targetData{:, :});
title('Monthly Prices');
xlabel('Time');
ylabel('Price');
% 只显示图例名称，去除多余的注释
legend(targetData.Properties.VariableNames, 'Location', 'best');
grid on;

% 绘制滚动波动率曲线
subplot(3, 1, 2);
plot(rollingVolatility{:, :});
title('Monthly Rolling Volatility');
xlabel('Time');
ylabel('Volatility');
% 只显示图例名称，去除多余的注释
legend(rollingVolatility.Properties.VariableNames, 'Location', 'best');
grid on;

% 绘制对数收益率曲线
subplot(3, 1, 3);
plot(logReturns{:, :});
title('Log Returns');
xlabel('Time');
ylabel('Percentage Log Returns');
% 只显示图例名称，去除多余的注释
legend(logReturns.Properties.VariableNames, 'Location', 'best');
grid on;

% 保存合并图
saveas(gcf, 'Combined_Charts.png');

disp('Combined chart has been generated and saved as Combined_Charts.png.');
%%
% Clear workspace and command window
clear;
clc;

% Define file path for the combined data file
combinedDataPath = 'E:/OneDrive - Lancaster University/Msc MBF/2. ECO420 Dissertation/Week 3 Data/Final data/Data for matlab/stationaryfilled.xlsx';

% Check if the file exists
if ~isfile(combinedDataPath)
    error('File does not exist. Please check the file path.');
end

% Read combined data
dataTable = readtable(combinedDataPath, 'ReadVariableNames', true);

% Assume the first three columns are the target log returns
dataArray = table2array(dataTable(:, 1:3));
variableNames = dataTable.Properties.VariableNames(1:3);

% Generate the plot
figure;
hold on;
plot(dataArray(:,1), 'LineWidth', 1.5);
plot(dataArray(:,2), 'LineWidth', 1.5);
plot(dataArray(:,3), 'LineWidth', 1.5);
hold off;

% Add title and labels
title('Combined Plot of Target Log Returns');
xlabel('Time');
ylabel('Log Returns');

% Add legend
legend(variableNames, 'Location', 'best');

% Display the grid
grid on;

% Save the figure as an image
saveas(gcf, 'CombinedLogReturnsPlot.png');
%%
% 读取图片
img1 = imread('E:\OneDrive - Lancaster University\Msc MBF\2. ECO420 Dissertation\Week 3 Data\1.png');
img2 = imread('E:\OneDrive - Lancaster University\Msc MBF\2. ECO420 Dissertation\Week 3 Data\2.png');
img3 = imread('E:\OneDrive - Lancaster University\Msc MBF\2. ECO420 Dissertation\Week 3 Data\3.png');
img4 = imread('E:\OneDrive - Lancaster University\Msc MBF\2. ECO420 Dissertation\Week 3 Data\4.png');
img5 = imread('E:\OneDrive - Lancaster University\Msc MBF\2. ECO420 Dissertation\Week 3 Data\5.png');

% 确保所有图片具有相同的宽度
width = size(img1, 2);
img2 = imresize(img2, [NaN width]);
img3 = imresize(img3, [NaN width]);
img4 = imresize(img4, [NaN width]);
img5 = imresize(img5, [NaN width]); % 这里调整第五张图片的宽度

% 将前四张图片拼接成2x2的网格
row1 = cat(2, img1, img2); % 拼接第一行
row2 = cat(2, img3, img4); % 拼接第二行
grid2x2 = cat(1, row1, row2); % 拼接成2x2网格

% 将第五张图片与2x2网格垂直拼接
combined_image = cat(1, grid2x2, img5);

% 显示组合后的图片
imshow(combined_image);

% 保存组合后的图片
imwrite(combined_image, 'E:\OneDrive - Lancaster University\Msc MBF\2. ECO420 Dissertation\Week 3 Data\6.png');
