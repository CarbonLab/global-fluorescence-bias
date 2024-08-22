% Example data (replace this with your actual data)
% x = [1, 2, 3, 4, 5, 6, 7, 8, 9];
% y = [2, 3, 4, 5, 5, 6, 7, 8, 9];
% fig_flag = 1 to make figure, 0 to skip
function [x_cleaned,y_cleaned,outliers] = chauvenet(x,y,fig_flag)

% Step 1: Fit a linear trend
% coefficients = polyfit(x, y, 1);
% linear_trend = polyval(coefficients, x);
clear mlr; mlr = fitlm(x,y,'Linear','Intercept',false);
linear_trend = mlr.Fitted;

% Step 2: Calculate residuals
residuals = y - linear_trend;

% Step 3: Standardize residuals
std_residuals = residuals ./ std(residuals,'omitnan');

% Step 4: Apply Chauvenet's Criteria
threshold = chauvenet_threshold(length(x));

% Step 5: Remove outliers
outliers = abs(std_residuals) > threshold;
x_cleaned = x(~outliers);
y_cleaned = y(~outliers);

% Step 6: Refit the model
% coefficients_cleaned = polyfit(x_cleaned, y_cleaned, 1);
% linear_trend_cleaned = polyval(coefficients_cleaned, x_cleaned);
clear mlr; mlr = fitlm(x_cleaned,y_cleaned,'Linear','Intercept',false);
linear_trend_cleaned = mlr.Fitted;

% Plotting
if fig_flag == 1
    figure;
    scatter(x, y, 'o', 'DisplayName', 'Original Data');hold on
    hold on;
    plot(x, linear_trend, '-r', 'DisplayName', 'Original Linear Trend');
    scatter(x(outliers), y(outliers), 'x', 'r', 'DisplayName', 'Outliers');
    plot(x_cleaned, linear_trend_cleaned, '-b', 'DisplayName', 'Cleaned Linear Trend');
    xlabel('X');
    ylabel('Y');
    legend('show');
    title('XY Plot with Linear Trend and Outlier Removal');
else
end
    function threshold = chauvenet_threshold(n)
        % Chauvenet's criteria threshold calculation
        p = 0.5; % Probability threshold (adjust as needed)
        z = erfinv(1 - p/(2*n));
        threshold = sqrt(2) * z;
    end

end

