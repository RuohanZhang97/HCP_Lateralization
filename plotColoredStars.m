function plotColoredStars(x_star, y_star, LH, RH)
    for i = 1:length(x_star)
        idx = x_star(i);
        if RH(idx) > LH(idx)
            starColor = [0.6350 0.0780 0.1840];   % Right > Left
        else
            starColor = [0 0.4470 0.7410];  % Left > Right
        end
        text(idx, y_star(i), '*', 'FontSize', 18, ...
             'Color', starColor, 'HorizontalAlignment', 'center');
    end
end