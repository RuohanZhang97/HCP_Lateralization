function [tVal, pVal, dz, mean_scene, mean_other, diff_mean] = run_within_hemi_contrast(X_scene, X_other, task_cov)
% 对单一半球进行 Scenes 与 Other 的直接比较（配对）
% X_scene / X_other: (Region x Subject)
% task_cov: (Subject x Covariates)
    NumRegion_local = size(X_scene,1);
    NumSub_local    = size(X_scene,2);
    tVal = zeros(NumRegion_local,1);
    pVal = zeros(NumRegion_local,1);
    dz   = zeros(NumRegion_local,1);
    mean_scene = mean(X_scene,2);
    mean_other = mean(X_other,2);
    diff_mean  = mean_scene - mean_other;   % 方向：Scenes - Other

    % 拼接两条件，一次回归取残差，再按条件拆分
    X_cov = [task_cov; task_cov];                  % (2N x K)
    X_cov_aug = [ones(size(X_cov,1),1) X_cov];     % 加截距

    for rr = 1:NumRegion_local
        y_concat = [X_scene(rr,:) X_other(rr,:)]'; % (2N x 1)
        % 回归协变量，取残差
        [~,~,r] = regress(y_concat, X_cov_aug);
        res_scene = r(1:NumSub_local);
        res_other = r(NumSub_local+1:end);

        % 配对 t 检验
        [~,p,~,stats] = ttest(res_scene, res_other);
        pVal(rr) = p;
        tVal(rr) = stats.tstat;

        % Cohen's dz（配对设计）：差值的均值 / 差值的标准差（N-1）
        d_vec = res_scene - res_other;
        dz(rr) = mean(d_vec) / std(d_vec, 0);
    end
end