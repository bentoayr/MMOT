% generate plots for hypergraph clustering
max_valid = number_of_repetitions;
histdata = nan(2,3,max_valid);
for method = 1:2
    for type = 1:3
        best_cut = inf;
        for th_id = 1:length(thresholds_lev_range)
            tmp = nan(max_valid,1);
            for i = 1:max_valid
                tmp(i) = all_errors_hyper{i}{type,th_id}(method);
            end
            if (best_cut > mean(tmp))
                best_cut = mean(tmp);
                histdata(method,type,:) = tmp;
            end
        end
    end
end

all_data = nan(max_valid,6);
all_data(:,1) = histdata(1,1,:);
all_data(:,2) = histdata(1,2,:);
all_data(:,3) = histdata(1,3,:);

all_data(:,4) = histdata(2,1,:);
all_data(:,5) = histdata(2,2,:);
all_data(:,6) = histdata(2,3,:);


figure;
hold on;
histogram(all_data(:,1),'NumBins',10,'FaceAlpha',0.5);
histogram(all_data(:,2),'NumBins',10,'FaceAlpha',0.5);
histogram(all_data(:,3),'NumBins',10,'FaceAlpha',0.5);

xlabel('Fraction of miss-classified graphs');
ylabel('Number of repetitions');
title('Hypergraph clustering via NH-Cut');
n_metric_NH_cut_mean = mean(all_data(:,1));
non_n_metric_NH_cut_mean = mean(all_data(:,2));
barycenter_NH_cut_mean = mean(all_data(:,3));
n_metric_NH_cut_std = std(all_data(:,1));
non_n_metric_NH_cut_std = std(all_data(:,2));
barycenter_NH_cut_std = std(all_data(:,3));

legend({['n-metric (mean ', num2str(round(n_metric_NH_cut_mean,3)), ')'],['non-n-metric (mean ', num2str(round(non_n_metric_NH_cut_mean,3)), ')'],['barycenter (mean ', num2str(round(barycenter_NH_cut_mean,3)), ')']},'Location','northwest');
box on;
set(gca,'fontname','times');
set(gca,'fontsize',14);

figure;
hold on;
histogram(all_data(:,4),'NumBins',10,'FaceAlpha',0.5);
histogram(all_data(:,5),'NumBins',10,'FaceAlpha',0.5);
histogram(all_data(:,6),'NumBins',10,'FaceAlpha',0.5);

xlabel('Fraction of miss-classified graphs');
ylabel('Number of repetitions');
title('Hypergraph clustering via TTM');
legend({'n-metric','non-n-metric'},'Location','northwest');
n_metric_TTM_mean = mean(all_data(:,4));
non_n_metric_TTM_mean = mean(all_data(:,5));
barycenter_TTM_mean = mean(all_data(:,6));
n_metric_TTM_std = std(all_data(:,4));
non_n_metric_TTM_std = std(all_data(:,5));
barycenter_TTM_std = std(all_data(:,6));


legend({['n-metric (mean ', num2str(round(n_metric_TTM_mean,3)), ')'],['non-n-metric (mean ', num2str(round(non_n_metric_TTM_mean,3)), ')'],['barycenter (mean ', num2str(round(barycenter_TTM_mean,3)), ')']},'Location','northwest');
box on;
set(gca,'fontname','times');
set(gca,'fontsize',14);

% generate plots for spectral clustering
all_data = nan(max_valid,1);
for i = 1:max_valid
    tmp = (all_errors_classic{i});
    all_data(i) = tmp(1);
end

figure;
hold on;
histogram(all_data,'NumBins',10,'FaceAlpha',0.5);
xlabel('Fraction of miss-classified graphs');
ylabel('Number of repetitions');
title('Spectral clustering');
spectral_mean = mean(all_data);
spectral_std = std(all_data);

legend({['WD (mean ', num2str(round(spectral_mean,3)),')']},'Location','northwest');
box on;
set(gca,'fontname','times');
set(gca,'fontsize',14);
