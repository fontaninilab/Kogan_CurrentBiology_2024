function [sig_neurons_f, f_stat] = calculate_f_stat_mixtures_fit(sum_of_squares_resid)

    for z=1:3
        if ismember(z,[1 3])
            p=2;
        else
            p=3;
        end
        x = 0:0.001:50;
        x=x(2:end);
        y = fpdf(x,p-1,6-p);
        f_cutoff(z) = x(min(find(y(2:end)<0.01))); %pval cutoff for fit

        for i=1:size(sum_of_squares_resid,2)
            f_stat(z,i) = ((sum_of_squares_resid(4,i) - sum_of_squares_resid(z,i))/(p-1))...
                /(sum_of_squares_resid(z,i)/(6-p));
        end
        sig_fits_f(z) = length(find(f_stat(z,:)>f_cutoff(z)));

    end
    [~,idx_max] = max(f_stat(1:3,:));
% 
    for z=1:3
        sig_neurons_f{1,z} = find(f_stat(z,:)>f_cutoff(z));
        sig_neurons_f{2,z} = intersect(find(idx_max == z), find(f_stat(z,:)>f_cutoff(z)));

%         sig_neurons_f{3,z} = find(idx_max == z);
    end
    sig_neurons_f{1,4} = setdiff(1:size(f_stat,2), [sig_neurons_f{2,1} sig_neurons_f{2,2} sig_neurons_f{2,3}]);
    sig_neurons_f{2,4} = setdiff(1:size(f_stat,2), [sig_neurons_f{2,1} sig_neurons_f{2,2} sig_neurons_f{2,3}]);
end