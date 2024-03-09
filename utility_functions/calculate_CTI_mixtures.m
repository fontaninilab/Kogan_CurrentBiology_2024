function CTI = calculate_CTI_mixtures(avg_neural_data)
    L_pairs = nchoosek(1:3,2);
    R_pairs = nchoosek(4:6,2);
    all_pairs = nchoosek(1:6,2);

    idx_diff_direction_pairs = ~ismember(all_pairs,[L_pairs;R_pairs],'rows');
    same_direction_pairs = [L_pairs;R_pairs];
    diff_direction_pairs = all_pairs(idx_diff_direction_pairs,:);
    CTI = zeros(1,size(avg_neural_data,2));
    for i=1:size(avg_neural_data,2)

        data_for_CTI = avg_neural_data(:,i);

        for p=1:size(same_direction_pairs,1) %loop through all pairs of mixtures for L trials
            same_direction_pairs_differences(p) = abs(data_for_CTI(same_direction_pairs(p,1),1) - data_for_CTI(same_direction_pairs(p,2),1));
        end

        for p=1:size(diff_direction_pairs,1) %loop through all pairs of mixtures for L trials
            diff_direction_pairs_differences(p) = abs(data_for_CTI(diff_direction_pairs(p,1),1) - data_for_CTI(diff_direction_pairs(p,2),1));
        end


        category_index_same = mean(same_direction_pairs_differences);
        category_index_difference = mean(diff_direction_pairs_differences);
        %calculate normalized category index
        CTI(1,i) = (category_index_difference - category_index_same )./(category_index_same + category_index_difference);
    end

end