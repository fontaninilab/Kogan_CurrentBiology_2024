function [performance,DECODING_RESULTS] = ndt_decoding_mixtures_v2(binned_format_file_name,shuffle, lbls2use,specific_label_name_to_use,num_neurons2use)
    toolbox_basedir_name = 'D:\Code_data_analysis\ndt.1.0.4';
    addpath(toolbox_basedir_name);
    % add the NDT paths using add_ndt_paths_and_init_rand_generator
    add_ndt_paths_and_init_rand_generator


    % will decode the identity of which object was shown (regardless of its position)
%     specific_label_name_to_use = 'trialID';
    %         if numneurons(qq)>=a(i)
    num_cv_splits = 10; %changed from 10 to 5
%     load_binned_data_and_convert_firing_rates_to_spike_counts 
    ds = basic_DS(binned_format_file_name, specific_label_name_to_use, num_cv_splits);
    
    %this will ignore neurons with too few trials
    ds.sites_to_use = find_sites_with_k_label_repetitions(ds.the_labels, num_cv_splits);
    ds.randomly_shuffle_labels_before_running = shuffle; %can shuffle labels by setting this parameter to 1
    ds.label_names_to_use = lbls2use; %only use specific tastes
    the_feature_preprocessors={};
    
    %can select thresholds for including/excluding neurosn in analysis
%     the_feature_preprocessors{1} = select_pvalue_significant_features_FP;
%     the_feature_preprocessors{1}.pvalue_threshold = 0.05;
%     the_feature_preprocessors{1}.save_extra_info = 0;
%     the_feature_preprocessors{2} = zscore_normalize_FP;
%     the_feature_preprocessors{1} = select_or_exclude_top_k_features_FP;
    
%     the_feature_preprocessors{1}.num_features_to_exclude = 10;
%     the_feature_preprocessors{1}.num_features_to_use = num_neurons2use;
%     the_feature_preprocessors{1}.save_extra_info = 1;
    the_feature_preprocessors{1} = zscore_normalize_FP;    

    the_classifier = max_correlation_coefficient_CL;
    the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
    the_cross_validator.test_only_at_training_times =1;
    
    % set how many times the outer 'resample' loop is run
    % generally we use more than 2 resample runs which will give more accurate results
    % but to save time in this tutorial we are using a small number.

    the_cross_validator.num_resample_runs =50; 
    the_cross_validator.confusion_matrix_params.create_confusion_matrix =1;

    DECODING_RESULTS = the_cross_validator.run_cv_decoding;
    performance = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results;
end
