function [sum_of_squares_resid, sig_CI, fit_vals_all] = calculate_fits_mixtures(avg_data_for_fit_all)
    
%     %fit functions
%     line2 = {@(b,x) (b(1)*x./100) + b(2),... %linear function
%     @(b,x) b(1)*(x>50)+b(2)*(x<=50),...%step function
%     @(b,x) (b(1)*abs(b(1)*x-b(1)*50)/100)+b(2)}; %V shape

    
    
    % sigfunc = @(B,x) minimum+(maximum-minimum)./(1+10.^((B(1)-x)*B(2)));
%     
    X = [100; 75; 60; 40; 25; 0]; %x values

    for i=1:size(avg_data_for_fit_all,2)
        indv_neuron_data = (squeeze(avg_data_for_fit_all(:,i,:)));
        %fit functions
        minimum = min(indv_neuron_data);
        maximum = max(indv_neuron_data);
        line2 = {@(b,x) (b(1)*x./100) + b(2),... %linear function
            @(b,x) b(1)+(b(2)-b(1))./(1+10.^((50-x)*b(3))),...%sig function
            @(b,x) (b(1)*abs(b(1)*x-b(1)*50)/100)+b(2)}; %V shape
        
        for z=1:3
            if z==2
                init_vals = [minimum maximum 0.01];
            else
                init_vals = [0.01 0.01];
            end
                [fit_vals,resid,~,CovB,~,~] = ...
                    nlinfit(reshape(repmat(X,1,size(indv_neuron_data,2)),1,[]),reshape(indv_neuron_data,1,[]),line2{z},init_vals);
                sum_of_squares_resid(z,i) = sum(resid.^2);
                fit_vals_all{z}(:,i) = fit_vals;
                %confidence interval of fit
                ci_fit = nlparci(fit_vals,resid,'covar',CovB); %each row is the CI for each parameter
                if z==3
                    ci_fit = diff(ci_fit);
                end
                if ismember(length(find(ci_fit(1,:)<=0)),[0 2]) %check if CI includes zero
                    sig_CI(z,i) = 1;
                else
                    sig_CI(z,i) = 0;
                end
                
%             end
            
        end

        null_fit = mean(mean(indv_neuron_data,2));
        sum_of_squares_resid(4,i) = sum((indv_neuron_data - null_fit).^2,'all');
    end
end