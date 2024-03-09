function [xvals, b1] = bar_custom_mixtures(varargin)



    b1 = bar(varargin{1});
    % for i=1:length(b1)
    if length(b1)>1
        b1(1).FaceColor = [0 0 0 ];
        b1(1).EdgeColor = [0 0 0];
        b1(1).LineWidth = .5;
        b1(1).CData = repmat([0 0 0],size(b1(1).CData,1),1);

        b1(2).FaceColor = [0 1 1 ];
        b1(2).EdgeColor = [0 0 0];
        b1(2).LineWidth = .5;
        b1(2).CData = repmat([0 0 0],size(b1(2).CData,1),1);
    else
        b1(1).FaceColor = 'flat';
        b1(1).EdgeColor = [0 0 0];
        b1(1).LineWidth = .5;
        b1(1).CData(1,:) = [1 1 1];
        b1(1).CData(2,:) = [.5 .5 .5];
        

    end
    % b1(2).LineStyle = {'-',':'};

    box off
    for i=1:length(b1)
        xvals(i,:) = b1(i).XData+b1(i).XOffset;
    end

end