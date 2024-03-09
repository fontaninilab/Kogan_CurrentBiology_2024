function sem = find_sem(A)
if size(A,1) == 1
    sem = std(A,'omitnan')/sqrt(size(A,2));
else
    sem = std(A,'omitnan')/sqrt(size(A,1));

end

end