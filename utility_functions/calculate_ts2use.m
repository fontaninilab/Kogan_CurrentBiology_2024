function ts2use = calculate_ts2use(time_stamps, window2use)
    ts2use = intersect(find(time_stamps>=window2use(1)), find(time_stamps<=window2use(2)));

end