function result = sum_of_squares(lvls)
    result = 0;
    for i = 1:size(lvls,1)
        for j = 1:size(lvls,2)
            result = result + lvls(i,j)^2;
        end
    end
end