function result = sum_of_squares_diff(lvls)
result = 0;
for i = 1:size(lvls,1)
	for j = 1:(size(lvls,2)-1)
		result = result + (lvls(i,j)-lvls(i,j+1))^2;
	end
end

end