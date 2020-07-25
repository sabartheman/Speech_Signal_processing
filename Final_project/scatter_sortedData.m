close all
[height,width] = size(sortedData);

for i = 1:height
   
   figure()
   hold on
    for j = 1:width
       temp = cell2mat(sortedData(i,j));
       if(~isempty(temp))
           temp;
           scatter(temp(1),temp(2))
           title(charlist{i,:})
   
       end
   end
end