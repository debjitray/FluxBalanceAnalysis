function array = GtZero(array)
for i=1:length(array)
    if (array(i)>0) array(i)=0;end
end;
