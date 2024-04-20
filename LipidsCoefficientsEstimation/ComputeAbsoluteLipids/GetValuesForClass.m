function absoluteClass_i = GetValuesForClass(absoluteNames, abundanceValues, species)

getIdx = ismember(absoluteNames(:,1), species);
absoluteClass_i = abundanceValues(getIdx, :);
absoluteClass_i = (min(min(absoluteClass_i)))/5;

end