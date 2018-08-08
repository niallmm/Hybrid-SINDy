function [II_Xi, II_ind, Xi_ind] = build_Xi_Library(Xicomb, II_Xi)

if ~iscell(Xicomb)
    XiMat = reshape(Xicomb, [],1);
else
    for kk = 1:length(Xicomb) % for each unique returned coefficient matrix
        XiMat(:,kk) = reshape(Xicomb{kk}, [],1);
    end
end

IItest = unique((abs(XiMat)>0)','rows')'; % make a structure vector of 0s and 1s

if isempty(II_Xi) % add the first vector to the structure library
    II_Xi = IItest;
else

    Lia = ismember(IItest', II_Xi', 'rows');

    II_Xi = [II_Xi IItest(:,~Lia')]; % if structurally unique add to library%
    
end

II_ind = find(ismember(II_Xi',IItest','rows'));
Xi_ind = find(ismember(IItest', II_Xi','rows'));
