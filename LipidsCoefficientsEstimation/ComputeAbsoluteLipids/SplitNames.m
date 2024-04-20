function listNames = SplitNames(listLipids)

sterolsList     = {'Campesterol';'Sitosterol';'Stigmasterol'};
sterols         = listLipids((ismember(listLipids(:,1), sterolsList)) == 1, 1);
absoluteNames   = listLipids((ismember(listLipids(:,1), sterolsList)) == 0, 1);
listNames       = split(absoluteNames, ' ');
listNames       = [listNames; [sterols,sterols]];

end