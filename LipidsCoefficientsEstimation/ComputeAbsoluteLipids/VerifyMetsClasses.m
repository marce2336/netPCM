function uniqueNamesNew = VerifyMetsClasses(uniqueNames)

listClasses = {'PC';'PE';'PG';'PI';'PS';'PA';'LPA';'LPC';'PGP';'MG';'DG';...
    'TG';'MGMG';'MGDG';'DGDG';'CL';'SQDG';'CDP-DG';'WE';...
    'FFA';'FA';'FAL';'FOH';'HC';'OHC';'Oxo';'DCA';'HFA';...
    'SPB';'Cer';'GlcCer';'GlcA-IPC';'Ara-GlcA-IPC';'Gal-GlcA-IPC';'Glc-GlcA-IPC';'Man-GlcA-IPC';'PI-Cer';...
    'Campesterol';'Sitosterol';'Stigmasterol';'ST';'Brassicasterol'};

[idxMiss,~] = ismember(listClasses, uniqueNames);

missingClasses = listClasses(idxMiss == 0);

uniqueNamesNew = [uniqueNames;missingClasses];

end