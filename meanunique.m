function [UniquePsthArray] = meanunique(files,varname);
UniquePsthArray = nan(length(files),401);
if( isempty( files ) )
        var = [];
        return;
 end
    for k=1:(length(files))
    var = load( files{k}, varname );
    var = var.(varname);
    var = unique(var,'rows','stable');
    UniquePsthArray(k,:) = nanmean(var);
   
    end

