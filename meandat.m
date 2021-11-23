function mns = meandat(files, varname)

mns = nan(length(files),1);
 if( isempty( files ) )
        var = [];
        return;
 end
    for k=1:(length(files))
    var = load( files{k}, varname );
    var = var.(varname);
    var = mean(var);
    mns(k,:) = var;
    end
    
