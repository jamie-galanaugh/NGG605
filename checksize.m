function varsize = checksize(file, var)
varsize = load(file, var);
varsize=cell2mat(struct2cell(varsize));
varsize = size(varsize,1);
