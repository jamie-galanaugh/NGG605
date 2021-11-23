function[Where] = Where_file(filename);

Where = strvcat(filename);
pos = strfind(Where, '\');
pos = max(pos);
Where = Where(1,1:pos);

end
