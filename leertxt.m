function [conductas2]=leertxt(filename);

%% loading the file
fid = fopen(filename);%load the file and open it
    C = textscan(fid, '%s', 'Delimiter', '\n'); C = C{1};%reads data from an open text file into an array
    original_file = C;  
    fclose(fid);
    n = length(C);
    data = cell(0, 2);
%% finding the boxes

     for i=1:n% is traversing the variable n
        line_i = C{i}; %is separating them by matrix
        pos = strfind(line_i, 'Box:');%  it searches within a box string and returns a 1 when it finds it
        if (pos == 1) % enter a conditional if pos is 1 create a matrix with the number of equals line_i 
            Box_number1{i,1} = line_i; %create a matrix that goes from i down
            nn(i,1) =1;
        end
     end
   
     nnn= find (nn==1); %find the variable nn and create a new variable box number1
     Box_number1 = Box_number1(nnn);
    
 %% Getting the C array    
    for i=1:n 
        line_i = C{i};
        % Obtain the position of the character :
        pos = strfind(line_i, ':');% now search inside the characters:
        if (length(pos) == 1) %enters a conditional and asks it when poss is 1 separates the values ​​and the variables.
            % Separate into variable and value
            vari = line_i(1:(pos - 1)); %inside the variable line_i from a position to pos-1
            vari = strrep(vari, ' ', '');   %here search and replace within the variable vari       
            % Check if vari is a number
            val_vari = str2double(vari);%to separate the values ​​if they are numbers
            val = line_i((pos + 1):end);
            % Check if variiable is a single letter
            if ((length(vari) == 1) || isfinite(val_vari)) % returns the value of the array of the size of vari
                [x] = str2num(val);% generates an array with the value of val
                k = size(data, 1) + 1;% generate a variable k that gives you the size of an array
                if (isfinite(val_vari))%asks if the matrix is ​​infinite and creates a matrix that is dependent de k ya sea de val_ vari o vari
                    data{k, 1} = val_vari;
                else
                    data{k, 1} = vari;
                end;
                data{k, 2} = x;
            end;
        end;
    end;
    
    % Process the matrix variiables
    new_data = {}; %create a new empty array
    cur_vari = '';
    temp = [];
    for i=1:size(data, 1) %It goes through data asked if it is a number, save it in temp or if it is empty 
        if (isnumeric(data{i, 1}))
            temp = [temp, data{i, 2}];
        elseif (isempty(data{i, 2})) 
            if (~isempty(temp)) %If it is empty, return k and add 1
                k = size(new_data, 1) + 1;
                new_data{k, 1} = cur_val;
                new_data{k, 2} = temp;
                cur_val = ''; temp = []; %it is obtaining the values ​​of all the measured variables
            end;
            cur_val = data{i, 1}; 
        else
            if (~isempty(temp))
                k = size(new_data, 1) + 1;
                new_data{k, 1} = cur_val;
                new_data{k, 2} = temp;
                cur_val = ''; temp = [];
            end;
            k = size(new_data, 1) + 1;
            new_data{k, 1} = data{i, 1};
            new_data{k, 2} = data{i, 2}; %It generates a matrix that separates the numbers from the letters and makes the data data available
        end;
    end;
    if (~isempty(temp)) % Flush the final data
        k = size(new_data, 1) + 1;
        new_data{k, 1} = cur_val;
        new_data{k, 2} = temp;
    end;

%% finding the C array from the new_data variable 
    
    CC = new_data(:,1);
    
     for i=1:length (CC)
        line_i = CC{i};
        pos = strfind(line_i, 'C'); %You search inside the matrix CC until you find C, which is the one with the behavior data
        if (pos == 1)
            Cs_number1{i,1} = line_i; %when it finds it creates another variable 
            qq(i,1) =1;
        end
     end
   
     qqq= find (qq==1);
     Cs =  qqq;


%% Doing the analyses
a= length(Cs);
test(1,1)= new_data(Cs(1,1),2);
test= reshape (test, 1,1);
box1=test{1,1}';
box1(10000,1)=NaN;
                              %%generates a variable that is test where
                              %%accommodates the behavior data and then
                              %%takes the transpose of that matrix
Boxes = [box1]; 

conductas=Boxes;
aaaa=find(conductas==0);
conductas2=conductas(1:aaaa(1,1)-1,1);
end