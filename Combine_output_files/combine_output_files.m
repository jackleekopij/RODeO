%% Combines output files from price-taker model
% 1. looks at every csv file in the selected directory
% 2. Combines the files according to the GAMS price-taker model format
%    (Storage_dispatch_FC-EY.gms)
% 3. Outputs a CSV file with combined outputs
clear all
close all
clc

dir1 = 'C:\Users\jeichman\Documents\gamsdir\projdir\RODeO\Projects\Test\';
year_val = 2017;    % Initialize year
path_name2 = [dir1,'Output\'];
Select_file = 0;    % 0=utility rates, 1=nodes

%%% Loads all csv files in a given folder
    files2load = dir(path_name2); files2load2={files2load.name}';  % Identify files in a folder    
    for i0=1:length(files2load2) % Remove items from list that do not fit criteria
        load_file1(i0)=~isempty(strfind(files2load2{i0},'.csv'));       % Find only CSV files
        load_file1(i0)=~isempty(strfind(files2load2{i0},'summary_'));   % Find only Summary files
        load_file2(i0)=~isempty(strfind(files2load2{i0},'.csv'));       % Find only CSV files
        load_file2(i0)=~isempty(strfind(files2load2{i0},'dispatch_results'));    % Find only Results files
        load_file3(i0)=~isempty(strfind(files2load2{i0},'.csv'));       % Find only CSV files
        load_file3(i0)=~isempty(strfind(files2load2{i0},'dispatch_inputs'));    % Find only Results files
    end 
    files2load2results=files2load2(load_file2); % Collect file names to load
    files2load2inputs=files2load2(load_file3);
    files2load2=files2load2(load_file1);    clear load_file1 load_file2 load_file3

    try
        [num1,txt1,raw1] = xlsread([path_name2,'Time_data_2015.xlsx']);
    catch
    end
    
%% Use following lists to set technologies each run
in_tech = {'EY','SMR_'};
out_tech = {'HYPS','Batt','FC','H2O2'};
AS_type = {'Eonly','All'};
AS_type2 = {'Eonly','All'};

for i2=1:length(in_tech)
    interim1 = strfind(files2load2,in_tech{i2});
    for i3=1:length(files2load2), in_tech_match(i3,i2) = ~isempty(interim1{i3}); end
end
for i2=1:length(out_tech)
    interim1 = strfind(files2load2,out_tech{i2});
    for i3=1:length(files2load2), out_tech_match(i3,i2) = ~isempty(interim1{i3}); end
end    
for i2=1:length(AS_type)
    interim1 = strfind(files2load2,AS_type{i2});
    for i3=1:length(files2load2), AS_type_match(i3,i2) = ~isempty(interim1{i3}); end
end 
H2_stor = ((in_tech_match(:,find(strcmp(in_tech,'EY')))+in_tech_match(:,find(strcmp(in_tech,'SMR_'))))>0);    %% +out_tech_match(:,find(strcmp(out_tech,'FC'))) (Removed b/c having FC doesn't mean you need storage)

for i0=1:length(files2load2)
    fid = fopen([path_name2,char(files2load2(i0))],'rt'); % Read and parse file
    C1 = textscan(fid, '%s %s','Delimiter',',','HeaderLines',0,'MultipleDelimsAsOne',true, 'CollectOutput',false);
    fclose(fid);
    category1 = [C1{1}];
    category2 = [{'Scenario'};category1(1:2,:);' ';category1(3:17,:);' ';category1(18:24,:);' ';category1(25:31,:);' ';category1(32:end,:)];     %;' ';{'Input Tech'};{'Output Tech'};{'H2 Storage'};{'Short Name'}];
%     category2 = [{'Scenario'};category1(1:2,:);' ';category1(3:17,:);' ';category1(18:24,:);' ';category1(25:31,:);' ';category1(32:end,:);' ';{'Input Tech'};{'Output Tech'};{'H2 Storage'};{'Short Name'}];

    data1 = C1{2};
    filename2 = files2load2{i0};    
    interim2 = {''};
    if isempty(out_tech(find(out_tech_match(i0,:)))),
        last_rows(2) = {' '};
    else last_rows(2) = out_tech(find(out_tech_match(i0,:)));
    end
    if isempty(in_tech(find(in_tech_match(i0,:)))),
        last_rows(1) = {' '};
    else last_rows(1) = in_tech(find(in_tech_match(i0,:))); 
        if strcmp(last_rows(2),{' '}), interim2 = {''};
        else                           interim2 = {'-'};
        end
    end
    last_rows(3)={[last_rows{2},interim2{:},last_rows{1},' ',AS_type2{AS_type_match(i0,:)}]};       
    data2(:,i0) = [filename2(26:end-4);data1(1:2,:);' ';data1(3:17,:);' ';data1(18:24,:);' ';data1(25:31,:);' ';data1(32:end,:)];    %;' ';last_rows{1};last_rows{2};num2str(H2_stor(i0));last_rows{3}];
%     data2(:,i0) = [filename2(26:end-4);data1(1:2,:);' ';data1(3:17,:);' ';data1(18:24,:);' ';data1(25:31,:);' ';data1(32:end,:);' ';last_rows{1};last_rows{2};num2str(H2_stor(i0));last_rows{3}];
    disp(['Summary ',num2str(i0),' of ',num2str(length(files2load2))])
end

%% Turn on/off parsing of results and inputs
%%%%%%%
% Might need to fix after adjusting files in GAMS
%%%%%%%

%%% Load Output
if 1==0
    header_data3 = {'Operation','Utility','Rate','Connection','Renewable Type','Capacity Factor','Installed renewables (% of cap)','Storage duration'};
    B2data_header = {'Interval','Input Pwr (MW)','Output Pwr (MW)','Storage Level (MW-h)','Input Reg Up (MW)','Output Reg Up (MW)','Input Reg Dn (MW)','Output Reg Dn (MW)','Input Spin Res (MW)','Output Spin Res (MW)','Input Nonspin (MW)','Output Nonspin (MW)','H2 Out (kg)','Renewable Input (MW)','Nonrenewable Input (MW)'};
    tic
    c5=0;
    files2load2results2 = {};
    for i1=1:length(files2load2results)
        if ~isempty(strfind(files2load2results{i1},'_4hrs'))
            c5=c5+1;
            files2load2results2(c5,1) = files2load2results(i1);
        elseif ~isempty(strfind(files2load2results{i1},'_Base_'))
            c5=c5+1;
            files2load2results2(c5,1) = files2load2results(i1);
        end        
    end
    B2data = zeros(8760,length(B2data_header),length(files2load2results2)); % Initialize data matrix (row=time, col=property, depth=scenario)
%     B2datacell = cell(8760,length(B2data_header),length(files2load2results)); % Initialize data matrix (row=time, col=property, depth=scenario)
%     B2data_add = cell(8760,length(header_data3),length(files2load2results)); % Initialize additional data matrix (row=time, col=property, depth=scenario)
    B2datacell = cell(8760*length(files2load2results2),length(B2data_header)); % Initialize data matrix (row=time, col=property, depth=scenario)
    B2data_add = cell(8760*length(files2load2results2),length(header_data3)); % Initialize additional data matrix (row=time, col=property, depth=scenario)

    
% % % %     
% % % %     
% % % %     fid = fopen([path_name2,'Combined_operation_data_matlab',num2str(i0),'.csv'],'wt'); % Write output to file
% % % %     s1 = fgetl(fid);
% % % %     fprintf(fid,[repmat(['%s,'],1,n1-1),'%s'],C2{:,:});
% % % %     fprintf(fid,'\n');
% % % %     fclose(fid);
% % % %     
    for i0=1:length(files2load2results2)
        B2 = csvread([path_name2,char(files2load2results2(i0))],28,0,[28 0 8787 14]);        
        
% % % % % %         fid = fopen([path_name2,char(files2load2results2(i0))],'rt'); % Read and parse results file
% % % % % %         A1 = fscanf(fid,'%c');                  %Read data and store in A
% % % % % %         fclose(fid);                            %Close file
% % % % % %         A1_line_break = regexp(A1,'\n');        %Determine where line breaks are (Removed \r\n for HPC runs)
% % % % % %         B1 = cell(length(A1_line_break),1);     %Initialize matrix to hold text file separated by lines
% % % % % %         B2 = cell(length(A1_line_break),13);    %Initialize matrix to hold text file separated by lines
% % % % % %         for i=1:length(A1_line_break)           %Separate text file by lines and save in B
% % % % % %             if i==1, B1(i) = cellstr(A1(1:A1_line_break(i)));
% % % % % %             else     B1(i) = cellstr(A1(A1_line_break(i-1):A1_line_break(i)));
% % % % % %             end
% % % % % %             Find_str_B1 = findstr(',',B1{i});
% % % % % %             i1 = 1;     % Initial input
% % % % % %             if isempty(Find_str_B1)
% % % % % %                 B2(i,i1) = B1(i);
% % % % % %             else
% % % % % %                 int_txt = strsplit(B1{i},',');
% % % % % %                 B2(i,1:length(int_txt)) = int_txt;
% % % % % %             end
% % % % % %         end
% % % % % %         if i0==1
% % % % % %             category_data = B2(28,:);
% % % % % %         end
% % % % % %         B2data(:,:,i0) = str2num(cell2mat(B2(29:end,1:end)));

        B2data(:,:,i0) = B2;
        file_name_int = files2load2results2{i0};
        file_name_int2 = file_name_int(26:end);
        if ~isempty(strfind(file_name_int,'_Base_'))    % Reduce first element in the baseload runs (as a result of model operation)
            B2data(1,2,i0)=B2data(1,2,i0)/5;
            B2data(1,15,i0)=B2data(1,15,i0)/5;
        end
%         B2datacell(:,:,i0) = num2cell(B2data(:,:,i0));  % Convert to cell for combination with other properties and export
        B2datacell((i0-1)*8760+1:i0*8760,:) = num2cell(B2data(:,:,i0));  % Convert to cell for combination with other properties and export
        
        break_title = findstr(file_name_int2,'_'); 
        Add_categories = {file_name_int2(break_title(1)+1:break_title(2)-1),...
                          file_name_int2(break_title(2)+1:break_title(3)-1),...
                          file_name_int2(break_title(3)+1:break_title(4)-1),...
                          file_name_int2(break_title(4)+1:break_title(5)-1),...
                          file_name_int2(break_title(5)+1:break_title(6)-1),...
                          file_name_int2(break_title(6)+1:break_title(7)-3),...
                          file_name_int2(break_title(7)+1:break_title(8)-3),...
                          file_name_int2(break_title(8)+1:end-7)}; 
%         B2data_add(:,:,i0) = repmat(Add_categories,8760,1);
        B2data_add((i0-1)*8760+1:i0*8760,:) = repmat(Add_categories,8760,1);
            

% % % % %         %%% Write all values to file
% % % % %         xlswrite([path_name2,'Combined_operation_data_matlab.xlsx'],[B2datacell(:,:,i0),B2data_add(:,:,i0)],'sheet1',['A',num2str(2+(i0-1)*8760)]);        
        time1 = toc;
        disp(['Results ',num2str(i0),' of ',num2str(length(files2load2results2)),' :: ',num2str(round((time1/i0*length(files2load2results2)-time1)/60,1)),' minutes remaining'])
    end
    clear B2data
%     B2datacell_long = reshape(B2datacell,[],length(B2data_header));
%     B2data_add_long = reshape(B2data_add,[],length(header_data3));

    
    %%% Print results to file
    [m0 n0] = size(B2datacell);
    [m00 n00] = size(B2data_add);
    
%     C2 = [B2data_header,header_data3];
%     xlswrite([path_name2,'Combined_operation_data_matlab.xlsx'],C2);
    
    C2 = [[B2data_header,header_data3];[B2datacell,B2data_add]];
    clear B2data B2data_add B2datacell

%     xlswrite([path_name2,'Combined_operation_data_matlab.xlsx'],C2,'sheet1','A2');
    [m1 n1] = size(C2);
    div_files = 10;
    for i2=1:div_files
        if i2==1 div_file_start(i2)=1;
        else     div_file_start(i2)=round(m1*(i2-1)/div_files,0);
        end
        if i2>1 div_file_end(i2-1)=div_file_start(i2)+1;
        end
    end
    div_file_end(div_files)=m1;  
    
    for i3=1:div_files
        start1 = div_file_start(i3);   %((i3-1)*m1/div_files+1);
        end1 = div_file_end(i3);     %(i3*m1/div_files);
        xlswrite([path_name2,'Combined_operation_data_matlab',num2str(i3),'.xlsx'],C2(start1:end1,:)); %%%,'sheet1',['A1:',char(n1+'A'-1),num2str(m1)]);        
    end
    
% % % % %     [m1 n1] = size(C2);
% % % % %     fid = fopen([path_name2,'Combined_operation_data_matlab',num2str(i0),'.csv'],'wt'); % Write output to file
% % % % %     s1 = fgetl(fid);
% % % % %     for i1=1:1000 %m1
% % % % %         fprintf(fid,[repmat(['%s,'],1,n1-1),'%s'],C2{i1,:});
% % % % %         fprintf(fid,'\n');
% % % % %         if mod(i1,100000)==0
% % % % %             disp([num2str(i1),' of ',num2str(m1)])
% % % % %         end
% % % % %     end
% % % % %     fclose(fid);
end


%% Load Input
if 1==0
    B2input = zeros(8760,11,length(files2load2inputs)); % Initialize input matrix (row=time, col=property, depth=scenario)
    B2input_header= {'Interval','Elec Purchase ($/MWh)','Elec Sale ($/MWh)','Reg Up ($/MW)','Reg Dn ($/MW)','Spin Res ($/MW)','Nospin Res ($/MW)','Nat Gas ($/MMBTU)','H2 ($/kg)','Renewable In (MW)','Meter ($/mth)'};
    for i0=1:length(files2load2inputs)
        fid = fopen([path_name2,char(files2load2inputs(i0))],'rt'); % Read and parse results file
        A1 = fscanf(fid,'%c');                  %Read data and store in A
        fclose(fid);                            %Close file
        A1_line_break = regexp(A1,'\n');        %Determine where line breaks are (Removed \r\n for HPC runs)
        B1 = cell(length(A1_line_break),1);     %Initialize matrix to hold text file separated by lines
        B2 = cell(length(A1_line_break),9);     %Initialize matrix to hold text file separated by lines
        for i=1:length(A1_line_break)           %Separate text file by lines and save in B
            if i==1, B1(i) = cellstr(A1(1:A1_line_break(i)));
            else     B1(i) = cellstr(A1(A1_line_break(i-1):A1_line_break(i)));
            end
            Find_str_B1 = findstr(',',B1{i});
            i1 = 1;     % Initial input
            if isempty(Find_str_B1)
                B2(i,i1) = B1(i);
            else
                int_txt = strsplit(B1{i},',');
                B2(i,1:length(int_txt)) = int_txt;
            end
        end    
        if i0==1
            category_inputs = B2(36,:);
        end
        B2input(:,:,i0) = str2num(cell2mat(B2(37:end,1:end)));
        disp(['Inputs ',num2str(i0),' of ',num2str(length(files2load2inputs))])
    end
    
    % Average between time periods for TOU rates (DO NOT USE "B2data2" for CAISO rates)
    if 1==0
        B2data2 = B2data;
        for i1=1:length(files2load2inputs)
            c0=1; 
            for i2=2:m1        
                if (B2input(i2,13,i1)-B2input(i2-1,13,i1)~=0);
                    B2data2(c0:i2-1,2,i1)= mean(B2data(c0:i2-1,2,i1));
                    c0=i2;
                end
            end
        end
    end
end

%% Quick load input files (NODES)
if 1==0
    Select_range = 1; % 0=record all columns, 1=record only one column
    Select_column = 2; % Select column to record (select only one if Select_range=1)
    B2input_header= {'Interval','Elec Purchase ($/MWh)','Elec Sale ($/MWh)','Reg Up ($/MW)','Reg Dn ($/MW)','Spin Res ($/MW)','Nospin Res ($/MW)','Nat Gas ($/MMBTU)','H2 ($/kg)','Renewable In (MW)','Meter ($/mth)'};
    if Select_range==0
        B2input = zeros(8760,11,length(files2load2inputs)); % Initialize input matrix (row=time, col=property, depth=scenario)        
    elseif Select_range==1
        B2input = zeros(8760,length(files2load2inputs)); % Initialize input matrix (row=time, col=property, depth=scenario)
        B2input_header = B2input_header(Select_column); % Shorten header to singel column
    end    
    for i0=1:length(files2load2inputs)
        if Select_range==0
            B2input(:,:,i0) = csvread([path_name2,char(files2load2inputs(i0))],36,1);   % Load entire input file
        elseif Select_range==1
            B2input(:,i0) = csvread([path_name2,char(files2load2inputs(i0))],36,Select_column-1,[36 Select_column-1 35+8760 Select_column-1]);  % Load a specific row
        end
        disp(['Inputs ',num2str(i0),' of ',num2str(length(files2load2inputs))])
    end    
    Elec_price = reshape(mean(B2input(:,:),1),[],1);
    header_data2 = {};
    for i1=1:length(files2load2inputs)
        int2 = files2load2inputs{i1};        
        int3 = int2(25:end-4);
        break_title = findstr(int3,'_'); 
        if 1==1     % NODES: Change how additional categories are separated
            header3 = {'Scenario','Property','Value','Operation','Utility','Rate','Connection','Renewable Type','Capacity Factor','Installed renewables (% of cap)','Storage duration','Market','Node Name'};    
            Add_categories = {int3(break_title(1)+1:break_title(2)-1),...
                              int3(break_title(2)+1:break_title(3)-1),...
                              int3(break_title(3)+1:break_title(4)-1),...
                              int3(break_title(5)+1:break_title(6)-1),...
                              int3(break_title(end-3)+1:break_title(end-2)-1),...
                              int3(break_title(end-2)+1:break_title(end-1)-3),...
                              int3(break_title(end-1)+1:break_title(end)-3),...
                              int3(break_title(end)+1:end-3),...
                              int3(break_title(7)+1:break_title(8)-1),...
                              int3(break_title(8)+1:break_title(end-3)-1)};
        else
            if length(break_title)>8
                break_title(5)=[];
            end
            header_data3 = {'Scenario','Property','Value','Operation','Utility','Rate','Connection','Renewable Type','Capacity Factor','Installed renewables (% of cap)','Storage duration'};
            Add_categories = {int1(break_title(1)+1:break_title(2)-1),...
                              int1(break_title(2)+1:break_title(3)-1),...
                              int1(break_title(3)+1:break_title(4)-1),...
                              int1(break_title(4)+1:break_title(5)-1),...
                              int1(break_title(5)+1:break_title(6)-1),...
                              int1(break_title(6)+1:break_title(7)-3),...
                              int1(break_title(7)+1:break_title(8)-3),...
                              int1(break_title(8)+1:end-3)};                          
        end 
        header_data2(i1,:)=Add_categories;
    end
    
    % Create name values for nodes for use in Excel and GAMS
    util_header = {'PGAE','SCE','SDGE'};
    node_header = {'Load','Generator','Hub','Aggregate','Zone'};
    file_ind = [];
    file_data = cell(length(files2load2),19);    
    c4=0;
    % Load file describing the nodal data to be analyzed
    [num9,txt9,raw9] = xlsread('GAMS_nodes_caliso_v4_data');
    for i6=1:length(header_data2)
        for i7=1:length(raw9)   % Matches the node from the cross reference file to the .csv file data
            files_name_wo_csv = header_data2{i6,end};
%             files_name_wo_csv = files_name_wo_csv(1:length(files_name_wo_csv)-4);
            match_node = strcmp(raw9{i7,1},files_name_wo_csv);
            if ((match_node == 1) && (strcmp(raw9{i7,19},util_header{1}) || strcmp(raw9{i7,19},util_header{2}) || strcmp(raw9{i7,19},util_header{3})))
                % REMOVED FROM IF statement ABOVE  (strcmp(raw9{i7,2},'Aggregate')) &&             
                c4 = c4+1;
                file_ind(c4,1) = i7;    % All files that have a CSV file and are in IOU areas
                file_data(c4,:) = raw9(i7,:);
                util1(c4,1) = strcmp(raw9{i7,19},util_header{1})*1+...
                              strcmp(raw9{i7,19},util_header{2})*2+...
                              strcmp(raw9{i7,19},util_header{3})*3;                
                break
            end        
        end
        disp([num2str(i6),' of ',num2str(length(files2load2))])
    end
    lat_long = raw9(file_ind,[17,18]);
    Elec_price_cell = num2cell(Elec_price);
    data_output1 = [[header3(:,end),header3(:,4:end-1),{'Elec Prices','Latitude','Longitude'}];...
                    [header_data2(:,end),header_data2(:,1:end-1),Elec_price_cell,lat_long]];
    xlswrite([path_name2,'Elec_prices.xlsx'],data_output1);
end


%% Print output data to file
[m0 n0] = size(data2);
file_num = ceil(n0/16000);
for i0=1:file_num
    C2 = [category2,data2(:,((i0-1)*(n0/file_num))+1:i0*n0/file_num)];
    [m1 n1] = size(C2);
    fid = fopen([path_name2,'Combined_results',num2str(i0),'.csv'],'wt'); % Write output to file
    s1 = fgetl(fid);
    for i1=1:m1
        fprintf(fid,[repmat(['%s,'],1,n1-1),'%s'],C2{i1,:});
        fprintf(fid,'\n');
    end
    fclose(fid);
end

%% Process output file to include cost calculations then pivot data and export
try
    [num1,txt1,raw1] = xlsread([path_name2,'Combined_results_adj.xlsx']);
    
% %     %%%% Combine two files (doesn't work)
% %     [num0A,txt0A,raw0A] = xlsread([path_name2,'Combined_results_adj1.xlsx']);
% %     [num0B,txt0B,raw0B] = xlsread([path_name2,'Combined_results_adj2.xlsx']);    
% %     num1  = [num0A,num0B];
% %     txt1  = [txt0A,txt0B(:,2:end)];
% %     raw1  = [raw0A,raw0B(:,2:end)];
   
    category3 = raw1(:,1);
    data3 = raw1(:,2:end);
    [m1,n1] = size(raw1);
    
    %%% Pivot data
    c0=0;
    Select_rows = [5:19,37:61,66,73:77,80:84,87,88,91,92];
    data4 = {};
    header_data3 = {'Scenario','Property','Value','Operation','Utility','Rate','Connection','Renewable Type','Capacity Factor','Installed renewables (% of cap)','Storage duration'};
    names1 = txt1(1,2:end);
    for i1=1:n1-1
        int1 = names1{1,i1};
        break_title = findstr(int1,'_'); 
        if 1==Select_file     % NODES: Change how additional categories are separated
            header_data3 = {'Scenario','Property','Value','Operation','Utility','Rate','Connection','Renewable Type','Capacity Factor','Installed renewables (% of cap)','Storage duration','Market','Node Name'};    
            Add_categories = {int1(break_title(1)+1:break_title(2)-1),...
                              int1(break_title(2)+1:break_title(3)-1),...
                              int1(break_title(3)+1:break_title(4)-1),...
                              int1(break_title(5)+1:break_title(6)-1),...
                              int1(break_title(end-3)+1:break_title(end-2)-1),...
                              int1(break_title(end-2)+1:break_title(end-1)-3),...
                              int1(break_title(end-1)+1:break_title(end)-3),...
                              int1(break_title(end)+1:end-3),...
                              int1(break_title(7)+1:break_title(8)-1),...
                              int1(break_title(8)+1:break_title(end-3)-1)};
        else        % UTILITY: Change how additional categories are separated
            if strcmp(int1(1:3),'SMR')
                break_title = [0,break_title];
            end
            if length(break_title)>8
                break_title(5)=[];
            end
            header_data3 = {'Scenario','Property','Value','Operation','Utility','Rate','Connection','Renewable Type','Capacity Factor','Installed renewables (% of cap)','Storage duration'};
            Add_categories = {int1(break_title(1)+1:break_title(2)-1),...
                              int1(break_title(2)+1:break_title(3)-1),...
                              int1(break_title(3)+1:break_title(4)-1),...
                              int1(break_title(4)+1:break_title(5)-1),...
                              int1(break_title(5)+1:break_title(6)-1),...
                              int1(break_title(6)+1:break_title(7)-3),...
                              int1(break_title(7)+1:break_title(8)-3),...
                              int1(break_title(8)+1:end-3)};
        end
        c0=c0+length(Select_rows);
        data4(c0-length(Select_rows)+1:c0,:) = [repmat(names1(1,i1),length(Select_rows),1),category3(Select_rows,1),data3(Select_rows,i1),repmat(Add_categories,length(Select_rows),1)];
        if (mod(i1,100)==0)
            disp([num2str(i1),' of ',num2str(n1-1)]);
        end
    end
    if length(data4)<300000
        xlswrite([path_name2,'Pivoted_data.xlsx'],[header_data3;data4]);    
    elseif length(data4)<600000
        data4_len = round((length(data4))/2);
        xlswrite([path_name2,'Pivoted_data1.xlsx'],[header_data3;data4(1:data4_len,:)]);    
        xlswrite([path_name2,'Pivoted_data2.xlsx'],[header_data3;data4(length(data4)-data4_len+1:end,:)]);    
    elseif length(data4)<1000000
        xlswrite([path_name2,'Pivoted_data1.xlsx'],[header_data3;data4(1:300000,:)]);    
        xlswrite([path_name2,'Pivoted_data2.xlsx'],[header_data3;data4(300001:600000,:)]);    
        xlswrite([path_name2,'Pivoted_data3.xlsx'],[header_data3;data4(600001:length(data4),:)]);    
    else
        error('Make new pivot data files: Length of data4 is too long')
    end
catch
end


%% Adjust combined CSV file
if 1==0
    fid = fopen(['C:\Users\jeichman\Documents\gamsdir\projdir\Output\TEST\COMBINE_CSV_FILES\all_data1.csv']); % Read and parse file
    C11 = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',',','HeaderLines',0,'MultipleDelimsAsOne',true, 'CollectOutput',false);
    fclose(fid);

    fid = fopen(['C:\Users\jeichman\Documents\gamsdir\projdir\Output\TEST\COMBINE_CSV_FILES\all_data11.csv'], 'wt') ;
    data_int1 = C11{1,1};
    data_int2 = C11{1,2};
    data_int3 = C11{1,3};
    data_int4 = C11{1,4};
    data_int5 = C11{1,5};
    data_int6 = C11{1,6};
    data_int7 = C11{1,7};
    data_int8 = C11{1,8};
    data_int9 = C11{1,9};
    data_int10 = C11{1,10};
    data_int11 = C11{1,11};
    data_int12 = C11{1,12};
    data_int13 = C11{1,13};
    clear C11
    for i1=1:length(data_int2)
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',[data_int1{i1}],[',',data_int2{i1}],[',',data_int3{i1}],...
                                                                           [',',data_int4{i1}],[',',data_int5{i1}],[',',data_int6{i1}],...
                                                                           [',',data_int7{i1}],[',',data_int8{i1}],[',',data_int9{i1}],...
                                                                           [',',data_int10{i1}],[',',data_int11{i1}],[',',data_int12{i1}],...
                                                                           [',',data_int13{i1}]);
        if (mod(i1,1000)==0)
            disp([num2str(i1),' of ',num2str(length(data_int2))])
        end
    end
    fclose(fid);
end

%% Calculations and adjustments of results and input data
%%%%%%%
% Not necessary since demand charges were put into GAMS
%%%%%%%
if 1==0 
    Vec_length = length(B2input(:,1,1));                % Calculate vector length for use later
    
    TOU_schedule = reshape(B2input(:,13,:),[],12);      % Get TOU schedule
    TOU_schedule_diff = diff(TOU_schedule);             % Calculate changes in TOU schedule
    Renewable_gen = reshape(B2data(:,14,:),[],12);      % Initialize generation signal
    Non_renewable_gen = reshape(B2data(:,15,:),[],12);  % Initialize generation signal    
    Non_renewable_gen_smooth = Non_renewable_gen;       % Initialize generation signal to smooth
    
    set1=1;   % Start of set identifier
    set2=1;   % End of set identifier
    for i2=1:length(files2load2)
        for i1=1:length(TOU_schedule(:,i2))-1
            if(i1==1), 
                if(Non_renewable_gen(i1,i2)>=5*Non_renewable_gen(i1+1,i2)),Non_renewable_gen(i1,i2)=Non_renewable_gen(i1,i2)/5; 
                end 
            end
            if(TOU_schedule_diff(i1,i2)==0),
            else
                set2=i1;
                Non_renewable_gen_smooth(set1:set2,i2)=mean(Non_renewable_gen(set1:set2,i2));
                set1=i1+1;
            end
        end
        set1=1;   % Reset set identifier
        set2=1;   % Reset set identifier
    end
    Total_gen_smooth = Non_renewable_gen_smooth+Renewable_gen;
    Total_gen_smooth_MW = Total_gen_smooth/1000;
    Non_renewable_gen_smooth_MW = Non_renewable_gen_smooth/1000;
    Renewable_gen_MW = Renewable_gen/1000;
    
    %Calculate Demand Charges
    Monthly_charges = zeros(12,length(files2load2),3);   % Initialize vector of monthly demand charges (meter, fixed demand, timed demand)
    Fix_dem_chg = reshape(B2input(:,11,:),[],length(files2load2));     % Draw out fixed demand charge
    Timed_dem_chg = reshape(B2input(:,12,:),[],length(files2load2));   % Draw out timed demand charge
        
    for i2=1:length(files2load2)
       
        for i1=1:12     % Monthly calculations
            % Meter charge
            if findstr('PGE',files2load2{i2})                
                 Monthly_charges(i1,i2,1) = B2input(1,14,i2)*eomday(year_val,i1);    % Calculate the monthly charge ($) (if statement used because PGE prescribes over month)
            else Monthly_charges(i1,i2,1) = B2input(1,14,i2);                        % Calculate the monthly charge ($)
            end
            % Fixed Demand Charge
            max_val = max(Non_renewable_gen_smooth_MW((num1(:,3)==i1),i2));                % Find max monthly value
            max_loc = find(Non_renewable_gen_smooth_MW((num1(:,3)==i1),i2) == max_val);    % Find location(s) of max value
            Monthly_charges(i1,i2,2) = max(Fix_dem_chg(max_loc,i2))*1000*max_val;             % Calculate the monthly charge ($)
            % Timed Demand Charge
            hrs_on_peak = (num1(:,3)==i1).*(TOU_schedule(:,i2)==1);
            hrs_mid_peak = (num1(:,3)==i1).*(TOU_schedule(:,i2)==2);
            hrs_off_peak = (num1(:,3)==i1).*(TOU_schedule(:,i2)==3);
            
            max_on_peak = max(Non_renewable_gen_smooth_MW(find(hrs_on_peak==1),i2));    % Find max monthly value
            max_mid_peak = max(Non_renewable_gen_smooth_MW(find(hrs_mid_peak==1),i2));   % Find max monthly value
            max_off_peak = max(Non_renewable_gen_smooth_MW(find(hrs_off_peak==1),i2));   % Find max monthly value
            
            if (find(hrs_on_peak>0))    % Calculates indices of max values
                 max_loc_on_peak = find((Non_renewable_gen_smooth_MW(:,i2) == max_on_peak).*hrs_on_peak);       % Find location(s) of max value
            else max_loc_on_peak = [];
            end
            if (find(hrs_mid_peak>0))
                 max_loc_mid_peak = find((Non_renewable_gen_smooth_MW(:,i2) == max_mid_peak).*hrs_mid_peak);	% Find location(s) of max value
            else max_loc_mid_peak = [];
            end
            if (find(hrs_off_peak>0))
                 max_loc_off_peak = find((Non_renewable_gen_smooth_MW(:,i2) == max_off_peak).*hrs_off_peak); 	% Find location(s) of max value
            else max_loc_off_peak = [];
            end
            Monthly_charges(i1,i2,3) = sum(max(Timed_dem_chg(max_loc_on_peak,i2).*Non_renewable_gen_smooth_MW(max_loc_on_peak,i2)))*1000+...
                                       sum(max(Timed_dem_chg(max_loc_mid_peak,i2).*Non_renewable_gen_smooth_MW(max_loc_mid_peak,i2)))*1000+...
                                       sum(max(Timed_dem_chg(max_loc_off_peak,i2).*Non_renewable_gen_smooth_MW(max_loc_off_peak,i2)))*1000;       % Calculate the monthly charge ($)            
        end
    end
    
    % Yearly charges (meter, fixed demand, timed demand)
    Yearly_charges = reshape(sum(Monthly_charges,1),[],3)';
    
end
%% Generate renewable H2 vector
if 1==0

clear all, close all, clc    

cd('C:\Users\jeichman\Documents\gamsdir\projdir');
C3=1;    %initialize counter
path_name2(C3) = {'Output\CAISO_0711_Base\'}; C3=C3+1;
path_name2(C3) = {'Output\CAISO_0711_RENx2\'}; C3=C3+1;
path_name2(C3) = {'Output\CAISO_0711_RENx3\'}; C3=C3+1;
path_name2(C3) = {'Output\CAISO_0711_RENx4\'}; C3=C3+1;
path_name2(C3) = {'Output\CAISO_0711_RENx5\'}; C3=C3+1;

RENPEN_hdr = {'2022','RENx2','RENx3','RENx4','RENx5'};
load('RENPEN_period') 
[M1,N1] = size(RENPEN_period);

for C4 = 1:length(path_name2)
%%% Loads all csv files in a given folder
    files2load = dir(path_name2{C4}); files2load2={files2load.name}';  % Identify files in a folder    
    for i0=1:length(files2load2) % Remove items from list that do not fit criteria
        load_file1(i0)=~isempty(strfind(files2load2{i0},'.csv')); 
        load_file1(i0)=~isempty(strfind(files2load2{i0},'dispatch_results')); 
    end 
    files2load2=files2load2(load_file1); clear load_file1
    
    if C4==1, Ren_portion_max=zeros(length(files2load2),N1);  % Sets up matrix for max renewable portion
              Ren_portion_min=zeros(length(files2load2),N1);  % Sets up matrix for min renewable portion
              Elec_gen_mat=zeros(length(files2load2),N1);  % Sets up matrix for min renewable portion
              Elec_consume_mat=zeros(length(files2load2),N1);  % Sets up matrix for min renewable portion
    end
    
    for i0=1:length(files2load2)
        i1=1;
        if findstr(char(files2load2(i0)),'Batt'), i1=0;
        elseif findstr(char(files2load2(i0)),'HYPS'), i1=0;
        elseif findstr(char(files2load2(i0)),'noH2sale'), i1=0;  
        end
        if i1==1
            [data1,text1,raw1] = xlsread([path_name2{C4},char(files2load2(i0))]);
            elec_consumed = data1(24:M1+23,2);  % Calculates the electricity that is consumed
            elec_gen = data1(24:M1+23,3);       % Calculates the electricity that is generated
            Elec_consume_mat(i0,C4) = sum(elec_consumed);
            Elec_gen_mat(i0,C4) = sum(elec_gen);
            Ren_portion_max(i0,C4) = sum(elec_consumed.*RENPEN_period(:,C4))/sum(elec_consumed);  % Calculates the renewable penetration if all 
            Ren_portion_min(i0,C4) = (sum(elec_consumed.*RENPEN_period(:,C4))-sum(elec_gen))/sum(elec_consumed);
        end
        disp(['Completed ',num2str(i0+(length(files2load2)*(C4-1))),' of ',num2str(length(files2load2)*length(path_name2))])
    end     
end      
end


