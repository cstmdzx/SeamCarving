% SeamCurving Algorithm
% authored by nklyp
clc;clear;

%设置要删除的行数和列数
RowNum = 101;
ColNum = 101;
entryT = zeros(RowNum,ColNum);
map = zeros(RowNum,ColNum);

IMAGE = cell(1,ColNum);

%读入图片
[fn,pn,fi]=uigetfile('*.jpg','选择图片');
I=imread([pn fn]);
tic;
IMAGE{1} = I;

%--------------先处理第一行-------------
for cn = 2:ColNum
    [size_r,size_c,n] = size(I);
    route = zeros(size_r,size_c);

    %取一个temp用来处理边界
    temp_row = I(2,:,:);
    temp_col = I(:,2,:);

    %求x方向偏导数
    I_row1 = [I;temp_row];
    I_row2 = [temp_row;I];
    I_row = abs(I_row1 - I_row2);
    I_rowr = I_row(1:size_r,:,:);

    %求y方向偏导数
    I_col1 = [I temp_col];
    I_col2 = [temp_col I];
    I_col = abs(I_col1 - I_col2);
    I_colr = I_col(:,1:size_c,:);

    %求和
    I_sum = I_rowr + I_colr;

    energy = sum(I_sum,3);

    %计算energy,同时保存路径
    e_col_f = energy (:,2);
    e_col_l = energy (:,size_c-1);
    
    energy1 = [e_col_f energy(:,1:size_c-1)];
    energy3 = [energy(:,2:size_c) e_col_l];
    
    energy_tri(:,:,1) = energy1;
    energy_tri(:,:,2) = energy;
    energy_tri(:,:,3) = energy3;

    for r = 2:size_r
        row_min(:,:,:) = energy_tri(r-1,:,:);
        [row2 index] = min(row_min,[],3);
        row2 = row2+energy (r,:);
        energy (r,:) = row2;
        row2_col_f = row2(2);
        row2_col_l = row2(size_c-2);
        energy_tri(r,:,1) = [row2_col_f row2(1:size_c-1)];
        energy_tri(r,:,2) = row2;
        energy_tri(r,:,3) = [row2(2:size_c) row2_col_l];
        route(r,:) = index + [1:size_c] - 2;
    end
    route(route == 0) = 2;
    route(route == size_c+1) = size_c-1;
    

    %遍历最后一行，找到最小的energy,同时定位   
    [energy_min index_min] = min (energy(size_r,:));
    entryT(1,cn) = entryT(1,cn-1) + energy_min;
    
    %删除掉energy最小的那一列
    for r = size_r:-1:1
       % temp = I(r,index_min,:);
        if(index_min < size_c)
            I(r,index_min:size_c-1,:) =I(r,index_min+1:size_c,:); 
        end
%        I(r,index_min,:) = [255 255 255];
        index_min = route(r,index_min);
    end
    I(:,size_c,:) = [];
    energy_tri = [];
    row_min = [];
    IMAGE{cn} = I;
end

%--------------循环的次数即是删除的Seam的个数-----------------
for cnx = 2:RowNum
    map(cnx,1) = 1;
    
    %处理每行第一个
    I = IMAGE{1};
    I = imrotate(I,90);
    
    [size_r,size_c,n] = size(I);
    route = zeros(size_r,size_c);

    %取一个temp用来处理边界
    temp_row = I(2,:,:);
    temp_col = I(:,2,:);

    %求x方向偏导数
    I_row1 = [I;temp_row];
    I_row2 = [temp_row;I];
    I_row = abs(I_row1 - I_row2);
    I_rowr = I_row(1:size_r,:,:);

    %求y方向偏导数
    I_col1 = [I temp_col];
    I_col2 = [temp_col I];
    I_col = abs(I_col1 - I_col2);
    I_colr = I_col(:,1:size_c,:);

    %求和
    I_sum = I_rowr + I_colr;

    energy = sum(I_sum,3);

    %计算energy,同时保存路径
    e_col_f = energy (:,2);
    e_col_l = energy (:,size_c-1);
    
    energy1 = [e_col_f energy(:,1:size_c-1)];
    energy3 = [energy(:,2:size_c) e_col_l];
    
    energy_tri(:,:,1) = energy1;
    energy_tri(:,:,2) = energy;
    energy_tri(:,:,3) = energy3;

    for r = 2:size_r
        row_min(:,:,:) = energy_tri(r-1,:,:);
        [row2 index] = min(row_min,[],3);
        row2 = row2+energy (r,:);
        energy (r,:) = row2;
        row2_col_f = row2(2);
        row2_col_l = row2(size_c-2);
        energy_tri(r,:,1) = [row2_col_f row2(1:size_c-1)];
        energy_tri(r,:,2) = row2;
        energy_tri(r,:,3) = [row2(2:size_c) row2_col_l];
        route(r,:) = index + [1:size_c] - 2;
    end
    route(route == 0) = 2;
    route(route == size_c+1) = size_c-1;
    

    %遍历最后一行，找到最小的energy,同时定位   
    [energy_min index_min] = min (energy(size_r,:));
    entryT(cnx,1) = entryT(cnx,1) + index_min;

    %删除掉energy最小的那一列
    for r = size_r:-1:1
       % temp = I(r,index_min,:);
        if(index_min < size_c)
            I(r,index_min:size_c-1,:) =I(r,index_min+1:size_c,:); 
        end
%        I(r,index_min,:) = [255 255 255];
        index_min = route(r,index_min);
    end
    I(:,size_c,:) = [];
    energy_tri = [];
    row_min = [];
    
    I = imrotate(I,270);
    
    IMAGE{1} = I;
    
    for cny = 2:ColNum
        
        I = IMAGE{cny-1};

        %-----------------------计算纵向seam的最低的energy值-----------------------
        
        [size_r_x,size_c_x,n_x] = size(I);

        ratio = size_r_x / size_c_x;

        route_x = zeros(size_r_x,size_c_x);

        %取一个temp用来处理边界
        temp_row_x = I(2,:,:);
        temp_col_x = I(:,2,:);

        %求x方向偏导数
        I_row1_x = [I;temp_row_x];
        I_row2_x = [temp_row_x;I];
        I_row_x = abs(I_row1_x - I_row2_x);
        I_rowr_x = I_row_x(1:size_r_x,:,:);

        %求y方向偏导数
        I_col1_x = [I temp_col_x];
        I_col2_x = [temp_col_x I];
        I_col_x = abs(I_col1_x - I_col2_x);
        I_colr_x = I_col_x(:,1:size_c_x,:);

        %求和
        I_sum_x = I_rowr_x + I_colr_x;

        energy_x = sum(I_sum_x,3);

        %计算energy,同时保存路径
        e_col_f_x = energy_x (:,2);
        e_col_l_x = energy_x (:,size_c_x-1);

        energy1_x = [e_col_f_x energy_x(:,1:size_c_x-1)];
        energy3_x = [energy_x(:,2:size_c_x) e_col_l_x];

        energy_tri_x(:,:,1) = energy1_x;
        energy_tri_x(:,:,2) = energy_x;
        energy_tri_x(:,:,3) = energy3_x;

        for r_x = 2:size_r_x
            row_min_x(:,:,:) = energy_tri_x(r_x-1,:,:);
            [row2_x index_x] = min(row_min_x,[],3);
            row2_x = row2_x+energy_x (r_x,:);
            energy_x (r_x,:) = row2_x;
            row2_col_f_x = row2_x(2);
            row2_col_l_x = row2_x(size_c_x-2);
            energy_tri_x(r_x,:,1) = [row2_col_f_x row2_x(1:size_c_x-1)];
            energy_tri_x(r_x,:,2) = row2_x;
            energy_tri_x(r_x,:,3) = [row2_x(2:size_c_x) row2_col_l_x];
            route_x(r_x,:) = index_x + [1:size_c_x] - 2;
        end
        route_x(route_x == 0) = 2;
        route_x(route_x == size_c_x+1) = size_c_x-1;

        %遍历最后一行，找到最小的energy,同时定位

        [energy_min_x index_min_x] = min (energy_x(size_r_x,:));

%         %把计算出的能量值乘一个权重
% 
%         energy_min_x = energy_min_x * ratio;

         %-----------------------计算横向seam的最低的energy值-----------------------

         I = IMAGE{cny};
         I = imrotate(I,90);

        [size_r_y,size_c_y,n_y] = size(I);

        route_y = zeros(size_r_y,size_c_y);

        %取一个temp用来处理边界
        temp_row_y = I(2,:,:);
        temp_col_y = I(:,2,:);

        %求x方向偏导数
        I_row1_y = [I;temp_row_y];
        I_row2_y = [temp_row_y;I];
        I_row_y = abs(I_row1_y - I_row2_y);
        I_rowr_y = I_row_y(1:size_r_y,:,:);

        %求y方向偏导数
        I_col1_y = [I temp_col_y];
        I_col2_y = [temp_col_y I];
        I_col_y = abs(I_col1_y - I_col2_y);
        I_colr_y = I_col_y(:,1:size_c_y,:);

        %求和
        I_sum_y = I_rowr_y + I_colr_y;

        energy_y = sum(I_sum_y,3);

        %计算energy,同时保存路径
        e_col_f_y = energy_y (:,2);
        e_col_l_y = energy_y (:,size_c_y-1);

        energy1_y = [e_col_f_y energy_y(:,1:size_c_y-1)];
        energy3_y = [energy_y(:,2:size_c_y) e_col_l_y];

        energy_tri_y(:,:,1) = energy1_y;
        energy_tri_y(:,:,2) = energy_y;
        energy_tri_y(:,:,3) = energy3_y;

        for r_y = 2:size_r_y
            row_min_y(:,:,:) = energy_tri_y(r_y-1,:,:);
            [row2_y index_y] = min(row_min_y,[],3);
            row2_y = row2_y+energy_y (r_y,:);
            energy_y (r_y,:) = row2_y;
            row2_col_f_y = row2_y(2);
            row2_col_l_y = row2_y(size_c_y-2);
            energy_tri_y(r_y,:,1) = [row2_col_f_y row2_y(1:size_c_y-1)];
            energy_tri_y(r_y,:,2) = row2_y;
            energy_tri_y(r_y,:,3) = [row2_y(2:size_c_y) row2_col_l_y];
            route_y(r_y,:) = index_y + [1:size_c_y] - 2;
        end
        route_y(route_y == 0) = 2;
        route_y(route_y == size_c_y+1) = size_c_y-1;

        %遍历最后一行，找到最小的energy,同时定位

        [energy_min_y index_min_y] = min (energy_y(size_r_y,:));


        %删除掉energy最小的
        value_x = entryT(cnx,cny-1) + energy_min_x;%纵向Seam
        value_y = entryT(cnx-1,cny) + energy_min_y;%横向Seam
 %      value_x = value_x * ratio;        
        if value_y > value_x
            entryT(cnx,cny) = value_x;
            I = IMAGE{cny-1};
            index_min = index_min_x;
            size_r = size_r_x;
            size_c = size_c_x;
            route = route_x;
            map(cnx,cny) = 0;
        else
            entryT(cnx,cny) = value_y;
            index_min = index_min_y;
            size_r = size_r_y;
            size_c = size_c_y;
            route = route_y;
            map(cnx,cny) = 1;
        end

        for r = size_r:-1:1
           % temp = I(r,index_min,:);
            if(index_min < size_c)
                I(r,index_min:size_c-1,:) =I(r,index_min+1:size_c,:); 
            end
    %        I(r,index_min,:) = [255 255 255];
            index_min = route(r,index_min);
        end
        I(:,size_c,:) = [];
        energy_tri_x = [];
        energy_tri_y = [];
        row_min_x = [];
        row_min_y = [];
        if ~(value_y > value_x)
            I = imrotate(I,270);
        end
        IMAGE{cny} = I;
    end
end
toc;
imshow(I);