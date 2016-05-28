% SeamCurving Algorithm
% authored by nklyp
clc;clear;

%设置要删除的行数和列数
cnx = 100;
cny = 50;
entryT = zeros(cnx,cny);
direction = zeros(cnx+cny);

%读入图片
[fn,pn,fi]=uigetfile('*.jpg','选择图片');
I=imread([pn fn]);
tic;

%--------------循环的次数即是删除的Seam的个数-----------------
for cn = 1:(cnx+cny)
    
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
    
    %把计算出的能量值乘一个权重
    
    energy_min_x = energy_min_x * ratio;
    
     %-----------------------计算横向seam的最低的energy值-----------------------
     
     I = rot90(I);
     
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

    %删除掉energy最小的那一列
    
    if energy_min_y > energy_min_x
        
    
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
end
toc;
imshow(I);