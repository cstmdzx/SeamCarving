% SeamCurving Algorithm
% authored by nklyp
clc;clear;

%����Ҫɾ��������������
cnx = 100;
cny = 50;
entryT = zeros(cnx,cny);
direction = zeros(cnx+cny);

%����ͼƬ
[fn,pn,fi]=uigetfile('*.jpg','ѡ��ͼƬ');
I=imread([pn fn]);
tic;

%--------------ѭ���Ĵ�������ɾ����Seam�ĸ���-----------------
for cn = 1:(cnx+cny)
    
    %-----------------------��������seam����͵�energyֵ-----------------------
    
    [size_r_x,size_c_x,n_x] = size(I);
    
    ratio = size_r_x / size_c_x;
    
    route_x = zeros(size_r_x,size_c_x);

    %ȡһ��temp��������߽�
    temp_row_x = I(2,:,:);
    temp_col_x = I(:,2,:);

    %��x����ƫ����
    I_row1_x = [I;temp_row_x];
    I_row2_x = [temp_row_x;I];
    I_row_x = abs(I_row1_x - I_row2_x);
    I_rowr_x = I_row_x(1:size_r_x,:,:);

    %��y����ƫ����
    I_col1_x = [I temp_col_x];
    I_col2_x = [temp_col_x I];
    I_col_x = abs(I_col1_x - I_col2_x);
    I_colr_x = I_col_x(:,1:size_c_x,:);

    %���
    I_sum_x = I_rowr_x + I_colr_x;

    energy_x = sum(I_sum_x,3);

    %����energy,ͬʱ����·��
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
    
    %�������һ�У��ҵ���С��energy,ͬʱ��λ
   
    [energy_min_x index_min_x] = min (energy_x(size_r_x,:));
    
    %�Ѽ����������ֵ��һ��Ȩ��
    
    energy_min_x = energy_min_x * ratio;
    
     %-----------------------�������seam����͵�energyֵ-----------------------
     
     I = rot90(I);
     
    [size_r_y,size_c_y,n_y] = size(I);
    
    route_y = zeros(size_r_y,size_c_y);

    %ȡһ��temp��������߽�
    temp_row_y = I(2,:,:);
    temp_col_y = I(:,2,:);

    %��x����ƫ����
    I_row1_y = [I;temp_row_y];
    I_row2_y = [temp_row_y;I];
    I_row_y = abs(I_row1_y - I_row2_y);
    I_rowr_y = I_row_y(1:size_r_y,:,:);

    %��y����ƫ����
    I_col1_y = [I temp_col_y];
    I_col2_y = [temp_col_y I];
    I_col_y = abs(I_col1_y - I_col2_y);
    I_colr_y = I_col_y(:,1:size_c_y,:);

    %���
    I_sum_y = I_rowr_y + I_colr_y;

    energy_y = sum(I_sum_y,3);

    %����energy,ͬʱ����·��
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
    
    %�������һ�У��ҵ���С��energy,ͬʱ��λ
   
    [energy_min_y index_min_y] = min (energy_y(size_r_y,:));

    %ɾ����energy��С����һ��
    
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