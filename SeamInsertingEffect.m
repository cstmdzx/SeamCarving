% SeamCurving Algorithm
% authored by nklyp
clc;clear;
%����ͼƬ
[fn,pn,fi]=uigetfile('*.jpg','ѡ��ͼƬ');
I=imread([pn fn]);
tic;
%����֮ǰɾȥ��������Ŀ
cn = 120;

%����energy����
[size_r,size_c,n] = size(I);
route = zeros(size_r,size_c);
routeInc = zeros(size_r,size_c);

%ȡһ��temp��������߽�
temp_row = I(2,:,:);
temp_col = I(:,2,:);

%��x����ƫ����
I_row1 = [I;temp_row];
I_row2 = [temp_row;I];
I_row = abs(I_row1 - I_row2);
I_rowr = I_row(1:size_r,:,:);

%��y����ƫ����
I_col1 = [I temp_col];
I_col2 = [temp_col I];
I_col = abs(I_col1 - I_col2);
I_colr = I_col(:,1:size_c,:);

%���
I_sum = I_rowr + I_colr;
I_sum = sum(I_sum,3);
energy = I_sum;
   
%����energy,ͬʱ����·��
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

%����ߴ�
fr = [I(:,:,1) zeros(size_r,cn)];
fg = [I(:,:,2) zeros(size_r,cn)];
fb = [I(:,:,3) zeros(size_r,cn)];
I = [];
I2(:,:,1) = fr;
I2(:,:,2) = fg;
I2(:,:,3) = fb;
%imshow(I2);
I = I2;

%ֻ��Ҫ�ҳ���С��n����������Ҫ֪���ĸ������С��������Ҫ��insert
% cmp = sort( energy(size_r,:) );
% cmpCn = cmp(cn) + 1;%��cn+1����ֻҪ����С�ľ���ǰcn���е�һ��
% energy_row_l = energy(size_r,:);
%index = find(energy_row_l < cmpCn);


for i = 1:cn
    [useless_min flag] = min( energy(size_r,:) );%ÿһ�������С֮��energy��û���ˣ���Ϊ��route
    indexI = find(route(r,:) > flag);
%    index(indexI) = index(indexI) + 1;

    
    for r = size_r:-1:1
        
        offsetFlag = routeInc(r,flag)+flag;
        if flag == 1
            offsetFlag = offsetFlag + 1;
        end
        I_sum(r,flag) = 5000;%����һ���ϴ��energyֵ
        I(r,(1+offsetFlag):(size_c+i),:) = I(r,offsetFlag:(size_c+i-1),:);
%        I(r,offsetFlag,:) = floor( I(r,offsetFlag-1,:) /2 + I(r,offsetFlag,:) /2 );

        I(r,offsetFlag,1) = 255;
        I(r,offsetFlag,2) = 0;
        I(r,offsetFlag,3) = 0;
         %���Ÿ���route,index��flag
        
        if r < size_r && r > 1
            routeI = find(route(r,:) + routeInc(r-1,:) > offsetFlag);
            routeInc(r,routeI) = routeInc(r,routeI) + 1;
        end
%         if flag > size_c
%             flag = size_c;
%         end
        if r > 1
            flag = route( r,flag );
        end    
    end
    energy(:,:) = I_sum;%��energy������Ϊƫ�����ĺͣ���ʱ�����ٳ��ֵĻᱻ��Ϊ�ܸߵ�Ȩֵ
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
    routeInc(size_r,indexI) = routeInc(size_r,indexI) + 1;
%    energy(:,:) = energy_tri(:,:,2);
end
toc;
imshow(I);