% SeamCurving Algorithm
% authored by nklyp
clc;clear;
%����ͼƬ
[fn,pn,fi]=uigetfile('*.jpg','ѡ��ͼƬ');
I=imread([pn fn]);
tic;
cn = 0;
while(true)
    [size_r,size_c,n] = size(I);
    route = zeros(size_r,size_c);
    
    %ȡ��ͼ���rgb���������ҵ�����ֵ��index
    fr = I(:,:,1);
    fg = I(:,:,2);
    fb = I(:,:,3);
    
    indexR = find(fr > 240);
    indexG = find(fg < 20);
    indexB = find(fb < 20);
    
    indexRGB = intersect( intersect(indexR,indexG),indexB );
    
    if(isempty(indexRGB))
        break;
    end

    cn = cn+1;
    
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
    energy = sum(I_sum,3);
    
    %������ֵ��energy��Ϊ��С
    energy(indexRGB) = -1000;
    
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
    

    [energy_min index_min] = min (energy(size_r,:));

    %ɾ����energy��С����һ��    
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
    indexR = [];
    indexG = [];
    indexB = [];
    indexRGB = [];
%     if mod(cn,50) == 0
%         imshow(I);
%     end
end
toc;
imshow(I);