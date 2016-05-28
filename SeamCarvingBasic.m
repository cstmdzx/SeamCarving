% SeamCurving Algorithm
% authored by nklyp
clc;clear;
%读入图片
[fn,pn,fi]=uigetfile('*.jpg','选择图片');
I=imread([pn fn]);
tic;
%I = imrotate(I,90);

%--------------循环的次数即是删除的Seam的个数-----------------
for cn = 1:100
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

    %计算偏导数
    % for r = 1:size_r
    %     for c = 1:size_c
    %         %此处还未做特殊处理,第一行第一列最后一列均做特殊处理
    %         if c == 1
    %             energy(r,c) = abs(I(r,c,1)-I(r,c+1,1))+abs(I(r,c,2)-I(r,c+1,2))+abs(I(r,c,3)-I(r,c+1,3));%x方向偏导数，第一列特殊处理
    %         else
    %             energy(r,c) = abs(I(r,c,1)-I(r,c-1,1))+abs(I(r,c,2)-I(r,c-1,2))+abs(I(r,c,3)-I(r,c-1,3));%x方向偏导数
    %         end
    %         
    %         if r == 1
    %             energy(r,c) = energy(r,c)+abs(I(r,c,1)-I(r+1,c,1))+abs(I(r,c,2)-I(r+1,c,2))+abs(I(r,c,3)-I(r+1,c,3));%y方向偏导数，第一行特殊处理
    %         else
    %             energy(r,c) = energy(r,c)+abs(I(r,c,1)-I(r-1,c,1))+abs(I(r,c,2)-I(r-1,c,2))+abs(I(r,c,3)-I(r-1,c,3));%y方向偏导数
    %         end
    %         
    % %         energy(r,c)=min( energy(r-1,c-1), min( energy(r-1,c),energy(r-1,c+1) ) )  + energy(r,c);%之前的最优解
    % %         
    % %         if energy(r-1,c-1) > energy(r-1,c+1) && energy(r-1,c) > energy(r-1,c+1)%保存路径
    % %             route(r) = c+1;
    % %         else if energy(r-1,c-1) > energy(r-1,c)
    % %                 route(r) = c;
    % %             else
    % %                 route(r) = c-1;
    % %             end
    % %         end
    %         
    %         
    %     end
    % end

    %计算energy,同时保存路径
    e_col_f = energy (:,2);
    e_col_l = energy (:,size_c-1);
    
    energy1 = [e_col_f energy(:,1:size_c-1)];
    energy3 = [energy(:,2:size_c) e_col_l];
    
    energy_tri(:,:,1) = energy1;
    energy_tri(:,:,2) = energy;
    energy_tri(:,:,3) = energy3;



%    tic;
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
    
%     for r = 2:size_r
%         for c = 1:size_c
%             if c == 1
%                 energy(r,c)=min( energy(r-1,c),energy(r-1,c+1) )+ energy(r,c);
%                 if energy(r-1,c) > energy(r-1,c+1)
%                     route(r,c) = c+1;
%                 else
%                     route(r,c) = c;
%                 end
%             end
%             if c == size_c
%                 energy(r,c)=min( energy(r-1,c),energy(r-1,c-1) )+ energy(r,c);
%                 if energy(r-1,c) > energy(r-1,c-1)
%                     route(r,c) = c-1;
%                 else
%                     route(r,c) = c;
%                 end
%             end
% 
%             if c > 1 && c < size_c
%                 energy(r,c)=min( energy(r-1,c-1), min( energy(r-1,c),energy(r-1,c+1) ) )  + energy(r,c);
%                 if energy(r-1,c-1) > energy(r-1,c+1) && energy(r-1,c) > energy(r-1,c+1)
%                     route(r,c) = c+1;
%                 else if energy(r-1,c-1) > energy(r-1,c)
%                         route(r,c) = c;
%                     else
%                         route(r,c) = c-1;
%                     end
%                 end 
%             end
% 
%         end
%     end
 %   toc;
    %遍历最后一行，找到最小的energy,同时定位
%     mie = energy(size_r,1);
%     col = 1;
%     for c = 1:size_c
%         if energy(size_r,c) < mie
%             mie = energy(size_r,c);
%             col = c;
%         end
%     end
    
    [energy_min index_min] = min (energy(size_r,:));

    %删除掉energy最小的那一列
%     res = [];
%     for r = size_r:-1:1
%         if col > 1
%             temp = I(r,1:col-1,:);
%         end
%         if col < size_c
%             temp = [temp I(r,col+1:size_c,:)];
%         end
% 
%         res = [temp;res];
%         %I(r,col,:) = [255,255,255];
%         col = route(r,col);%更新col
%         temp = [];
%     end
%     I = res;
    
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
%I = imrotate(I,-90);
imshow(I);