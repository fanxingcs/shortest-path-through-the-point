close,clear,clc
% 计算距离
sj = load('data.txt'); % 读取中间需要经过点的位置的x,y坐标
num = length(sj); % 中间经过位置的数量
d_start=[70,40]; % 起始位置坐标
d_end = [70,40]; % 终点位置坐标

% 组成坐标的矩阵，第一行是起始位置不会改变，最后一行是终点位置也不会改变
% 中间的坐标会根据算法调整位置，得到最短的消耗距离
sj=[d_start;sj;d_end];

% 利用坐标计算距离，由于这里是经纬坐标，所以使用了以下的方法
% 如果是平面坐标，可以根据自己的需要改成欧式距离等
% 最终得到距离矩阵d,即d(i,j)是第i个位置与第j个位置的距离
 sj=sj*pi/180;
d = zeros(length(sj));
for i=1:num+1
    for j=i+1:num+2
        temp=cos(sj(i,1)-sj(j,1))*cos(sj(i,2))*cos(sj(j,2))+sin(sj(i,2))*sin(sj(j,2));
        d(i,j)=6370*acos(temp);
    end
end
d=d+d';

%% 初始化工作
% 子代的数量
% 可以自定义，没有更好的选择可以先使用64 32等，为大于2的偶数
n = 64; 

% 初始化父代
gen = zeros(n,num+2);
gen(:,1) = 1;   gen(:,end)=num+2;
cost = zeros(n,1);
for geni=1:n
    gen(geni,2:end-1) = randperm(num)+1;
    for i = 1:num+1
        cost(geni) = cost(geni) + d(gen(geni,i),gen(geni,i+1));
    end
    [gen(geni,:),cost(geni)] = get_best(gen(geni,:),cost(geni),d);
end

%% 遗传核心部分
gen_best = gen(1,:);    cost_best = cost(1);
% count代表遗传多少代，一般1000多就可以了
% 多了会提高经度，但是也会增加程序运行时间
% 第一次运行不要取太高
for count = 1:100
    % 生成子代
    for geni = 1:floor(n/2)
        if geni<=n/6
            a = [gen(geni*2-1,:); gen(geni*2,:)];  a = unique(a(:));
            b = [gen(geni*2,:); gen(geni*2-1,:)];  b = unique(b(:));
            costa = 0;  costb = 0;
            for i = 1:num+1
                costa = costa + d(a(i),a(i+1));
                costb = costb + d(b(i),b(i+1));
            end
            [a,costa] = get_best(a,costa,d);
            [b,costb] = get_best(b,costb,d);
            
            if costa<cost(geni*2-1)
                gen(geni*2-1,:) = a;    cost(geni*2-1) = costa;
            end
            if costb<cost(geni*2)
                gen(geni*2,:) = b;    cost(geni*2) = costb;
            end
            
        else
            c = [gen(geni*2,:); gen(geni*2-1,:)];  c = unique(c(:));
            costc = 0;
            for i = 1:num+1
                costc = costc + d(c(i),c(i+1));
            end
            [c,costc] = get_best(c,costc,d);
            
            if costc<cost(geni*2-1)
                gen(geni*2-1,:) = c;    cost(geni*2-1) = costc;
            end
            
            g = zeros(1,num+2); g(1)=1;   g(end)=num+2;
            g(2:end-1)=randperm(num)+1;
            costg = 0;
            for i = 1:num+1
                costg = costg + d(g(i),g(i+1));
            end 
            [gen(geni*2,:),cost(geni*2)] = get_best(g,costg,d);
        end
    end
    
    % 种群排序
    [~,location] = sort(cost);
    gen = gen(location,:);
    cost = cost(location);

    % 显示
    if cost(1)<cost_best
        clf
        p1=plot(sj(gen(1,:),1),sj(gen(1,:),2),'-b','LineWidth',1.5);
        hold on
        p2 = plot(sj(gen_best,1),sj(gen_best,2),'-.r','LineWidth',4);
        p2.Color(4) = 0.3;
        hold on
        legend([p1,p2],["当下最优结果","前一次最优结果"],'location','best');
        xlabel(sprintf("第%d次结果有提高，cost=%.6f",count,cost(1)));
        pause(0.1)
        
        cost_best = cost(1);
        gen_best = gen(1,:);
    end
    
end
fprintf("最好的结果，第i个点在第gen_best(i)次经过时距离最短");
gen_best
%% 获取局部最优
function [x,cost] = get_best(x,cost,d)
num = 100;
for t = 1:8
    
    
    % 互换位置
    for i = 2:num
        for j = i+2:num+1
            a = d(x(i),x(i+1))+d(x(i),x(i-1))+d(x(j),x(j+1))+d(x(j),x(j-1));
            b =  d(x(j),x(i+1))+d(x(j),x(i-1))+d(x(i),x(j+1))+d(x(i),x(j-1));
            temp = a-b;
            if temp>0
                c = x(i);
                x(i) = x(j);    x(j) = c;
                cost = cost-temp;
            end
        end
    end
    
    % 调换位置
    for i = 2:num
        for j = i+1:num+1
            a = d(x(i),x(i-1))+d(x(j),x(j+1));
            b = d(x(j),x(i-1))+d(x(i),x(j+1));
            temp = a-b;
            if temp>0
                x(i:j) = x(j:-1:i);
                cost = cost-temp;
            end
        end
    end
    
    
end

end
