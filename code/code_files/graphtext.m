function [] = graphtext(level,z_vals, tile, num, varargin)
%takes a row vector of column vectors
%each column vector is a standard address with a z-value appended on the
%end
g = GraphApprox(level);
plotting_points = zeros([level+1, (3^(level+1)+3)/2]);
for i = 1:(3^(level+1)+3)/2
    plotting_points(1:end,i) = g.vertices(i).address';
end
%gasketdata = [plotting_points; z_vals];
gasketdata = plotting_points;
%declaring empty output
X = zeros(length(gasketdata),1);
Y = zeros(length(gasketdata),1);
%defining IFS
f2 = @(x) 1/2*(x-[0.5,sqrt(3)]) + [0.5,sqrt(3)];
f1 = @(x) 1/2*(x-[1,0]) + [1,0];
f0 = @(x) 1/2*(x-[0,0]);
%for each point...
for i = 1:length(gasketdata)
    %find the proper qi...
    if gasketdata(1,i) == 0
        cartesian = [0,0];
    end
    if gasketdata(1,i) == 1
        cartesian = [1,0];
    end
    if gasketdata(1,i) == 2
        cartesian = [0.5,sqrt(3)];
    end
    %...and iterate the funciton maps to get the actual point
    for j = length(gasketdata(:,i)):-1:2
        if gasketdata(j,i) == 0
            cartesian = f0(cartesian);
        end
        if gasketdata(j,i) == 1
            cartesian = f1(cartesian);
        end
        if gasketdata(j,i) == 2
            cartesian = f2(cartesian);
        end
    end
    %insert this x,y,z triple into the cartesian vectors for output
    X(i,1) = cartesian(1);
    Y(i,1) = cartesian(2);
    Z(i,1) = z_vals(i);
end
hold on
scatter(X,Y,'.r')
for i = tile
    tileX = X(tile);
    tileY = Y(tile);
    scatter(tileX,tileY, 'r', 'filled');
end

%Z = 1:(3^(level+1)+3)/2;
if num
    nonzeroX = X(Z ~= 0);
    nonzeroY = Y(Z ~= 0);
    nonzeroZ = Z(Z ~= 0);
    text(nonzeroX,nonzeroY+0.05,num2str(double(nonzeroZ)), 'FontSize', 15)
else
    text(X+0.01,Y+0.06,Z,'FontSize', 15)
end
[numedges, ~] = size(g.edges);
for  i = 1:numedges
    Xedge = X(g.edges(i,1:2));
    Yedge = Y(g.edges(i,1:2));
    plot(Xedge, Yedge, 'k');
end

if nargin > 4
    for i = 5:nargin
        data = varargin{i-4};
        if data{1} == "p"
            plot(data{2:end});
        else
            text(data{2:end}, 'FontSize', 15);
        end
    end
end 
    ax = gca;
    ax.XColor = 'none';
    ax.YColor = 'none';
    
    hold off
end