function [] = gasketgraph( plotting_points,z_vals )
%takes a row vector of column vectors
%each column vector is a standard address with a z-value appended on the
%end
dims = size(plotting_points);
gasketdata = [plotting_points; z_vals'];
%declaring empty output
X = zeros(length(gasketdata),1);
Y = zeros(length(gasketdata),1);
Z = zeros(length(gasketdata),1);
%defining IFS
f0 = @(x) 1/2*(x-[0.5,tan(60)*2]) + [0.5,tan(60)*2];
f1 = @(x) 1/2*(x-[1,0]) + [1,0];
f2 = @(x) 1/2*(x-[0,0]) ;
%for each point...
for i = 1:length(gasketdata)
    %find the proper qi...
    if gasketdata(1,i) == 0
        cartesian = [0.5,tan(60)*2];
    end
    if gasketdata(1,i) == 1
        cartesian = [1,0];
    end
    if gasketdata(1,i) == 2
        cartesian = [0,0];
    end
    %...and iterate the funciton maps to get the actual point
    for j = 2:length(gasketdata(:,i))-1
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
    Z(i,1) = gasketdata(end,i);
end
plot3(X,Y,Z,'.r')
end
