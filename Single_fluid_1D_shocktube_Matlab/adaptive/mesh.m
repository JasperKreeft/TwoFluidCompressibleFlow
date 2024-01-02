clear all
clf%ose all
clc

L = 1.0;
N = 10;

x_nodes = linspace(0,L,N+1);
dx(1) = L/N;
dx(2) = dx(1)/2;
dx(3) = dx(2)/2;
x_cellcenter = (x_nodes(1:end-1)+x_nodes(2:end))/2;
plot(x_nodes,zeros(N+1,1),'+-','linewidth',2)
hold on
plot(x_cellcenter,zeros(N,1),'ro','linewidth',2)
grid on

%               cellnr  positie       level     neighbours
celldata = [ (3:N+2)' x_cellcenter' ones(N,1) [1 3:N+1]' [4:N+2 2]'];

nr = N+2;
cellstosplit = sort([9]);

i=0;
for celltosplit = cellstosplit(end:-1:1)

    i=i+1;

    nrofcells = length(celldata);
    rowtosplit = (1:nrofcells)*(celldata(:,1)==celltosplit);
    celldata = [celldata;
                nr+1 x_cellcenter(rowtosplit)-dx(celldata(rowtosplit,3)+1)/2 celldata(rowtosplit,3)+1 celldata(rowtosplit,4) nr+2;
                nr+2 x_cellcenter(rowtosplit)+dx(celldata(rowtosplit,3)+1)/2 celldata(rowtosplit,3)+1 nr+1 celldata(rowtosplit,5)];
    celldata(rowtosplit-1,5) = nr+1;
    celldata(rowtosplit+1,4) = nr+2;

    celldata(rowtosplit,:) = [];
    nr = nr+2;

    plot(celldata(:,2),i*ones(length(celldata),1),'go','linewidth',2)
    ylim([-1 length(cellstosplit)+1]);
end