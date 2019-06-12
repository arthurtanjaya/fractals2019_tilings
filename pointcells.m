function [ cell1 cell2 ] = pointcells( address )
%given an address, we want the two adjacent cells
%with the cell address acting on K to give the cell

%we can actually get this just using the same address and cutting off the
%first bit, and then using the secondary address for the other side
cell1 = primary(address);
cell2 = secondary(address);
cell1 = cell1(2:end);
cell2 = cell2(2:end);
end
