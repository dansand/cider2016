function [ error ] = junk( A )
%SRAY Fill in G matrix according to surface wave ray theory

size(A)
A(1,1)=A(1,1)+1;
assignin('base','G',A);
end

