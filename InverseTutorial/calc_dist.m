function [ delta ] = calc_dist( blat,blon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dim=size(blat,1);
delta=zeros(dim,dim);
%dist=spalloc(dim,dim,round(0.5*(dim-1)*dim)); %sparse is much slower?!
for i=1:dim-1
    delta(i,i+1:dim)=m_idist(blon(i),blat(i),blon(i+1:dim),blat(i+1:dim),'sphere')*180/(pi*6370997.0); %radius used in m_idist
end
delta=delta+delta'; %fill in whole matrix
end

