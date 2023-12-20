function [D, alpha]=euclid_dist(a,varargin)
%compute the Euclidian distance (D) and the angle relative to the horizontal (alpha) between points in a 2 column vector (if one input)
% or sets of points a and b if two inputs
%rows are individual observations, columns are x and then y value of the
%point

if isempty(varargin)
    %compute differences between points in 2col vector
    if size(a,1)<=2
        disp('error, vector contains only 1 point')
    end
    
    D=zeros(size(a,1)-1,1);
    alpha=zeros(size(a,1)-1,1);
    
    for i=1:size(a,1)-1
        D(i)= sqrt((a(i,1)-a(i+1,1))^2 + (a(i,2)-a(i+1,2))^2);
%         alpha(i)=-atand((a(i,2)-a(i+1,2))/(a(i,1)-a(i+1,1)));
        alpha(i)=atan2d(a(i,1)-a(i+1,1),a(i,2)-a(i+1,2));
    end
    
else
    b=varargin{1};
    %if one of the inputs is only a point, the other one is multiple,
    %repmat the single point to the same size
    if size(a,1)~=size(b,1) || size(a,1)==1  || size(b,1)==1
        if size(a,1)==1 
            a=repmat(a,size(b,1),1);
        elseif size(b,1)==1 
            b=repmat(b,size(a,1),1);
        end
    end
    
    D=zeros(size(a,1),1);
    alpha=zeros(size(a,1),1);
    
    for i=1:size(a,1)
        D(i)= sqrt((a(i,1)-b(i,1))^2 + (a(i,2)-b(i,2))^2);
%         alpha(i)=-atand((a(i,2)-b(i,2))/(a(i,1)-b(i,1)));
        alpha(i)=atan2d(a(i,1)-b(i,1),a(i,2)-b(i,2));
    end
end

end