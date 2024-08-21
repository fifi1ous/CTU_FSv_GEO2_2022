clc; clear; format long G
fid=fopen("body.txt",'r');
body=fscanf(fid,'%*d %f %f',[2,inf])';
r=size(body,1);
%% podle Y
py=body(1,2)*(body(2,1)-body(end,1));
for n=2:(r-1)
    py1=body(n,2)*(body(n+1,1)-body(n-1,1));
    py=py+py1;
end
py2=body(end,2)*(body(1,1)-body(end-1,1));
py=(1/2)*(py+py2)
%% podle X
px=body(1,1)*(body(end,2)-body(2,2));
for n=2:(r-1)
    px1=body(n,1)*(body(n-1,2)-body(n+1,2));
    px=px+px1;
end
px2=body(end,1)*(body(end-1,2)-body(1,2));
px=(1/2)*(px+px2)