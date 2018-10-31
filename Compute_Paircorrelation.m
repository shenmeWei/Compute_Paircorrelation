function Compute_Paircorrelation(varargin)%varargin表示变输入参数列表，其中参数的个数用nargin表示
switch nargin
case 2%输入参数个数为2时，进行自相关分析
    Im=varargin{1}; 
    maxrad1=varargin{2};%第一个参数为图像内容，第二个参数为进行自相关分析的最大半径
    [gr,gr_err,radius1,gr_2D,Npart,Area_Im]=Compute_Autocorrelation(Im,maxrad1);%调用Compute_Autocorrelation函数
    assignin('base','gr',gr);
    assignin('base','gr_err',gr_err);     
    assignin('base','radius1',radius1);
    assignin('base','gr_2D',gr_2D);     
    assignin('base','Npart',Npart);
    assignin('base','Area_Im',Area_Im);%assignin：实现不同m函数的变量的共享
case 3%输入参数个数为3时，进行互相关分析
    Im1=varargin{1}; 
    Im2=varargin{2};
    maxrad1=varargin{3};%第一二个参数为两图像内容。第三个参数为进行互相关分析的最大半径
    [cr,cr_err,radius1,cr_2D,Npart_Im1,Npart_Im2,Area_Im]=Compute_Crosscorrelation(Im1,Im2,maxrad1);%调用Compute_Crosscorrelation函数
    assignin('base','cr',cr);
    assignin('base','cr_err',cr_err);
    assignin('base','radius1',radius1);
    assignin('base','cr_2D',cr_2D);
    assignin('base','Npart_Im1',Npart_Im1); 
    assignin('base','Npart_Im2',Npart_Im2);       
    assignin('base','Area_Im',Area_Im);
end
%自相关分析函数：Compute_Autocorrelation函数如下：
function [gr,gr_err,radius1,gr_2D,Npart,Area_Im]=Compute_Autocorrelation(Im,maxrad1)
win1=ones(size(Im));%生成图像尺寸的全1矩阵
radius1=0:maxrad1;
Imsize1=min(size(Im));
if Imsize1<1.25*maxrad1     
    maxrad1=round(Imsize1/1.25);%向上取整
end
Area_Im=sum(sum(win1));%sum：列求和；Area_Im是全1矩阵总点数
Npart=sum(sum(Im));%Npart是图片的总点数
Im =double(Im);        
dim1 = size(Im,1)+maxrad1;%图像行数+最大半径
dim2 = size(Im,2)+maxrad1;%图像列数+最大半径
NormF = real(fftshift(ifft2(abs(fft2(win1,dim1,dim2)).^2)));
gr_2D_full= Area_Im^2/Npart^2*real(fftshift(ifft2(abs(fft2(Im,dim1,dim2)).^2)))./NormF;
gr_2D= imcrop(gr_2D_full,[floor(dim1/2+1)-maxrad1,floor(dim2/2+1)-maxrad1,2*maxrad1,2*maxrad1]);
x1=ones(1, 2*maxrad1+1)'*(-maxrad1:maxrad1);
y1=(-maxrad1:maxrad1)'*ones(1, 2*maxrad1+1);
[~,rad1,hgh1] = cart2pol(x1,y1,gr_2D);
CorrRad1=reshape(rad1,1,(2*maxrad1+1)^2);
AutoCorr1=reshape(hgh1,1,(2*maxrad1+1)^2);
[CorrRad,index1] = sort(CorrRad1); AutoCorr=AutoCorr1(index1);
rad2=0:floor(max(CorrRad));
[np1,bin_indx]=histc(CorrRad,rad2-.5);
gr=zeros(1,maxrad1+1);
gr_err=zeros(1,maxrad1+1);
for ii=1:maxrad1+1
    npixc=np1(ii);
    indx_c=bin_indx==ii;
    if npixc>0
        gr(ii)=sum(indx_c.*AutoCorr)/npixc;
        gr_err(ii)=sqrt(sum(indx_c.*(AutoCorr-gr(ii)).^2))/npixc;
    end
end
end
%互相关函数：
function [cr,cr_err,radius1,cr_2D,Npart_Im1,Npart_Im2,Area_Im] =Compute_Crosscorrelation(Im1,Im2,maxrad1)
win1=ones(size(Im1));%图1的全1矩阵
radius1=0:maxrad1;
Imsize1=min(size(Im1));
if Imsize1<1.25*maxrad1
    maxrad1=round(Imsize1/1.25);
end
Area_Im=sum(sum(win1));%总点数
Npart_Im1=sum(sum(Im1));%图1总点数
Npart_Im2=sum(sum(Im2));%图2总点数
Im1 =double(Im1);Im2=double(Im2);%转成double，因为做傅里叶变换需要double型数据
dim11=size(Im1,1)+maxrad1; dim12 = size(Im1,2)+maxrad1;%图1行列+maxrad1
dim21=size(Im2,1)+maxrad1; dim22 = size(Im2,2)+maxrad1;
NormF =real(fftshift(ifft2(abs(fft2(win1,dim11,dim12)).^2)));%实现全1矩阵的卷积
%用0填充win1再做二维傅里叶变换，abs求振幅，平方是做频域的乘积，逆傅里叶变换转换成卷积
%填充maxrad1做傅里叶变换才等于卷积
cr_2D_full=Area_Im^2/(Npart_Im1*Npart_Im2)*real(fftshift(ifft2((fft2(Im1,dim11,dim12).*conj(fft2(Im2,dim21,dim22))))))./NormF;
cr_2D= imcrop(cr_2D_full,[floor(dim11/2+1)-maxrad1,floor(dim12/2+1)-maxrad1,2*maxrad1,2*maxrad1]);
x1=ones(1, 2*maxrad1+1)'*(-maxrad1:maxrad1);
y1=(-maxrad1:maxrad1)'*ones(1, 2*maxrad1+1);
[~,rad1,hgh1] = cart2pol(x1,y1,cr_2D);
CorrRad1=reshape(rad1,1,(2*maxrad1+1)^2);
CrossCorr1=reshape(hgh1,1,(2*maxrad1+1)^2);
[CorrRad,index1] = sort(CorrRad1); CrossCorr=CrossCorr1(index1);
rad2=0:floor(max(CorrRad));
[np1,bin_indx]=histc(CorrRad,rad2-.5);
cr=zeros(1,maxrad1+1);
cr_err=zeros(1,maxrad1+1);
for ii=1:maxrad1+1
    npixc=np1(ii);
    indx_c=bin_indx==ii;%bin_indx中每个元素与ii比较，相同则返回1
    if npixc>0
        cr(ii)=sum(indx_c.*CrossCorr)/npixc;
        cr_err(ii)=sqrt(sum(indx_c.*(CrossCorr-cr(ii)).^2))/npixc;
    end
end
end
end
