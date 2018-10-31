function Compute_Paircorrelation(varargin)%varargin��ʾ����������б����в����ĸ�����nargin��ʾ
switch nargin
case 2%�����������Ϊ2ʱ����������ط���
    Im=varargin{1}; 
    maxrad1=varargin{2};%��һ������Ϊͼ�����ݣ��ڶ�������Ϊ��������ط��������뾶
    [gr,gr_err,radius1,gr_2D,Npart,Area_Im]=Compute_Autocorrelation(Im,maxrad1);%����Compute_Autocorrelation����
    assignin('base','gr',gr);
    assignin('base','gr_err',gr_err);     
    assignin('base','radius1',radius1);
    assignin('base','gr_2D',gr_2D);     
    assignin('base','Npart',Npart);
    assignin('base','Area_Im',Area_Im);%assignin��ʵ�ֲ�ͬm�����ı����Ĺ���
case 3%�����������Ϊ3ʱ�����л���ط���
    Im1=varargin{1}; 
    Im2=varargin{2};
    maxrad1=varargin{3};%��һ��������Ϊ��ͼ�����ݡ�����������Ϊ���л���ط��������뾶
    [cr,cr_err,radius1,cr_2D,Npart_Im1,Npart_Im2,Area_Im]=Compute_Crosscorrelation(Im1,Im2,maxrad1);%����Compute_Crosscorrelation����
    assignin('base','cr',cr);
    assignin('base','cr_err',cr_err);
    assignin('base','radius1',radius1);
    assignin('base','cr_2D',cr_2D);
    assignin('base','Npart_Im1',Npart_Im1); 
    assignin('base','Npart_Im2',Npart_Im2);       
    assignin('base','Area_Im',Area_Im);
end
%����ط���������Compute_Autocorrelation�������£�
function [gr,gr_err,radius1,gr_2D,Npart,Area_Im]=Compute_Autocorrelation(Im,maxrad1)
win1=ones(size(Im));%����ͼ��ߴ��ȫ1����
radius1=0:maxrad1;
Imsize1=min(size(Im));
if Imsize1<1.25*maxrad1     
    maxrad1=round(Imsize1/1.25);%����ȡ��
end
Area_Im=sum(sum(win1));%sum������ͣ�Area_Im��ȫ1�����ܵ���
Npart=sum(sum(Im));%Npart��ͼƬ���ܵ���
Im =double(Im);        
dim1 = size(Im,1)+maxrad1;%ͼ������+���뾶
dim2 = size(Im,2)+maxrad1;%ͼ������+���뾶
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
%����غ�����
function [cr,cr_err,radius1,cr_2D,Npart_Im1,Npart_Im2,Area_Im] =Compute_Crosscorrelation(Im1,Im2,maxrad1)
win1=ones(size(Im1));%ͼ1��ȫ1����
radius1=0:maxrad1;
Imsize1=min(size(Im1));
if Imsize1<1.25*maxrad1
    maxrad1=round(Imsize1/1.25);
end
Area_Im=sum(sum(win1));%�ܵ���
Npart_Im1=sum(sum(Im1));%ͼ1�ܵ���
Npart_Im2=sum(sum(Im2));%ͼ2�ܵ���
Im1 =double(Im1);Im2=double(Im2);%ת��double����Ϊ������Ҷ�任��Ҫdouble������
dim11=size(Im1,1)+maxrad1; dim12 = size(Im1,2)+maxrad1;%ͼ1����+maxrad1
dim21=size(Im2,1)+maxrad1; dim22 = size(Im2,2)+maxrad1;
NormF =real(fftshift(ifft2(abs(fft2(win1,dim11,dim12)).^2)));%ʵ��ȫ1����ľ��
%��0���win1������ά����Ҷ�任��abs�������ƽ������Ƶ��ĳ˻����渵��Ҷ�任ת���ɾ��
%���maxrad1������Ҷ�任�ŵ��ھ��
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
    indx_c=bin_indx==ii;%bin_indx��ÿ��Ԫ����ii�Ƚϣ���ͬ�򷵻�1
    if npixc>0
        cr(ii)=sum(indx_c.*CrossCorr)/npixc;
        cr_err(ii)=sqrt(sum(indx_c.*(CrossCorr-cr(ii)).^2))/npixc;
    end
end
end
end
