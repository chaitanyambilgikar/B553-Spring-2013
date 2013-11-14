%{
 THIS IS THE MAIN CODE FILE
%}

%nmoiz and cbilgika
%---------------PLEASE CHANGE USERNAME/DIRECTORY IN LINES 6 and 64 IF USING
%FROM IU MACHINE-----------------------------------------------------------
%This is the directory where the training images are stored.
%Change as required. 
Directory = fullfile ('C:','Users','nmoiz','Documents','MATLAB','small');
Images = dir(fullfile(Directory,'*.png'));
Names = {Images.name}';
num_images = numel(Names);
images = [];
for k=1: numel(Names)
       I = imread( fullfile( Directory, Names{k} ) );
       grayI = rgb2gray(I);
       [rows_large,cols_large] = size(grayI);
       grayI = double(grayI);
       grayI = reshape(grayI,rows_large*cols_large,1);
       images = horzcat(images,grayI);
       %images(:,k)=grayI(:);
       
         %operations
end

%size of vector of high res images - N = rows_large*cols_large
N = rows_large.*cols_large;
%Step 1
%{
%Old method for finding eigen vectors - changed on 04/25/2013
Ref: http://blog.cordiner.net/2010/12/02/eigenfaces-face-recognition-matlab/
%}
%compute the mean image - size (Nx1) = (12228x1)
mean_faces = mean(images,2);
%not needed here - done as per ref
shifted_images = images - repmat(mean_faces, 1, num_images);
%compute the eigenvectors - PC [NxN] and eigenvalues - evalues [Nx1]
[PC, score, evalues] = princomp(images');
%making lambda as a diagonal matrix of evalues
lambda = diag(evalues);
%setting l - we will take the l best eigenvectors
l=15;
%taking the top l eigenvectors in B - [Nxl]
B = PC(:,1:l);
%take the l best eigenvalues as well - lambda is now [lxl]
lambda = lambda(1:l,1:l);


%{
New method for finding eigen vectors - writen on 04/25/2013 - image size
Ih changed from 265x192 to 128x96 to save memory. N = 12288

%calcualte size of the entire image matrix - should be 12288 X 16 
[m,n] = size(images);
%find the mean of this matrix
mean_faces = mean(images,2);
%subtract off the mean for each dimension
new_images = images - repmat(mean_faces,1,n);
covariance = new_images*new_images';
covariance = covariance./(n-1);
%}



%getting the low resolution image in Il - [Mx1]
% This is the path of the test image. Change as required
imgl =  imread('C:\Users\nmoiz\Documents\MATLAB\test\test_1.png');
imgl = rgb2gray(imgl);
Il = imgl(:);
[rows_small,cols_small] = size(imgl);
%size of image vecotr of low res images = M = rows_small*cols_small
M = rows_small.*cols_small;
Il = double(Il);
%initialize small lambda  - as done in paper - its a scalar
small_lambda = 0.3;
%the downscaling factor A - [MxN]
A = ones(M,N).*0.5;

%Now compute X* - should be [lx1]


X = (B'*(A'*A)*B + small_lambda.*inv(lambda))\B'*A'*(Il - (A*mean_faces));

%Compute IgH
Ig = B*X + mean_faces;
%convert Ig from a vector to a matrix with cols_large columns - see for loop above
%for cols_large

Ig = uint8(Ig);
IG = reshape(Ig,rows_large,cols_large);  %resulting image from step 1, global modeling
imshow(IG); 

% Step 2 begins here, local modeling
%loop to iterate through training images, to compare each patch in each
%training image 




H = fspecial('laplacian'); %initialize lapacian filter
externalPhi = cell(16,16); %initialize matrix with 16x16 zeros to store phi values

for i=1:16,  %initialize cells to be 8x6 matrices to represent patches
    for j=1:16,
        externalPhi{i,j} = zeros(8,6); %fill these cells with 8x6 matrix of 0.0s where each matrix is a patch
    end
end

for k=1: numel(Names) %iterate through images
       I = imread( fullfile( Directory, Names{k} ) );
       grayI = rgb2gray(I); %convert to grayscale
       [rows,cols] = size(grayI);
       grayI = double(grayI);
       %for each image we want to iterate through its patches which are 8x6, and over lap each other by 2 pixels
       for m=1:8:128,
           for n=1:6:96,
               trainingpatch=[]; %variable to store patch from current training image i
               globalpatch=[]; %variable to store patch from image from step 1
               %handle cases where patches are at the edge and dont need to
               %overlap
               if (m==1 && n==1) || (m==121 && n==1) || (m==1 && n==91) || (m==121 && n==91) %handle edge patch cases
                   temp1 = m+7;
                   temp2 = n+5;
                   trainingpatch = grayI(m:temp1, n:temp2); %patch from training image
                   globalpatch = IG(m:temp1, n:temp2);  %patch from image acquired from step 1
               elseif (m==1 || m==121) && (n~=1 || n~=91)
                   temp1 = m+7;
                   temp2 = n-2; %case where m is edge values, and n is not
                   temp3 = n+3;                   
                   trainingpatch = grayI(m:temp1, temp2:temp3);
                   globalpatch = IG(m:temp1,temp2:temp3);
               elseif (m~=1 || m~=121) && (n==1 || n==91)
                   temp1 = m-2;
                   temp2 = m+5; %case where n is edge values, and m is not
                   temp3 = n+5;                   
                   trainingpatch = grayI(m:temp1, temp2:temp3);
                   globalpatch = IG(temp1:temp2,n:temp3);
               else
                   trainingpatch = grayI(m-2:m+3, n-2:n+5); %-2 and +4 because we 
                   globalpatch = IG(m-2:m+3, n-2:n+5);      %want patches to overlap
               end
               if isempty(trainingpatch) == 0  && isempty(globalpatch) == 0
                   lapGlobalPatch = imfilter(globalpatch,H); %laplacian of global image patch
                   lapTrainPatch = imfilter(trainingpatch,H); %laplacian of training image patch
                   patch = double(lapGlobalPatch) - lapTrainPatch;
                   patch = norm(patch);
                   patch = patch.^2; %squared differences of two laplacians
                                 %next we compute the dirac delta function. 
                   globalpatch = double(globalpatch);
                   patch2 = dirac(globalpatch - trainingpatch);
                   PhiRowIndex = (m+7)/8;  %some arithmetic done to make sure we store 
                   PhiColIndex = (n+5)/6;  %to correct index in the large matrix externalPhi
                   %error mesages while performing next addition, hence
                   %commented out. These next two lines were supposed to
                   %sum up to complete the external component of the Gibbs
                   %function as described in the paper. 
                   %tempmatrix = externalPhi{PhiRowIndex, PhiColIndex} + (patch .* patch2);
                   %externalPhi{PhiRowIndex, PhiColIndex} = (patch .* patch2) ; %Measures how likely patch measures patch on other training image.
               end       
           end
       end
       
end  
               
               
       




