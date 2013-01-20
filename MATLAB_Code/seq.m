function x=seq(n,whattype)
%Generate a nxn matrix of 0s and 1s in some different types of patterns
if nargin<1
    n=5;
    whattype='check'; %default is a small 5x5 checkerboard
end
if nargin<2
    whattype='check';  %if n specified, default is a nxn checkerboard
end
switch whattype
    case 'check'
        x=checkseq(n);  %generate a checkerboard of 0s and 1s    
    case 'egg'  
       
           eggcrate= [0 1 0 1 0 1 0 1 0 1 0 1;
            1 1 1 1 1 0 1 0 0 0 0 0;
            0 1 1 1 1 1 0 0 0 0 0 1;
            1 1 1 1 1 0 1 0 0 0 0 0;
            0 1 1 1 1 1 0 0 0 0 0 1;
            1 0 1 0 1 0 1 0 1 0 1 0;
            0 1 0 1 0 1 0 1 0 1 0 1;
            1 0 0 0 0 0 1 1 1 1 1 0;
            0 0 0 0 0 1 0 1 1 1 1 1;
            1 0 0 0 0 0 1 1 1 1 1 0;
            0 0 0 0 0 1 0 1 1 1 1 1;
            1 0 1 0 1 0 1 0 1 0 1 0;];%egg crate sequence
        
        x=repmat(eggcrate,ceil(n/12),ceil(n/12));
        x=x(1:n,1:n);%chop off remaining corners
    case 'egg2'%shift off so there is a block of zeros in the center
       
   newegg=[0 0 0 1 1 1 1 1 0 1 0 0;
           0 0 1 0 1 1 1 1 1 0 0 0;
           0 1 0 1 0 1 0 1 0 1 0 1;
           1 0 1 0 1 0 1 0 1 0 1 0;
           1 1 0 1 0 0 0 0 0 1 1 1;
           1 1 1 0 0 0 0 0 1 0 1 1;
           1 1 0 1 0 0 0 0 0 1 1 1; 
           1 1 1 0 0 0 0 0 1 0 1 1; 
           0 1 0 1 0 1 0 1 0 1 0 1; 
           1 0 1 0 1 0 1 0 1 0 1 0; 
           0 0 0 1 1 1 1 1 0 1 0 0;
           0 0 1 0 1 1 1 1 1 0 0 0];
        x=repmat(newegg,ceil(n/12),ceil(n/12));
        x=x(1:n,1:n);%chop off remaining corners
        
    case 'egg3' %make symmetric along diags
   newegg=[0 0 1 0 1 1 1 1 1 0 1 0 0;
           0 0 0 1 1 1 1 1 1 1 0 0 0;
           1 0 1 0 1 0 1 0 1 0 1 0 1;
           0 1 0 1 0 1 0 1 0 1 0 1 0;
           1 1 1 0 0 0 0 0 0 0 1 1 1;
           1 1 0 1 0 0 0 0 0 1 0 1 1;
           1 1 1 0 0 0 0 0 0 0 1 1 1; 
           1 1 0 1 0 0 0 0 0 1 0 1 1;
           1 1 1 0 0 0 0 0 0 0 1 1 1; 
           0 1 0 1 0 1 0 1 0 1 0 1 0;
           1 0 1 0 1 0 1 0 1 0 1 0 1;
           0 0 0 1 1 1 1 1 1 1 0 0 0;
           0 0 1 0 1 1 1 1 1 0 1 0 0];
        x=repmat(newegg,ceil(n/13),ceil(n/13));
        x=x(1:n,1:n);%chop off remaining corners
        
    case 'peaks'%make a predetermined 20x20 array with some flipped tiles
        x=checkseq(n) %Start with a checkerboard
        %Flip a few tiles to see local peaks
        x(5,15)=1-x(5,15);
        x(15,5)=1-x(15,5);
        x(10,11)=1-x(10,11);
        x(15,15)=1-x(15,15);
        x(5,5)=1-x(5,5);
    case 'diag'
        %Create diagonal folds
        %starting from a checkerboard
        x=checkseq(n);
        for i=2:(n-1)
            x(i-1,i)=1;
            x(i+1,i)=1;
            x(i,i)=1;
        end
        %for i=(n-1):-1:2
        %     x(i, n-i)=1;
        %     x(i+1,n-i)=1;
        %end
    case 'bigline'
        x=checkseq(n);
        x(ceil(n/2)+1,:)=1;
        x(ceil(n/2),:)=1;
        x(ceil(n/2)-1,:)=1;
        x(ceil(n/2)-2,:)=1;
        
    case 'line'
        %Create a horizontal fold
        x=checkseq(n);
        x(ceil(n/2)+1,:)=1;
        %x(ceil(n/2)+2,:)=0;
        %x(ceil(n/2)+3,:)=0;
        % Create vertical fold
        %x(:,ceil(n/2))=1;
        %x(:,ceil(n/2)+1)=1;
        %x(:,ceil(n/2)+2)=1;
        %x(:,ceil(n/2)+3)=1;
    case 'dogear'
        x=checkseq(n);
        x(1:ceil(n/4),1:ceil(n/4))=1;
    case 'hill'
        x=zeros(n);
    case 'bowl'
        x=ones(n);
    case 'waterbomb'
        x=checkseq(n);
        x(ceil(n/2)-2,:)=0;
        x(ceil(n/2)+1,:)=0;
        x(ceil(n/2),:)=0;
        x(ceil(n/2)-1,:)=0;
        x(ceil(n/2)-2,:)=0;
        for i=1:n
            x(i,i)=1;
            x(i,n+1-i)=1;
        end
        for i=2:n
            x(i-1,i)=1;
            x(i-1,n+1-i)=1;
        end
        for i=1:(n-1)
            x(i+1,i)=1;
            x(i+1,n+1-i)=1;
        end
        for i=1:(n-2)
            x(i+2,i)=1;
            x(i+2,n+1-i)=1;
        end
    case 'randbowl'
        x=round(rand(n,n)+0.3);
    case 'randstaples'
        %staple=[1 1 1 1 1 1 1];
        x=checkseq(n)
        staploc=logical(round(rand(n,n)-0.45));
        for i=1:7
         %staploc=reshape(staploc,size(x));
             x(circshift(staploc,i))=1;
        end;     
        x=x(1:n,1:n);
    case 'striped'
        x=checkseq(n);
        x(1:4:n,:)=0;
        x(2:4:n,:)=0;
        
    otherwise
        x=checkseq(n);%default is a checkerboard               
end