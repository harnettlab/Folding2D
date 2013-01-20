function kerchief
%roll up a 2D shape from a curvature map, starting in the middle
inputx=40; %Array size, if it's even--add 1
NX=2*floor(inputx/2)+1;%now there will be a point in the middle of the array
x=linspace(-0.5,0.5,NX);
y=linspace(-0.5,0.5,NX);
[X Y]=meshgrid(x,y);
alpha=0.2;%how much to scale down z-displacement compared to tile width
diffx=x(2)-x(1); %tile width
Z=zeros(size(X));
%X Y and Z coords of vertices

XYZ(:,:,1)=X;
XYZ(:,:,2)=Y;
XYZ(:,:,3)=Z; %try to keep coords together like this for compactness
XYZBase=XYZ;%Store initial guesses for X,Y, and Z

%Generate a grid of tiles, each vertex is the corner of 4 tiles
myseq=seq(NX,'randbowl');
dummyseq=padarray(myseq,[1 1],'replicate','pre');
dzloc=(filter2([1 1;1 1]/4,dummyseq,'valid')-0.5)*diffx*alpha
%the above gives a NX by NX array of local delta z for each vertex
%local delta z means, how far that vertex is perpendicular from the plane
%that is the best fit to its four nearest neighbors on the grid

%get it started by generating center, edge, and then corner points around the 3x3 center of the XYZ array
rowmid=ceil(NX/2);
colmid=ceil(NX/2);%currently it's a square

XYZ(rowmid,colmid,:)=[0 0 0]; %set midpoint to be fixed here
XYZ(rowmid,colmid+1,3)=+dzloc(rowmid,colmid+1);%displace nearest neighbors along z by dzloc at their address.
XYZ(rowmid,colmid-1,3)=+dzloc(rowmid,colmid-1);
XYZ(rowmid+1,colmid,3)=+dzloc(rowmid+1,colmid);
XYZ(rowmid-1,colmid,3)=+dzloc(rowmid-1,colmid);
%flipped the sign of dzloc in the initial 3x3 square edges to get it started the
%right direction, since there's no farther out points to guess with, but it
%doesn't seem to matter in the end

%now fill in corners so it's a 3x3 square in the middle of the array
%upper left corner of 3x3 square
A=squeeze(XYZ(rowmid-1,colmid,:)-XYZ(rowmid,colmid,:));
B=squeeze(XYZ(rowmid,colmid-1,:)-XYZ(rowmid,colmid,:));
mynormal=cross(A,B)/norm(A)/norm(B);
XYZBase(rowmid-1,colmid-1,:)=XYZ(rowmid,colmid-1,:)+XYZ(rowmid-1,colmid,:)-XYZ(rowmid,colmid,:);
XYZ(rowmid-1,colmid-1,:)=XYZBase(rowmid-1,colmid-1,:)+shiftdim((dzloc(rowmid-1,colmid-1)*mynormal)',-1);

%lower left corner of 3x3 square
A=squeeze(XYZ(rowmid,colmid-1,:)-XYZ(rowmid,colmid,:));
B=squeeze(XYZ(rowmid+1,colmid,:)-XYZ(rowmid,colmid,:));
mynormal=cross(A,B)/norm(A)/norm(B);
XYZBase(rowmid+1,colmid-1,:)=XYZ(rowmid,colmid-1,:)+XYZ(rowmid+1,colmid,:)-XYZ(rowmid,colmid,:);
XYZ(rowmid+1,colmid-1,:)=XYZBase(rowmid+1,colmid-1,:)+shiftdim((dzloc(rowmid+1,colmid-1)*mynormal)',-1);

%upper right corner of 3x3 square
A=squeeze(XYZ(rowmid,colmid+1,:)-XYZ(rowmid,colmid,:));
B=squeeze(XYZ(rowmid-1,colmid,:)-XYZ(rowmid,colmid,:));
mynormal=cross(A,B)/norm(A)/norm(B);
XYZBase(rowmid-1,colmid+1,:)=XYZ(rowmid,colmid+1,:)+XYZ(rowmid-1,colmid,:)-XYZ(rowmid,colmid,:);
XYZ(rowmid-1,colmid+1,:)=XYZBase(rowmid-1,colmid+1,:)+shiftdim((dzloc(rowmid-1,colmid+1)*mynormal)',-1);

%lower right corner of 3x3 square
A=squeeze(XYZ(rowmid+1,colmid,:)-XYZ(rowmid,colmid,:));
B=squeeze(XYZ(rowmid,colmid+1,:)-XYZ(rowmid,colmid,:));
mynormal=cross(A,B)/norm(A)/norm(B);
XYZBase(rowmid+1,colmid+1,:)=XYZ(rowmid,colmid+1,:)+XYZ(rowmid+1,colmid,:)-XYZ(rowmid,colmid,:);
XYZ(rowmid+1,colmid+1,:)=XYZBase(rowmid+1,colmid+1,:)+shiftdim((dzloc(rowmid+1,colmid+1)*mynormal)',-1);

%stencil=[0 1 0;1 0 1;0 1 0]/4;%NSEW nearest neighbor stencil
stencil=[0.707 1 0.707 ;1 0 1; 0.707 1 0.707];%this stencil seemed to keep squares from going diagonal better
stencil=stencil/sum(sum(stencil));

 for myextent=2:(floor(NX/2))%this is the main loop. grow the array outward by one each time, iterate, repeat
    rowmin=rowmid-myextent;
    rowmax=rowmid+myextent;
    colmin=colmid-myextent;
    colmax=colmid+myextent;
    figure()
    surf(XYZ((rowmin+1):(rowmax-1),(colmin+1):(colmax-1),1),XYZ((rowmin+1):(rowmax-1),(colmin+1):(colmax-1),2),XYZ((rowmin+1):(rowmax-1),(colmin+1):(colmax-1),3)-XYZ(rowmid,colmid,3),myseq((rowmin+1):(rowmax-1),(colmin+1):(colmax-1)))
    curfigname=strcat('randbowl',int2str(myextent),'.tif');
    axis([-0.6 0.6 -0.6 0.6 0 0.5])
    axis off
    colormap spring
    curfighandl=gcf;%grab the latest plot
    saveas(curfighandl,curfigname,'tif');
   
   %create new edge points--except corners--by extrapolation, this means change XYZ array
    for r=(rowmin+1):(rowmax-1)
          A=XYZ(r,colmin+1,:)-XYZ(r,colmin+2,:);%difference vector between previous two points
          XYZ(r,colmin,:)=XYZ(r,colmin+1,:)+A/norm(squeeze(A))*diffx;%left edge points--try scaling by diffx so they don't stretch
          A=XYZ(r,colmax-1,:)-XYZ(r,colmax-2,:);          %though they don't stretch, angles can still go non-90 degrees
          XYZ(r,colmax,:)=XYZ(r,colmax-1,:)+A/norm(squeeze(A))*diffx;%right edge points
    end
    for c=(colmin+1):(colmax-1)
          A=XYZ(rowmin+1,c,:)-XYZ(rowmin+2,c,:);
          XYZ(rowmin,c,:)=XYZ(rowmin+1,c,:)+A/norm(squeeze(A))*diffx;%top edge points
          A=XYZ(rowmax-1,c,:)-XYZ(rowmax-2,c,:);
          XYZ(rowmax,c,:)=XYZ(rowmax-1,c,:)+A/norm(squeeze(A))*diffx;%bottom edge points
    end    

   %do nearest neighbor averaging and normal-finding on the inner array,
   %iterating enough times for convergence 
   %was going up to 4*myextent^2
   for i=1:(2*myextent)%scale number of iterations by approximate number of new array elements
       XYZBase((rowmin+1):(rowmax-1),(colmin+1):(colmax-1),1)=filter2(stencil,XYZ(rowmin:rowmax,colmin:colmax,1),'valid');%it is going to be 2 smaller than the rowmin-rowmax
       XYZBase((rowmin+1):(rowmax-1),(colmin+1):(colmax-1),2)=filter2(stencil,XYZ(rowmin:rowmax,colmin:colmax,2),'valid');%it is going to be 2 smaller than the rowmin-rowmax
       XYZBase((rowmin+1):(rowmax-1),(colmin+1):(colmax-1),3)=filter2(stencil,XYZ(rowmin:rowmax,colmin:colmax,3),'valid');%it is going to be 2 smaller than the rowmin-rowmax
       %then find normals and add dzloc, maybe there is a faster way to do
       %it but I will go point by point
       for myrow=(rowmin+1):(rowmax-1)
           for mycol=(colmin+1):(colmax-1)
               A=squeeze(XYZBase(myrow-1,mycol,:)-XYZBase(myrow+1,mycol,:));
               B=squeeze(XYZBase(myrow,mycol-1,:)-XYZBase(myrow,mycol+1,:));
               mynormal=cross(A,B)/norm(A)/norm(B);
               XYZ(myrow,mycol,:)=XYZBase(myrow,mycol,:)+shiftdim((dzloc(myrow,mycol)*mynormal)',-1);
           end
       end  
       %and recalc edge points--except corners--by extrapolation
       for r=(rowmin+1):(rowmax-1)
          A=XYZ(r,colmin+1,:)-XYZ(r,colmin+2,:);
          XYZ(r,colmin,:)=XYZ(r,colmin+1,:)+A/norm(squeeze(A))*diffx;%left edge points--try scaling by diffx so they don't stretch
          A=XYZ(r,colmax-1,:)-XYZ(r,colmax-2,:);
          XYZ(r,colmax,:)=XYZ(r,colmax-1,:)+A/norm(squeeze(A))*diffx;%right edge points
       end
       for c=(colmin+1):(colmax-1)
          A=XYZ(rowmin+1,c,:)-XYZ(rowmin+2,c,:);
          XYZ(rowmin,c,:)=XYZ(rowmin+1,c,:)+A/norm(squeeze(A))*diffx;%top edge points
          A=XYZ(rowmax-1,c,:)-XYZ(rowmax-2,c,:);
          XYZ(rowmax,c,:)=XYZ(rowmax-1,c,:)+A/norm(squeeze(A))*diffx;%bottom edge points
       end 
   
   end
   %finally guess new corner points by using three neighboring points--put in
   %xyz array

       %upper left corner
       A=squeeze(XYZ(rowmin,colmin+1,:)-XYZ(rowmin+1,colmin+1,:));
       B=squeeze(XYZ(rowmin+1,colmin,:)-XYZ(rowmin+1,colmin+1,:));
       mynormal=cross(A,B)/norm(A)/norm(B);
       XYZBase(rowmin,colmin,:)=XYZ(rowmin,colmin+1,:)+XYZ(rowmin+1,colmin,:)-XYZ(rowmin+1,colmin+1,:);
       XYZ(rowmin,colmin,:)=XYZBase(rowmin,colmin,:);%+shiftdim((dzloc(rowmin,colmin)*mynormal)',-1);
       
       %lower left corner
       A=squeeze(XYZ(rowmax-1,colmin,:)-XYZ(rowmax-1,colmin+1,:));
       B=squeeze(XYZ(rowmax,colmin+1,:)-XYZ(rowmax-1,colmin+1,:));
       mynormal=cross(A,B)/norm(A)/norm(B);
       XYZBase(rowmax,colmin,:)=XYZ(rowmax-1,colmin,:)+XYZ(rowmax,colmin+1,:)-XYZ(rowmax-1,colmin+1,:);
       XYZ(rowmax,colmin,:)=XYZBase(rowmax,colmin,:);%+shiftdim((dzloc(rowmax,colmin)*mynormal)',-1);
       
       %upper right corner 
       B=squeeze(XYZ(rowmin,colmax-1,:)-XYZ(rowmin+1,colmax-1,:));
       A=squeeze(XYZ(rowmin+1,colmax,:)-XYZ(rowmin+1,colmax-1,:));
       mynormal=cross(A,B)/norm(A)/norm(B);
       XYZBase(rowmin,colmax,:)=XYZ(rowmin,colmax-1,:)+XYZ(rowmin+1,colmax,:)-XYZ(rowmin+1,colmax-1,:);
       XYZ(rowmin,colmax,:)=XYZBase(rowmin,colmax,:);%+shiftdim((dzloc(rowmin,colmax)*mynormal)',-1);
       
       %lower right corner 
       B=squeeze(XYZ(rowmax-1,colmax,:)-XYZ(rowmax-1,colmax-1,:));
       A=squeeze(XYZ(rowmax,colmax-1,:)-XYZ(rowmax-1,colmax-1,:));
       mynormal=cross(A,B)/norm(A)/norm(B);
       XYZBase(rowmax,colmax,:)=XYZ(rowmax,colmax-1,:)+XYZ(rowmax-1,colmax,:)-XYZ(rowmax-1,colmax-1,:);
       XYZ(rowmax,colmax,:)=XYZBase(rowmax,colmax,:);%+shiftdim((dzloc(rowmax,colmax)*mynormal)',-1);
   %myextent will get incremented the next time through
 end

%XYZ(:,:,3)=XYZ(:,:,3)-XYZ(rowmid,colmid,3);% shift everything else by the local dz at the origin
%hold off
%surf(XYZ(:,:,1),XYZ(:,:,2),XYZ(:,:,3),myseq)
%mymap=[0.75 0.75 0.25; 0.5 0.25 0.25];
%colormap jet
%shading faceted
%hold on
%contour(XYZ(:,:,1),XYZ(:,:,2),XYZ(:,:,3))
beep
myseq