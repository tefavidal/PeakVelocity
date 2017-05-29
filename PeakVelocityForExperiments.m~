
%% read file
close all;
clear all;

% dir='20170508-New-S2h36m-6mmpmin';   
% name ='Filtered_top-270-800';
% 
% path=['/data.lfpn/eckstein/Torsten/' dir '/'];
% filename = [path name '.tif']

filename ='Filtered_bot-200-700.tif';
path='ble';
name='bla';

info = imfinfo(filename);
A=imread(filename, 'Info', info);
figure()
imagesc(A)
colormap(gray)

%% Stuff


if ~exist('A','var')
    fprintf('Please call your data A \n');
    return
end

S=double(flip(A,1));
Nx=length(S(1,:));
T=length(S(:,1));
%20 seconds dt
dt=1/3;
t=0:dt:dt*(T-1);
t=t';
prompt = 'How many pixels per 2mm? ';
aa = input(prompt);
dx=2/aa;
x=0:dx:dx*(Nx-1);
Lx=max(x);

clearvars tout
clearvars xout
clearvars xout2


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%      Peaks over time detection---> The are saved in out.                     %%
        %%      Takes discrete values (multiples of dx).                                %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_out=[path name '_peaks'];


if exist([data_out '.mat'], 'file') == 2

choice = questdlg(['Data for ' name 'already seems to exist'], ...
	'Peakfinder choice', ...
    'Calculate anyway','Load data','Load data');
% Handle response
switch choice
case 'Calculate anyway'

[xout2,tout]=DetectPeaksAlongX(x,t,S,0.0);
save(data_out,'-v7.3');

case 'Load data'
      calculate=0; 
      fprintf('DO NOT Calculate = %d \n',calculate)
     load(data_out);       
end
else

    [xout2,tout]=DetectPeaksAlongX(x,t,S,0.0);
    save(data_out,'-v7.3');
end

'part 1 done'

%% Part 2
close all;
figure(1)
imagesc(x,t,S)
 set(gca,'YDir','normal')
colormap(gray), hold on 
plot(xout2,tout,'o')
hold off

%% read line

clearvars tout1
clearvars xout1

input('Zoom the picture in the area you want, then press enter \n');

display('Click on the map the position from which you want to start measuring');

[auxX, auxY]=ginput(1);
distance=(auxX*ones(size(xout2))-xout2).^2+(auxY*ones(size(tout))-tout).^2;
[ignore, index]=min(distance);
tstart=tout(index);
xstart=xout2(index);


    forward=1;
    k=1;
    tout1(k)=tstart;
    xout1(k)=xstart;


Wiggle=0.003;

while tstart<t(T-1) && xout1(k) < Lx && xout1 (k) > dx
    k=k+1;
      
        
	  Jump=1;
      while 1          
          [idx, testt]=find_peak(tout,xout2, xout1,Jump,Wiggle, xstart, dx, k, tout1, dt);
          if Jump>10
              break
          elseif idx==-1
              Jump=Jump+1;
          elseif idx==-1 && Jump>=3
              break
          elseif idx ~= -1 && abs(tout1(k-1)-testt(idx,2))>3
              Jump=Jump+1;                  
          else
              break
          
          end
          
      end
      


        if idx==-1 || abs(tout1(k-1)-testt(idx,2))>3 || Jump>10
            xout1=xout1(1:(k-1));
            tout1=tout1(1:(k-1));
            break
        end    
                    
    xstart=xstart+dx*idx;
    xout1(k)=xstart;
   	tout1(k)=testt(idx,2);
    
   
end

figure(2)
imagesc(x,t,S)
set(gca,'YDir','normal')
colormap(gray), hold on 
plot(xout2,tout,'o',xout1,tout1,'*r')
hold off

p=polyfit(tout1,xout1,1);
pause

figure(3)
plot(tout1,xout1,'o',tout1,p(1)*tout1+p(2),'k')

fprintf('Peak velocity = %d \n',p(1))

%% For saving
if saving==1

Result=struct;

%% save
nr=nr+1
Result.(['line_' int2str(nr)])={[tout1;xout1];p(1)};
save([data_out '_Lines'], 'Result', '-v7.3');

%% load
load([data_out '_Lines']);  

%% create overlay
close all
fnames=fieldnames(Result);
velocities=zeros(numel(fnames),1);
figure(2)
imagesc(x,t,S)
set(gca,'YDir','normal')
colormap(gray), hold on 
for index=1:numel(fnames)
line=getfield(Result, fnames{index});
plot(line{1}(2,:),line{1}(1,:),'*r')
velocities(index)=line{2};
end
hold off
avg_v=mean(velocities);

savefig(figure(2), [data_out '_Lines.fig']);

% dlmread(liste{2})
% 
% for index=1:numel(fnames)
%     
% fprintf( 'a field:')
% fprintf( fnames{index})
% fprintf( '\n')
% end

end
% end
