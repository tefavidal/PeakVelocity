
%% read file
close all;
clear all;

dir='20170508-New-S2h36m-6mmpmin';   
name ='Filtered_top-270-800';

path=['/data.lfpn/eckstein/Torsten/' dir '/'];
filename = [path name '.tif']

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
%g0=0.5309*ones(T,Nx);
Nx=length(S(1,:));
T=length(S(:,1));
g0=zeros(T,Nx);
S2=(S-g0)<0;
out=NaN(size(S));
k=1;
clearvars tout
clearvars xout
clearvars xout2
for i=1:Nx
    j1=1;
    j2=1;
    while j1<T && j2<T
        while j1< T && S2(j1,i)==0 
            j1=j1+1;
        end
        %%% j1 marks the first one
        j2=j1;
        while j2< T && S2(j2,i)==1
            j2=j2+1;
        end
        %%% j2 mark the first zero
        if j2~= j1
            [aux1,aux2]=min(S(j1:j2,i));
          out(j1+aux2-1,i)=aux1;
          if j1+aux2<=T && aux2>1
            xout2(k)=x(i);
          %xout(k)=x(j1+aux2-1);
          
          %%%%%%%%
            laux=t((j1+aux2-2):(j1+aux2));
            laux2=S((j1+aux2-2):(j1+aux2),i);
            %laux2=laux2';
            lauxp=polyfit(laux,laux2,2);
            %%%% p(1)x*x + p(2)*x +p(3)
            tout(k)=-lauxp(2)/(2*lauxp(1));
            %out2(i)=x(j);
          %%%%%
          k=k+1;
          end
        end
        j1=j2;
    end
end
save(data_out,'-v7.3');

    case 'Load data'
          calculate=0; 
           fprintf('DO NOT Calculate = %d \n',calculate)
          load(data_out);       
   end
else

%g0=0.5309*ones(T,Nx);
Nx=length(S(1,:));
T=length(S(:,1));
g0=zeros(T,Nx);
S2=(S-g0)<0;
out=NaN(size(S));
k=1;
clearvars tout
clearvars xout
clearvars xout2
for i=1:Nx
    j1=1;
    j2=1;
    while j1<T && j2<T
        while j1< T && S2(j1,i)==0 
            j1=j1+1;
        end
        %%% j1 marks the first one
        j2=j1;
        while j2< T && S2(j2,i)==1
            j2=j2+1;
        end
        %%% j2 mark the first zero
        if j2~= j1
            [aux1,aux2]=min(S(j1:j2,i));
          out(j1+aux2-1,i)=aux1;
          if j1+aux2<=T && aux2>1
            xout2(k)=x(i);
          %xout(k)=x(j1+aux2-1);
          
          %%%%%%%%
            laux=t((j1+aux2-2):(j1+aux2));
            laux2=S((j1+aux2-2):(j1+aux2),i);
            %laux2=laux2';
            lauxp=polyfit(laux,laux2,2);
            %%%% p(1)x*x + p(2)*x +p(3)
            tout(k)=-lauxp(2)/(2*lauxp(1));
            %out2(i)=x(j);
          %%%%%
          k=k+1;
          end
        end
        j1=j2;
    end
end

%         for i=1:T
%             plot(x,S(i,:),'o-',x,out(i,:),'*r')
%             %axis([0 6 0 0.1])
%             %axis([0 Lx 0 7])
%             title(num2str(t(i)))
%             pause
%         end

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


% while 1
%%%%% For V>2.0
%tstart=15;
%%%%% For V<=2.0
prompt = 'Start t? ';
tstart = input(prompt);

forward = 1;


k=1;
tout1(k)=tstart;


prompt = 'Inital X? ';
aux = input(prompt);
if max(size(aux))==0
        xout1(k)=min(xout2(tout<=(tstart+0.002) & (tstart-0.002)<=tout));
else
    xstart=aux;
    xout1(k)=aux;
end


Wiggle=0.003;

while tstart<t(T-1) && xout1(k) < Lx && xout1 (k) > dx
    k=k+1;
      
%     while 1
%         
%     testt=zeros(3,2);
%     
%     for dxs=1:Jump 
%         xtest=xstart+dx*dxs;
%         
%     if dxs==3 && isempty(tout(xout2<=(xtest+0.002) & (xtest-0.002)<=xout2))
%         xout1=xout1(1:(k-1));
%         break
%     end
%     
%     testrange=tout(xout2<=(xtest+0.002) & (xtest-0.002)<=xout2);
%     
%     [~, testt(dxs,2)]=min(abs(testrange-tout1(k-1)));
%     testt(dxs,2)=testrange(testt(dxs,2));
%     testt(dxs,1)=sqrt(((testt(dxs,2)-xstart)/dx)^2+((testt(dxs,2)-tout1(k-1))/dt)^2);
%     
%     end
%     
%     end
%     
%     [~, idx]=min(testt(:,1));
        
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
      

%     Jump=1;
% 	result=0;
%     testt=zeros(Jump,2);
%     for dxs=1:Jump 
%         xtest=xstart+dx*dxs;        
%     if dxs==Jump && isempty(tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2))
%         result=-1;
%         break        
%     end    
%     testrange=tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2);    
%     [~, testt(dxs,2)]=min(abs(testrange-tout1(k-1)));
%     testt(dxs,2)=testrange(testt(dxs,2));
%     testt(dxs,1)=sqrt(((testt(dxs,2)-xstart)/dx)^2+((testt(dxs,2)-tout1(k-1))/dt)^2);   
%     end    
%     [~, idx]=min(testt(:,1));     
% 
%     if result==-1
%         Jump=Jump+1;
%     testt=zeros(Jump,2);    
%     for dxs=1:Jump 
%         xtest=xstart+dx*dxs;        
%     if dxs==Jump && isempty(tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2))
%         result=-1;
%         break
%     end    
%     testrange=tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2);    
%     [~, testt(dxs,2)]=min(abs(testrange-tout1(k-1)));
%     testt(dxs,2)=testrange(testt(dxs,2));
%     testt(dxs,1)=sqrt(((testt(dxs,2)-xstart)/dx)^2+((testt(dxs,2)-tout1(k-1))/dt)^2);   
%     end    
%     [~, idx]=min(testt(:,1));  
%         if result ==-1
%         Jump=Jump+1;
%       testt=zeros(Jump,2);    
%     for dxs=1:Jump 
%         xtest=xstart+dx*dxs;        
%     if dxs==Jump && isempty(tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2))
%         result=-1;
%         break
%     end    
%     testrange=tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2);    
%     [~, testt(dxs,2)]=min(abs(testrange-tout1(k-1)));
%     testt(dxs,2)=testrange(testt(dxs,2));
%     testt(dxs,1)=sqrt(((testt(dxs,2)-xstart)/dx)^2+((testt(dxs,2)-tout1(k-1))/dt)^2);   
%     end    
%     [~, idx]=min(testt(:,1));  
%         end
%     end
           
%     [idx, testt]=result;
%       
      
%     tout1(k)=testt(idx,1);
%     j=testt(idx,2);
%     [tout1(k), j]=testt(:,1);

        if idx==-1 || abs(tout1(k-1)-testt(idx,2))>3 || Jump>10
            xout1=xout1(1:(k-1));
            tout1=tout1(1:(k-1));
            break
        end    
                    
    xstart=xstart+dx*idx;
    xout1(k)=xstart;
   	tout1(k)=testt(idx,2);
    
%     aux1=tout(xout2<=(xtest+0.002) & (xtest-0.002)<=xout2);
% 	tout1(k)=aux1(j);
   
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
