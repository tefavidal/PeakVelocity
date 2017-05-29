if ~exist('A','var')
    fprintf('Please call your data A \n');
    return
end

S=A;
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

        for i=1:T
            plot(x,S(i,:),'o-',x,out(i,:),'*r')
            %axis([0 6 0 0.1])
            %axis([0 Lx 0 7])
            title(num2str(t(i)))
            pause
        end

%% Part 2


plot(xout2,tout,'o')

clearvars tout1
clearvars xout1



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




while tstart<t(T-1) && xout1(k) < Lx && xout1 (k) > dx
    k=k+1;
    xstart=xstart+dx;
    xout1(k)=xstart;
    if isempty(tout(xout2<=(xstart+0.002) & (xstart-0.002)<=xout2))
        xout1=xout1(1:(k-1));
        break
    end
        
        [tout1(k), j]=min(abs(tout(xout2<=(xstart+0.002) & (xstart-0.002)<=xout2)-tout1(k-1)));
        
        if tout1(k)>3
            xout1=xout1(1:(k-1));
            tout1=tout1(1:(k-1));
        break
        end
        
        aux1=tout(xout2<=(xstart+0.002) & (xstart-0.002)<=xout2);
        tout1(k)=aux1(j);
   
end



plot(xout2,tout,'o',xout1,tout1,'*r')

p=polyfit(tout1,xout1,1);
pause
plot(tout1,xout1,'o',tout1,p(1)*tout1+p(2),'k')

fprintf('Peak velocity = %d \n',p(1))


