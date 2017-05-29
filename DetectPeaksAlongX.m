function [xout2,tout]=DetectPeaksAlongX(x,t,S,cut)
%[xout2,tout]=DetectPeaksAlongX(x,t,S,cut)
% Detect the minimum along each x, when the valley goes below cut

Nx=length(S(1,:));
T=length(S(:,1));
g0=cut*ones(T,Nx);
S2=(S-g0)<0;
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
          if j1+aux2<=T && aux2>1
            xout2(k)=x(i);
          
          %%%%%%%%
            laux=t((j1+aux2-2):(j1+aux2));
            laux2=S((j1+aux2-2):(j1+aux2),i);
            lauxp=polyfit(laux,laux2,2);
            %%%% p(1)x*x + p(2)*x +p(3)
            tout(k)=-lauxp(2)/(2*lauxp(1));
          %%%%%
          k=k+1;
          end
        end
        j1=j2;
    end
end