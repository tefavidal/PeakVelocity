function [result, ts]=find_peak(tout,xout2, xout1,Jump,Wiggle, xstart, dx, k, tout1, dt)

    testt=zeros(Jump,2);
    
    for dxs=1:Jump 
        xtest=xstart+dx*dxs;
        
    if isempty(tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2))
        result=-1;
        break
    end
    
    testrange=tout(xout2<=(xtest+Wiggle) & (xtest-Wiggle)<=xout2);
    
    [~, testt(dxs,2)]=min(abs(testrange-tout1(k-1)));
    testt(dxs,2)=testrange(testt(dxs,2));
    testt(dxs,1)=sqrt((100*(testt(dxs,2)-xstart)/dx)^2+((testt(dxs,2)-tout1(k-1))/dt)^2);
   
    end
    
    while 1        
        [~, idx]=min(testt(:,1));
        if abs(tout1(k-1)-testt(idx,2))>3 
            testt(:,:)=[];
        else
            break
        end
        if size(testt,1)==0
            idx=-1;
            break
        end
    end
     result=idx;
     ts=testt;
    
end

    
