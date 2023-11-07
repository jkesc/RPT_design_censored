function [r,theta]=angularTransform(X,Y,Px,Py)
    r=sqrt((X-Px).^2+(Y-Py).^2);
    theta=atan((Y-Py)./(X-Px)).*(X>=Px)+(atan((Y-Py)./(X-Px))-pi).*(X<Px);
end