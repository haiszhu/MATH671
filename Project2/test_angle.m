 
k = 7;
x = [3,0,0,-3,2,2,-2];
y = [0,3,-3,0,0,0,0];
z = [0,0,0,0,1,-1,1];
vx = atan(y(k)/x(k)); 
        if x(k)<0,
            vx = vx + pi; 
        else if y(k)<0, 
                 vx = vx + 2*pi; 
             end
        end
        ux = asin(z(k));
        if sqrt(x(k)^2+y(k)^2)-2 < 0
            ux = ux + pi;
        else if z(k)<0,
                ux = ux + 2*pi;
            end
        end        
        uvprime = [ux,vx];
