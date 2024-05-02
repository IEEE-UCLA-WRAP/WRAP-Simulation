% 17 Nov 2021, EA
% Erdal Arikan reserves all rights to this code.

function x = polar_encode(u)
    N = size(u,2); % N must be a power of 2
    n = log2(N);
    if n==1
        x = [mod(u(:,1)+u(:,2),2) u(:,2)];
        return;
    else
        x1 = polar_encode(mod(u(:,1:N/2)+u(:,N/2+1:N),2));
        x2 = polar_encode(u(:,N/2+1:N));
        x = [x1 x2];
    end
return
 
