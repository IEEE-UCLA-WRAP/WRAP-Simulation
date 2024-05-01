% 17 Nov 2021 EA
% Erdal Arikan reserves all rights to this code.
% We thank him for sharing these codes with us for the purpose of the
% course project and should not be used otherwise.

function uh = polar_decode(Ly,FI,FV)
% Ly is a log likelihood ratio
% FI: 1xN 0-1 vector with FI(i) = 1 if i is a frozen index, 0 o.w.
% FV: 1xN 0-1 vector with FV(i) = value of frozen bit i (0 if i not frozen)
% LR is the likelihood ratio at all levels
% uh is the decoder decision vector

N= length(Ly);
n = log2(N);
uh = FV;
LR = zeros(1,2*N-1);
ALo = zeros(1,2*N-1); % initially all influences at all levels are zero (for odd indices)
ALe = zeros(1,N-1); % initially all influences at all levels are zero (for even indices)
oddeven = zeros(1,n+1);  % keeps track of which AL array to write to or read from 

LR(1,N:2*N-1) = Ly;

% Determine the decoder machine actions  
A = [1; 1; 1]; % Address variable for the arrays LR and AL
L = [1; 1; 1]; % L(t) denotes the level activated at time t; level n+1 is the channel level, level 1 is the message level
S = [1; 2; 2]; % S(t) denotes the smallest node index activated at time t at level L(t) 
G = [0; 1; 2];  % instruction type G combines D and E into a single variable
C = 3; % length of schedule for N=2


for i=2:n
    A = [2^(i-1); A ; 2^(i-1); A; 2^(i-1)];
    L = [i; L; i; L; i];
    S = [1; 2*S-1; 2; 2*S; 2];
    G = [0; G; 1; G; 2];
    C = 2*C+3;
end

for t=1:C  % scheduler time index
  if (G(t) == 0) % left move, odd index
    a1 = LR(1,2*A(t):2:4*A(t)-1);
    a2 = LR(1,2*A(t)+1:2:4*A(t)-1);
    LR(1,A(t):2*A(t)-1) = max(0,a1+a2)+log(1+exp(-abs(a1+a2)))-(max(a1,a2)+log(1+exp(-abs(a1-a2))));  % Exact 
  elseif (G(t) == 1) % left move, even index
    a1 = LR(1,2*A(t):2:4*A(t)-1);
    a2 = LR(1,2*A(t)+1:2:4*A(t)-1);
    if oddeven(L(t)) == 1
        LR(1,A(t):2*A(t)-1) =  (1-2*ALo(1,A(t):2*A(t)-1)).*a1 + a2;
    else
        LR(1,A(t):2*A(t)-1) =  (1-2*ALe(1,A(t):2*A(t)-1)).*a1 + a2;
    end
  elseif (G(t) == 2) % right move
      if oddeven(L(t)+1) == 0
        ALo(1,2*A(t)+1:2:4*A(t)-1) = ALe(1,A(t):2*A(t)-1);
        ALo(1,2*A(t):2:4*A(t)-1) = mod(ALo(1,A(t):2*A(t)-1)+ALe(1,A(t):2*A(t)-1),2);
        oddeven(L(t)+1) = 1;
      else
        ALe(1,2*A(t)+1:2:4*A(t)-1) = ALe(1,A(t):2*A(t)-1);
        ALe(1,2*A(t):2:4*A(t)-1) = mod(ALo(1,A(t):2*A(t)-1)+ALe(1,A(t):2*A(t)-1),2);
        oddeven(L(t)+1) = 0;
      end    
  end
  if (A(t) == 1) && (G(t) < 2) % decision level
      if FI(S(t)) == 0
          uh(S(t)) = (LR(1,1) <=0);
      end    
      if oddeven(1) == 0
          ALo(1,1) = uh(S(t));
          oddeven(1) = 1;
      else
          ALe(1,1) = uh(S(t));  
          oddeven(1) = 0;
      end    

  end
end

