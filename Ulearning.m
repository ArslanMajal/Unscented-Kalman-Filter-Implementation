%This is my implementation for the Unscented Kalman Filter 
%Please try to implement your own code and use this only for support
%Our Project involves Tracking a Person through a path (Our Ground Truth)
%We have data coming from four anchors which give us the range for the 
%Tag being tracked and we are estimating the co-ordinates of the tag
%This code wont Run For you becasue i am not including the data in the 
%inidat but the this code works believe me
%

q = 0.1; 
um = [0;0];
R = 0.64*eye(4);
n = 2;
Q = q*eye(2);
P = Q; %Process Covaraince
inidat;% REMOVE THIS, i am Importing my data from this file
y = C; %This is a matrix that contains all the ranges coming from the anchors
Mt = zeros(size(um,1),length(y)); %this matrix contains all the measurements we collect
Con = zeros(size(P,1),size(P,2),length(y));%this matrix keeps track of the covariance we calculate in each step


a = 1; %alpha for lambda
kappa = 0; %for later calculations
betta = 1;
lambda = (a^2)*(n+kappa)-n;
N = size(M,1); %for the number of steps in our Experiment
for t = 1:N

U = chol(P,'lower'); %cholsky matrix for the algorithm
X = zeros(n,(2*n+1)); %for the sigma points

X(:,1) = um;
for i = 2:n+1 %this forms a total of n points + the first point equal to m
    
    X(:,i) = um + sqrt(n + lambda)*U(:,i-1);

end

for i  = 1:n % this gives us n further points making a totla of 2n + 1
    
    X(:,i + n + 1) = um - sqrt(n + lambda)*U(:,i);
    
end




%%%---------------------------THE SIGMA POINTS ARE CORRECT----------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    Wm = zeros(2*n+1,1);
    Wc = zeros(2*n+1,1);
    
    Wm(1) = lambda/(n + lambda);
    Wc(1) = Wm(1) + (1-a^2 + betta);
    
for i = 2:2*n + 1    
    Wm(i) = 1/(2*(n + lambda));
    Wc(i) = 1/(2*(n + lambda));
end

Wm;
Wc;



   
    
%%%%----------the weights are also correct---------%%%

%in this entore portion we are taking data from the Matrix M 
    K1 = M(t,:);
    A1 = K1(1:3); %setting the Co-ordinates of the anchors
    A2 = K1(5:7);
    A3 = K1(9:11);
    A4 = K1(13:15);
    An=[A1; A2; A3;A4]; %vector containing the co-ordinates of the vectors
    
    P = real(P);
  
    
    
U1 = zeros(2*n+1,1);
U2 = zeros(2*n+1,1);


for i = 1:2*n+1
    U1(i) = Wm(i)*X(1,i);
    U2(i) = Wm(i)*X(2,i);
    
end
um = [sum(U1);sum(U2)];




G = zeros(size(X,1));

for i = 1:2*n+1
    G = G + Wc(i)*(X(:,i)-um)*(X(:,i)-um)' ;
end

G = G + Q;
G = real(G);
P = G;


U = chol(P,'lower');
X = zeros(n,(2*n+1));
X(:,1) = um; %setting the first values equalt to the mean we just calculated
for i = 2:n+1 %this forms a total of n points + the first point equal to m
    
    X(:,i) = um + sqrt(n + lambda)*U(:,i-1);

end

for i  = 1:n % this gives us n further points making a totla of 2n + 1
    
    X(:,i + n + 1) = um - sqrt(n + lambda)*U(:,i);
    
end


%%%----------Perfect till the second calculation of sigma
%%%points-----------%%%


       
Y = zeros(4 , size(X,2));


for i = 1 : size(X,2) % this will equal 2n + 1
    % we wil pass it through the equation for our non-linear function y and
    % store the points in a matrix
    Y(1:4,i)=[sqrt((X(1,i)-An(1,1))^2 + (X(2,i)-An(1,2))^2 + (0.84-An(1,3))^2);...
                    sqrt((X(1,i)-An(2,1))^2 + (X(2,i)-An(2,2))^2 + (0.84-An(2,3))^2);...
                    sqrt((X(1,i)-An(3,1))^2 + (X(2,i)-An(3,2))^2 + (0.84-An(3,3))^2);...
                    sqrt((X(1,i)-An(4,1))^2 + (X(2,i)-An(4,2))^2 + (0.84-An(4,3))^2)] ;
    %where An are the co-ordinates for our Anchor matrix (this is our non-linear function)
    % you just need to write anything like y(i) = g(X(i)) for each sigma point and
    % store it in the matrix Y
end






U1 = zeros(2*n+1,1);
U2 = zeros(2*n+1,1);
U3 = zeros(2*n+1,1);
U4 = zeros(2*n+1,1);

for i = 1:2*n+1
    U1(i) = Wm(i)*Y(1,i);
    U2(i) = Wm(i)*Y(2,i);
    U3(i) = Wm(i)*Y(3,i);
    U4(i) = Wm(i)*Y(4,i);
end
u = [sum(U1);sum(U2);sum(U3);sum(U4)];





        
S = zeros(size(Y,1));
for i = 1:2*n+1
    S = S + Wc(i)*(Y(:,i)-u)*(Y(:,i)-u)';
end
S = S + R;




C  = zeros(size(X,1),size(Y,1));

for i = 1:2*n+1
    C = C + Wc(i)*(X(:,i)-um)*(Y(:,i)-u)' ;
end








K = C/S; %filter gain
um = um + K*(y(:,t) - u); %new mean for the X
P = P - K*S*K';


        

Mt(:,t) = um;
Con(:,:,t) = P;



end

    X2 = H(:,1);
    Y2 = H(:,2);
    figure
    subplot(2,1,1)
    plot(X2(1:59),'r-')
    hold on 
    plot(Mt(1,:),'o')
    legend('Ground Truth','UKF')
    title('UKF estimate x-dimension');
    
    
    subplot(2,1,2)
    plot(Y2(1:59),'r-')
    hold on 
    plot(Mt(2,:),'o')
    legend('Ground Truth','UKF')
    title('UKF estimate y-dimension');
    















