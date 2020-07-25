
clc
clear
close all
warning('off','all')

% Number of clusters
M = 6

% Load data file
load("data.mat")

% Initialize variables from dataset
x = data;
[N, d] = size(x);
mean = zeros(M,2);
pi_val = zeros(1,M);
pi_val(:) = 1/M;

% Plot data
scatter(data(:,1),data(:,2),3)
hold on


% Determine initial guess using K-means method
[idx1,C1] = kmeans(data,M);

%getting the mean for each point.
for i = 1:M
  mu(i,:) = C1(i,:);
end

% Plot kmeans values
scatter(C1(:,1),C1(:,2),100)

% Initialize covariance matrices for each Gaussian (M) scalable
for m = 1:M
   Em(:,:,m) = eye(2);
end 

% Expectation Step (E)
%cond_prob = zeros(N,M);
cond_prob_total = 0;
cond_sum = zeros(18000,1);
for g = 1:15
for i = 1:N
    
    for m=1:M
        cond_prob(i,m) = mvnpdf(x(i,:),mu(m,:),Em(:,:,m));
    end
    
    for m =1:M
        cond_sum(i) = cond_prob(i,m)*pi_val(m) + cond_sum(i);
    end

    for m=1:M
        z(i,m) = cond_prob(i,m)*pi_val(m)/cond_sum(i);
    end
end

for i = 1:N
    for m = 1:M
        if isnan(z(i,m))
            z(i,m) = 0;
        end
    end
end
%% M step
cov_num = zeros(2,2,6);
for m = 1:M
    z_sum = sum(z(:,m));
    for i = 1:N
    cov_num(:,:,m)  = z(i,m) * (x(i,:) - mu(m,:))' * (x(i,:) - mu(m,:)) + cov_num(:,:,m); 
    end
    Em(:,:,m) = cov_num(:,:,m)/sum(z(:,m));
    
    for i = 1:N
      mu_num(i,:) = z(i,m) * x(i,:);
    end
    sum_num = sum(mu_num);
    mu(m,:) = sum_num/(z_sum);
    
    pi_val(m) = z_sum/N;
end
scatter(mu(:,1),mu(:,2),100)
end