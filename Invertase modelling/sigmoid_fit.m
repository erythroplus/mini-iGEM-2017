function sigmoid_fit

data = csvread('/PATH/TO/SOME/FOLDER/phic31.csv');

t = data(1:1000:end,1);
s = data(1:1000:end,2)./100;


k = 50; n = 3;
err = mse(s,sigmoid(t,k,n));

for i=1:1000000
knew = k + random('Uniform',-1,1)/1000;
nnew = n + random('Uniform',-1,1)/1000;

errnew = mse(s,sigmoid(t,knew,nnew));

if errnew < err
    k = knew;
    n = nnew;
    err = errnew;
end

end

fprintf('The MSE for k=%0.4f and n=%0.4f is %0.4f\n',k,n,err);

end


function s = sigmoid(x,k,n)
s = (x.^n)./((k^n)+(x.^n));
end

function err = mse(x,y)
    err = mean((x-y).^2);
end
