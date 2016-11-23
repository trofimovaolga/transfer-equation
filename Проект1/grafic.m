function g = grafic(n)
    X = [0:0.1:0.8];
    Y = dlmread('result.txt',  ';')
    length(Y)
    length(X)
    plot(X, Y);
end