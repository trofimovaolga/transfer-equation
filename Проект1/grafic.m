function g = grafic(n)
    h = 3.5/20;
    X = [0:h:(3.5-2*h)];
    Y = dlmread('result.txt',  ';')
    length(Y)
    length(X)
    plot(X, Y);
end