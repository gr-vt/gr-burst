function y = qint(alpha,beta,gamma)
    p = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma);
    y = beta - 1/4*(alpha-gamma)*p;
end