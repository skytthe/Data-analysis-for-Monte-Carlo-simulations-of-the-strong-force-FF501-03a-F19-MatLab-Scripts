function [b, a, chi] = linRegress(x, y)
    xa = sum(x)/length(x);
    ya = sum(y)/length(y);
    ssxx = sum((x-xa).^2);
    ssxy = sum((x-xa).*(y-ya));
    
    b = ssxy/ssxx;
    a = ya - b*xa;
    
    chi = sum( (y - (b*x+a)).^2 );
end