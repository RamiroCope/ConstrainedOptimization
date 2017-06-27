function [Spx,Splambda] = EqualityQPSensitivity(Wxp,cp,Wxx,cx)

[m,n] = size(cx);
Sp = -[Wxp -cp]/[Wxx -cx; -cx' zeros(n,n)];
Spx = Sp(:,1:m);
Splambda = Sp(:,m+1:end);

end