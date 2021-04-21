y=[0 1 1 1 3 4];
r=[0 0 1 1 2 2;0 1 0 1 1 2];
funcs={@(x)x(1); @(x)x(2); @(x)x(1).^2; @(x)x(2).^2; @(x)x(1).*x(2)};
[P, sgP]=LinApproximator(y, r, funcs)

y=[2 3 1.368 2.368 2.135 3.135];
r=[0 0 1 1 2 2;0 1 0 1 1 2];
fun = @(r,P)P(1).*exp(-P(2).*r(1))+P(3)+r(2);
[P, sgP]=NonLinApproximator(y, r, fun, 3)