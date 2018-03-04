% Function solves for the norm of a matrix given as a 6-element voigt
% vector
function [norm] = voigtnorm(A)
norm = sqrt(A(1)^2 + A(2)^2 + A(3)^2 +...
       2*( (A(4)/2)^2 + (A(5)/2)^2 + (A(6)/2)^2 ));
end