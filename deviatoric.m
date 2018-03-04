% Function solves for deviatoric component of a 2nd-order
% tensor provided in voigt notation
function [dev] = deviatoric(A)
dev = A - 1/3*sum(A(1:3))*voigtI2;
end