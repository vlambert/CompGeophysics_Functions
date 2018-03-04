% Rank 2 identity in voigt notation
function [I] = voigtI4()
I = [ 1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 0.5 0 0;
      0 0 0 0 0.5 0;
      0 0 0 0 0 0.5];
end