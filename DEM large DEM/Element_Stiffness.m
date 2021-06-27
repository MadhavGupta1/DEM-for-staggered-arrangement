% Evaluation of Truss Element Stifffness matrix

function Element_Stiff_Matrix = Element_Stiffness(L, E, A)
% This function evaluates 2x2 element stiffness Matrix
% This Matrix must be in GLOBAL Coordinates
Le=L;
% Element Stiffness Matrix
Element_Stiff_Matrix = ((E*A)/(L))*[1 -1 ; -1 1];
end


