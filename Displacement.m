% Evalustion of Displacements at each node

function Displacements = Displacement(p)
% Active Dof = subtract Prescribed Dof from Set of all Dof 
active_Dof=setdiff([1:p.Num_Nodes]', [p.Prescribed_Dof]);

% Modify forces at active nodes to account for homogenous/non-homogenous BC
p.Force_Modified=p.Force(active_Dof)-p.Stiffness(active_Dof,p.Prescribed_Dof)*p.Displacement(p.Prescribed_Dof);

%Solution U = Inverse of modified[K] x F for Active Dof
U=p.Stiffness(active_Dof,active_Dof)\p.Force_Modified;  

%Substituting Displacement at Active_Dof
p.Displacement(active_Dof)=U;

%Returning Displacement obtained at all Dof
Displacements = p.Displacement;
end