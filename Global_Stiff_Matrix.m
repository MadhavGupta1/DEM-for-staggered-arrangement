
function Global_Matrix = Global_Stiff_Matrix(p)
   for e=1:p.Num_Elements 
       %Nodes of each Element 
       el_node = p.Element_Nodes(e, 1:2);
       %x-coordinate of Node
       node_xx = p.Node_coordinate(el_node);  
       %Element_Stiff_Matrix: Element stiffness Matrix
       Element_Stiff_Matrix=Element_Stiffness(node_xx, p.Element_E(e), p.Element_Area(e));
       %Element_Dof: Nodes for the element   
       temp = p.NDof*(el_node(1)-1)+1:p.NDof*(el_node(1)-1)+p.NDof; 
       Element_Dof = [temp p.NDof*(el_node(2)-1)+1:p.NDof*(el_node(2)-1)+p.NDof];
       p.Stiffness(Element_Dof,Element_Dof) = p.Stiffness(Element_Dof,Element_Dof) + Element_Stiff_Matrix; 
   end 
 Global_Matrix = p.Stiffness; 
end
