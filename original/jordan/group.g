# Faithful affine action of the verified regular group on its 256-element normal subgroup.
basis:=[5,14,1,3,99,134,254,169];;
vec:=GF(2)^8;; base:=Basis(vec);;
vals:=List([0..255],i->Sum([1..8],j->(QuoInt(i,2^(j-1)) mod 2)*base[j]));;
R:=List(basis,j->vals[j+1]);;
linear:=PermList(List(vals,v->Position(vals,List(v)*R)));;
gens:=List(base,e->PermList(List(vals,v->Position(vals,v+e))));;
G:=Group(Concatenation(gens,[linear]));;
Print("order ",Size(G),"  SmallGroup ",IdGroup(G),"  center ",Size(Center(G)),"\n");
QUIT;
