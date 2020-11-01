% Test pour la multiplication d'entiers.
x = 262122154512154212
y = 314134464653513212
% la base
b = 20
xx = number2vector(x,b)
yy = number2vector(y,b)
zz = mult_entiers(xx,yy,b)
z = vector2number(zz,b)
% le résultat doit valoir zéro
z - x*y