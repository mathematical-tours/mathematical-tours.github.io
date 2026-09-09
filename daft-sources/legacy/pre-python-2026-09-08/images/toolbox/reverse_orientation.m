function face1 = reverse_orientation(face)

nface = size(face);
face1 = face;

for i=1:nface
    f = face1(i,:);
    face1(i,:) = f(3:-1:1);
end
