function [ new_mat ] = plusmat(matrix)
% PLUSMAT: crea una nuova matrice facendo la media tra le varie colonne/righe

% controllo che sia una matrice NxN
dim=size(matrix);
assert(dim(1)==dim(2),'ERROR: size of Matrix not correct!');

dime=dim(1);

% ridimensionamento della matrice - aumento delle COLONNE
mat=zeros(2*dime-1,2*dime-1);

for r=1:dime
for c=1:dime
    if c==1;
       mat(r,c)=matrix(r,c); 
    else 
        mat(r,2*c-1)=matrix(r,c);
    end
end
end

% medie colonne
for r=1:(dime)
for c=1:(dime-1)        
    mat(r,(2*c))=(matrix(r,c)+matrix(r,c+1))/2;    
end
end

% nuova matrice
mat2=zeros(2*dim-1,2*dim-1);

% spostamento delle righe
for i=1:dim
mat2((2*i-1),1:(2*dim-1))=mat(i,1:(2*dim-1));
end

% media righe
for c=1:(2*dime-1)
for r=1:(dime-1)
    mat2(2*r,c)=(mat2(2*r-1,c)+mat2(2*r+1,c))/2;
end
end

%% print new_mat
new_mat=mat2;

end

