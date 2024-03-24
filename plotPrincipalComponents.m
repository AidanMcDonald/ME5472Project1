function plotPrincipalComponents(A, gridCoords)
X = [];
Y = [];
U = [];
V = [];

for i=1:size(A,1)
    for j=1:size(A,2)   
        [eigvecs, lamda] = eig(reshape(A(i,j,:,:),[2,2]));
        [d, ind] = sort(diag(lamda));
        lamdaSorted = lamda(ind,ind);
        lamdas(i,j) = lamdaSorted(1,1);
        VSorted = eigvecs(:,ind);     
        
        V1 = lamdaSorted(1,1)*VSorted(:,1);
        U = [U,V1(1)];
        V = [V,V1(2)];
        X = [X,gridCoords(i,j,1)];
        Y = [Y,gridCoords(i,j,2)];

    end
end
figure()
quiver(X,Y,U,V)
title("First principal strain vector")
xlabel("x (pixels)")
ylabel("y (pixels)")

figure()
title("First principal strain magnitude")
contourf(gridCoords(:,1,1),gridCoords(1,:,2),lamdas, "ShowText",true)