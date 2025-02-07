########### Non-self-referential Gram-Schimdt Process #############
rm(list=ls()); 
##----------------
U_dot_k=function(U,k) ##U is a matrix with at least k+1 columns
{
U_k_plus1=U[,k+1];
U_1tok=U[,1:k];
u_dot_k=t(U_k_plus1)%*%U_1tok;
return(u_dot_k);
}
##-----------------
normMatrixfromV=function(V_kminus1)
{
n=dim(V_kminus1)[2];
normMatV=matrix(0,n,n);
diag(normMatV)=-1/(colSums(V_kminus1*V_kminus1));
return(normMatV);
}

##-----------------
V_k=function(b_plus_Delta_k,U,k)
{
Vk=U[,1:k]%*%t(b_plus_Delta_k);
return(Vk);
}

##------------------
B_plus_Delta_k=function(b_plus_Delta_k,U,k,Vk)
{
m=dim(b_plus_Delta_k)[1];
n=dim(b_plus_Delta_k)[2];
b_plus_Delta_k2=matrix(0,m+1,n+1);
b_plus_Delta_k2[1:m,1:n]=b_plus_Delta_k;
b_plus_Delta_k2[1,n+1]=0;
b_plus_Delta_k2[m+1,n+1]=1;
u_dot_k=U_dot_k(U,k);   
normMatV=normMatrixfromV(Vk);
a1=u_dot_k%*%t(b_plus_Delta_k);
a2=a1%*%normMatV;
a3=a2%*%b_plus_Delta_k;
b_plus_Delta_k2[m+1,1:n]=a3;
b_plus_Delta_k=b_plus_Delta_k2;
return(b_plus_Delta_k);
}

##-(main program)---Non-self-referential Gram-Schimdt Process----
##U is a m-by-n matrix consisting column vectors, where n<=m
##kk is the number of generated orthogonal vectors, where k<=n.
NsrGSP=function(U,kk)
{
b_plus_Delta_k=matrix(c(1),1,1); ##the initial B^+_{Delta_1}=[1];
Vk=as.matrix(U[,1],default=0);  ##the initial Vk=U[,1];
 for(k in 1:(kk-1))
 {
 b_plus_Delta_k=B_plus_Delta_k(b_plus_Delta_k,U,k,Vk);
 Vk=V_k(b_plus_Delta_k,U,k+1);  ##indeed this is V_{k+1}
 }
return(Vk);
}


######======= Experiment 1
U=matrix(c(-2,5,1,0,4,7,-1,-1,2,-4,10,0,-3,2,1,2,5,1,7,-1,4,3,11,0,6),5,5);
kk=dim(U)[2];
V=NsrGSP(U,kk);
check=t(V)%*%V;
print("InputU=");U;
print("Output V=NsrGSP(U)=");V;
print("check the orthogonality of V via t(V)%*%V=");t(V)%*%V;

######======== Experiment 2
U=matrix(runif(15*9),15,9);
kk=dim(U)[2];
V=NsrGSP(U,kk);
check=t(V)%*%V;
print("Input U=");U;
print("Output V=NsrGSP(U)=");V;
print("check the orthogonality of V via t(V)%*%V=");t(V)%*%V;






