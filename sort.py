L = [[1,'A'],[10,'A'],[4,'B'],[8,'A'],[5,'B'],[12,'A'],[2,'A'],[6,'A'],[11,'A'],[3,'A'],[9,'A'],[7,'A']]

def triFusion(L):
    if len(L) == 1:
        return L
    else:
        return fusion( triFusion(L[:len(L)//2]) , triFusion(L[len(L)//2:]) )
    
def fusion(A,B):
    if len(A) == 0:
        return B
    elif len(B) == 0:
        return A
    elif A[0] <= B[0]:
        return [A[0]] + fusion( A[1:] , B )
    else:
        return [B[0]] + fusion( A , B[1:] )