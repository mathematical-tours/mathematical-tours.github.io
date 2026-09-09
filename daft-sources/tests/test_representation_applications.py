import sys,itertools,math,unittest
from pathlib import Path
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'listings'))
from invariants import matrix_group,exponents_of_degree,reynolds_generators

class RepresentationApplicationChecks(unittest.TestCase):
    def test_character_columns_need_conjugation(self):
        for n in [3,4,7]:
            T=np.exp(2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)
            np.testing.assert_allclose(T.conj().T@T,n*np.eye(n),atol=2e-14)
        T=np.array([[1,1,1,1,1],[1,-1,1,-1,1],[3,1,0,-1,-1],[3,-1,0,1,-1],[2,0,-1,0,2]])
        size=np.array([1,6,8,6,3]);U=T*np.sqrt(size/24)
        np.testing.assert_allclose(U@U.T,np.eye(5),atol=3e-16)
        for values,multiplicity in [([6,0,0,2,2],[1,0,0,1,1]),([8,0,2,0,0],[1,1,1,1,0]),([12,2,0,0,0],[1,0,2,1,1])]:
            np.testing.assert_array_equal((T*size)@values//24,multiplicity)

    def test_s3_survey_exact_reconstruction(self):
        order=[(1,2,3),(2,1,3),(3,2,1),(1,3,2),(2,3,1),(3,1,2)]
        B=sp.Matrix([[1/sp.sqrt(2),1/sp.sqrt(6)],[-1/sp.sqrt(2),1/sp.sqrt(6)],[0,-2/sp.sqrt(6)]])
        P=[sp.eye(3)[:,[i-1 for i in p]] for p in order]
        R=[(B.T*M*B).applyfunc(sp.simplify) for M in P]
        for i in range(6):
            for j in range(6):
                product=tuple(order[i][order[j][k]-1] for k in range(3))
                self.assertEqual((R[i]*R[j]-R[order.index(product)]).applyfunc(sp.simplify),sp.zeros(2))
        basis=sp.Matrix([[1]*6,[M.det() for M in P]]+[[M[i,j] for M in R] for i,j in itertools.product(range(2),repeat=2)])
        weights=sp.diag(1,1,2,2,2,2)
        f=sp.Matrix([242,170,359,28,12,628]);coefficients=weights*basis*f
        self.assertEqual(coefficients,sp.Matrix([1439,325,-109,285*sp.sqrt(3),-947*sp.sqrt(3),-203]))
        self.assertEqual((basis.T*coefficients/6-f).applyfunc(sp.simplify),sp.zeros(6,1))
        self.assertEqual((basis*basis.T).applyfunc(sp.simplify),6*weights.inv())
        self.assertEqual([sum(int(f[i]) for i,p in enumerate(order) if p[j]==3) for j in range(3)],[987,40,412])

    def test_quaternion_matrices_and_basis(self):
        I=sp.eye(2);i=sp.diag(sp.I,-sp.I);j=sp.Matrix([[0,-1],[1,0]]);k=i*j
        self.assertEqual(i*i,j*j);self.assertEqual(i*i,-I);self.assertEqual(k*k,-I)
        self.assertEqual(j*k,i);self.assertEqual(k*i,j);self.assertEqual(j*i,-k)
        group=[I,-I,i,-i,j,-j,k,-k]
        for A in group:self.assertEqual(A.conjugate().T*A,I)
        chars=sp.Matrix([[1]*8,[1,1,1,1,-1,-1,-1,-1],[1,1,-1,-1,1,1,-1,-1],[1,1,-1,-1,-1,-1,1,1]])
        entries=sp.Matrix([[sp.sqrt(2)*M[a,b] for M in group] for a,b in itertools.product(range(2),repeat=2)])
        basis=chars.col_join(entries)
        self.assertEqual((basis*basis.conjugate().T).applyfunc(sp.simplify),8*sp.eye(8))
        a,b,c,d=sp.symbols('a b c d',real=True);M=a*I+b*i+c*j+d*k
        self.assertEqual((M.conjugate().T*M).applyfunc(sp.simplify),(a*a+b*b+c*c+d*d)*I)

    def test_cube_rotation_invariants(self):
        X,Y,Z=sp.symbols('X Y Z');v=sp.Matrix([X,Y,Z]);rotations=[]
        for perm in itertools.permutations(range(3)):
            for signs in itertools.product([-1,1],repeat=3):
                A=sp.eye(3)[:,perm]*sp.diag(*signs)
                if A.det()==1:rotations.append(A)
        self.assertEqual(len(rotations),24)
        f=(X+Y+Z)*(X+Y-Z)*(X-Y+Z)*(-X+Y+Z)
        g=(X*X-Y*Y)*(X*X-Z*Z)*(Y*Y-Z*Z)
        for A in rotations:
            substitution=dict(zip(v,A*v))
            for P in [X*X+Y*Y+Z*Z,(X*Y*Z)**2,f,g*g,X*Y*Z*g]:
                self.assertEqual(sp.expand(P.subs(substitution,simultaneous=True)-P),0)

    def test_exact_reynolds_and_molien(self):
        X,Y,t=sp.symbols('X Y t');A=sp.Matrix([[1,1],[1,-1]])/sp.sqrt(2);D=sp.diag(1,-1)
        for generators,expected_order,expected_series in [([A],2,1/((1-t)*(1-t*t))),([A,D],16,1/((1-t*t)*(1-t**8)))]:
            group=matrix_group(generators);self.assertEqual(len(group),expected_order)
            molien=sum(1/(sp.eye(2)-t*B).det() for B in group)/len(group)
            self.assertEqual(sp.factor(molien-expected_series),0)
            invariants=reynolds_generators(group,[X,Y]);self.assertTrue(invariants)
            for P in invariants:
                for B in generators:
                    self.assertEqual(sp.simplify(sp.expand(P.subs(dict(zip([X,Y],B*sp.Matrix([X,Y]))),simultaneous=True)-P)),0)
                if len(group)==16:
                    degree=sp.Poly(P,X,Y).total_degree()
                    monomials=[sp.expand((X*X+Y*Y)**((degree-8*j)//2)*(X*X*Y*Y*(X*X-Y*Y)**2)**j)
                               for j in range(degree//8+1)]
                    def vector(Q):return sp.Matrix([sp.Poly(Q,X,Y).coeff_monomial(X**i*Y**(degree-i)) for i in range(degree+1)])
                    M=sp.Matrix.hstack(*map(vector,monomials));self.assertEqual(M.rank(),M.row_join(vector(P)).rank())
        P=X*X+Y*Y;Q=X*X*Y*Y*(X*X-Y*Y)**2
        self.assertEqual(sp.expand(P**4-4*Q),X**8+14*X**4*Y**4+Y**8)
        self.assertNotEqual(sp.det(sp.Matrix([[sp.diff(P,v) for v in [X,Y]],[sp.diff(Q,v) for v in [X,Y]]])),0)
        self.assertEqual(list(exponents_of_degree(2,2)),[(0,2),(1,1),(2,0)])

if __name__=='__main__':unittest.main()
