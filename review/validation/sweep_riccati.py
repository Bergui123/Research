"""Parameterized exact-rational Riccati residual certificate for sweep data.

The result is conditional on the positive-definite Riccati comparison; see
sweep-riccati.md. All checked quantities are rational Bernstein bounds.
"""
from fractions import Fraction as F

from verify_riccati_nodes import bernstein, quadratic


def verify(data,min_hat,expected_initial=(F(2),F(3,10),F(1))):
    n,q=data['segments'],data['denominator']
    rho,lam=F(data['rho']),F(data['lambda'])
    assert n>0 and q>0 and 0<=rho<1 and lam>0 and min_hat>0
    assert data['components']==['P11','P12','P22'] and len(data['nodes'])==2*n+1
    nodes=[[F(v,q) for v in row] for row in data['nodes']]
    assert nodes[0]==list(expected_initial)
    c=[[F(1),rho],[rho,F(1)]]
    integral=F(0);max_row=F(0)
    for i in range(n):
        entries=[quadratic(*v) for v in zip(*nodes[2*i:2*i+3])]
        b=[bernstein(v) for v in entries]
        for k in range(3):
            a,x,d=b[0][k],b[1][k],b[2][k]
            assert x>=0 and a>min_hat and d>min_hat
            assert (a-min_hat)*(d-min_hat)>x*x
            max_row=max(max_row,a+x,d+x)
        p=[[entries[0],entries[1]],[entries[1],entries[2]]]
        residual=[]
        for a,b_index in ((0,0),(0,1),(1,1)):
            coeff=[F(0)]*5
            coeff[0]=n*p[a][b_index][1]-2*lam*c[a][b_index]
            coeff[1]=2*n*p[a][b_index][2]
            for k in (0,1):
                for j in range(3):
                    for ell in range(3):
                        coeff[j+ell]+=p[a][k][j]*p[k][b_index][ell]/2
            residual.append(max(map(abs,bernstein(coeff))))
        integral+=(residual[0]+2*residual[1]+residual[2])/n
    assert integral<min_hat
    return {'segments':n,'rho':str(rho),'lambda':str(lam),
            'interpolant_min_eigen_lower':str(min_hat),
            'interpolant_max_row_upper':str(max_row),
            'residual_integral_upper':str(integral),
            'true_min_eigen_lower':str(min_hat-integral),
            'true_row_upper':str(max_row+2*integral)}
