from sympy  import *

eta_cst, Eii, eta_ve, eta_el = symbols('eta_cst, Eii, eta_ve, eta_el')
C_pwl, n_pwl = symbols('C_pwl, n_pwl')
C_exp, n_exp, ST = symbols('C_exp, n_exp, ST')
C_lin, n_lin, m_lin, d1 = symbols('C_lin, n_lin, m_lin, d1')
G, gam, pg, lam, cg= symbols('G, gam, pg, lam, cg')
Tii          = 2.0 * eta_ve * Eii
Eii_cst      = Tii/2.0/eta_cst
Eii_pwl      = C_pwl * pow(Tii, n_pwl    );
Eii_exp      = C_exp * pow(Tii, ST+n_exp );
f            = Eii_cst + Eii_pwl + Eii_exp + Eii_lin;
Eii_lin      = C_lin * pow(Tii, n_lin) * pow(d1,-m_lin);

f.diff(eta_ve).subs(d1,'d1').subs(Eii_pwl, 'Eii_pwl').subs(Eii_exp, 'Eii_exp').subs(Eii_lin, 'Eii_lin')


eta_pwl     = (2*C_pwl)**(-1)*Tii**(1-n_pwl)
eta_exp     = (2*C_exp)**(-1)*Tii**(1-(ST+n_exp))
# eta_lin     = (2*C_lin)**(-1))_* pow( Eii, 1.0/n_lin - 1.0 ) * pow(d, m_lin/n_lin);
display(eta_pwl)
eta_ve  = (eta_el**-1 + eta_pwl**-1 + eta_exp**-1 + eta_lin**-1)**-1
Exx,Ezz,Exz=symbols('Exx,Ezz,Exz')
Txx,Tzz,Txz=symbols('Txx,Tzz,Txz')
Tii = sqrt(1/2*Txx**2 + 1/2*Tzz**2 + Txz**2)

print(eta_ve)

detadTxx = eta_ve.diff(Txx).subs(eta_ve, 'eta_ve').subs(Tii, 'Tii').subs(Tii**2, 'Jii')
print(detadTxx)

# d1           = exp(log( G *gam/(lam*(1.0/cg)* Tii *(Eii_pwl + Eii_exp)*pg))/(1.0+pg))
