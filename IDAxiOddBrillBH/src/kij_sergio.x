        o1 = a**2
        o2 = a*o1
        o3 = -2.00000000000000d0*mass*rBL
        o4 = rBL**2
        o5 = o1 + o3 + o4
        o6 = sqrt(o5)
        o7 = cos(phi(i,j,k))
        o8 = cos(q(i,j,k))
        o9 = o8**2
        o10 = o1*o9
        o11 = o10 + o4
        o12 = o11**2
        o13 = 1/o12
        o14 = -3.00000000000000d0*eta(i,j,k)
        o15 = exp(o14)
        o16 = sin(phi(i,j,k))
        o17 = sin(q(i,j,k))
        o18 = o17**2
        o19 = -o1
        o20 = mass**2
        o21 = o19 + o20
        o22 = sqrt(o21)
        o23 = 1/o4
        o24 = 2.50000000000000d-1*o21*o23
        o25 = 1.00000000000000d0 + o24
        o26 = o1 + o4
        o27 = 2.00000000000000d0*o26*o4
        o28 = o19 + o4
        o29 = o11*o28
        o30 = o27 + o29
        o31 = -2.00000000000000d0*eta(i,j,k)
        o32 = exp(o31)
        o33 = o7**2
        o34 = o16**2
        o35 = o17*o18
        kxx(i,j,k) = -(a*mass*o13*o16*o18*o22*o25*o30*o32*o7) + 4.000000
     &  0000000d0*mass*o13*o15*o16*o18*o2*o6*o7*o9*rBL
        kxy(i,j,k) = 5.0000000000000d-1*a*mass*o13*o18*o22*o25*o30*o32*o
     &  33 - 5.0000000000000d-1*a*mass*o13*o18*o22*o25*o30*o32*o34 - 2.0
     &  0000000000000d0*mass*o13*o15*o18*o2*o33*o6*o9*rBL + 2.0000000000
     &  0000d0*mass*o13*o15*o18*o2*o34*o6*o9*rBL
        kxz(i,j,k) = -5.0000000000000d-1*a*mass*o13*o16*o17*o22*o25*o30*
     &  o32*o8 - 2.00000000000000d0*mass*o13*o15*o16*o2*o35*o6*o8*rBL
        kyy(i,j,k) = a*mass*o13*o16*o18*o22*o25*o30*o32*o7 - 4.000000000
     &  0000d0*mass*o13*o15*o16*o18*o2*o6*o7*o9*rBL
        kyz(i,j,k) = 5.0000000000000d-1*a*mass*o13*o17*o22*o25*o30*o32*o
     &  7*o8 + 2.00000000000000d0*mass*o13*o15*o2*o35*o6*o7*o8*rBL
        kzz(i,j,k) = 0
