        mass = sqrt(2.00000000000000d0 + sqrt(4.0000000000000d0 + byJ**2
     &  ))
        a = byJ/mass
        rBL = 5.0000000000000d-1*exp(etagrd(i))*(1.00000000000000d0 + (1.d0*(a
     &   + mass)*exp(-etagrd(i)))/sqrt(-a**2 + mass**2))*sqrt(-a**2 + mass**2)
     &  *(1.00000000000000d0 + (1.d0*(a + mass)*exp(-etagrd(i)))/sqrt(a**2 + m
     &  ass**2))
