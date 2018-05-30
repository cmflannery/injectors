from __future__ import division, absolute_import, print_function
import numpy as np
import warnings


class SwirlElement():
    def __init__(self, elementType, initialConditions, version=1):
        """ SwirlElement() aids in the design of coaxial swirl elements for bipropellant
        injection.

        Parameters:
            elementType: [str] (required)
                internal mixing
                external mixing
            initialConditions: [dict] (required)
                alpha1, alpha2, deltaP1, deltaP2, mdot1, mdot2, rho1, rho2
            version: [int] (optional)
                1,2
                Note: if version is not specified, calculation version 1 will be used

        Units:
            Meters, Kilograms, Degrees, Kelvin
        """
        if version:
            version = int(version)
        if version != 1 and version != 2:
            raise ValueError('version is an optional parameter and can only be 1 or 2, see Bazarov for help')
        self.version = version
        elementType = elementType.lower()
        validElementTypes = ['internal', 'external']
        if elementType not in validElementTypes:
            raise TypeError('elementType, %s, is not a valid type' % elementType)
        self.elementType = elementType
        self.checkParams(initialConditions)
        self.initialConditions = initialConditions
        # assign ICs to variables
        self.alpha1 = self.initialConditions['alpha1']
        self.alpha2 = self.initialConditions['alpha2']
        self.deltaP1 = self.initialConditions['deltaP1']
        self.deltaP2 = self.initialConditions['deltaP2']
        self.mdot1 = self.initialConditions['mdot1']
        self.mdot2 = self.initialConditions['mdot2']
        self.rho1 = self.initialConditions['rho1']
        self.rho2 = self.initialConditions['rho2']
        self.nu1 = self.initialConditions['nu1']
        self.nu2 = self.initialConditions['nu2']
        self.n1 = self.initialConditions['n1']
        self.n2 = self.initialConditions['n2']
        self.checkInitialConditions()

    def checkParams(self, initialConditions):
        """ verify that all the required initial conditions were thrown to during initalization """
        requiredInitialConditions = [
            'alpha1',   # spray cone angle (degrees)
            'alpha2',  
            'deltaP1',  # delta P (Pascals)
            'deltaP2',
            'mdot1',    # mass flow rate (kg/s)
            'mdot2',
            'rho1',     # fluid density (kg/m^3)
            'rho2',
            'nu1',      # kinematic viscosity
            'nu2',
            'n1',       # number of inlet passages
            'n2'
        ]
        missingParams = []
        for ric in requiredInitialConditions:
            if ric not in initialConditions:
                missingParams.append(ric)
        if missingParams:
            missingParamsStr = ''
            for mp in missingParams:
                missingParamsStr = missingParamsStr + '\'%s\' ' % mp
            raise TypeError('initialConditions missing %i required arguments: %s' % (len(missingParams), missingParamsStr))

    def checkInitialConditions(self):
        """ checkInitialConditions verifies that the initial conditions thrown to the
        instance are physically valid."""
        alpha1 = self.initialConditions['alpha1']
        alpha2 = self.initialConditions['alpha2']

        if (2*alpha1 - 2*alpha2 < 10) or (2*alpha1 - 2*alpha2 > 15):
            warnings.warn('Spray Cone Angles not valid. See Bazarov pg. 73')

    def design(self):
        """ Start the design process """
        if self.elementType == 'external' and self.version == 1:
            self.version1internal()

    
    def version1internal(self):
        """ Internal mixing version 1 design process """
        # determine the Geometric Characteristic Parameter for each stage
        print('Use empirical figures 34a and b for the following calculations')
        while True:
            A1 = float(input('Enter the value for A1 given 2*alpha1 = %.2f: ' % (2*self.alpha1)))
            A2 = float(input('Enter the value for A2 given 2*alpha2 = %.2f: ' % (2*self.alpha2)))

            mu1 = float(input('Enter the value for mu1 given A1 = %.2f: ' % A1))
            mu2 = float(input('Enter the value for mu2 given A2 = %.2f: ' % A2))

            # Calculate the nozzle radii
            Rn1 = self.Rn(self.mdot1,self.rho1,self.deltaP1,mu1)
            Rn2 = self.Rn(self.mdot2,self.rho2,self.deltaP2,mu2)
            Rin1 = 1.2*Rn1  # note: 1.2 is semiarbitrary, see Bazarov for further explanation
            Rin2 = 1.2*Rn2

            rin1 = self.rin(Rn1,Rin1,self.n1,A1)
            rin2 = self.rin(Rn2,Rin2,self.n2,A2)

            Rein1 = self.Rein(self.mdot1,self.n1,rin1,self.rho1,self.nu1)
            Rein2 = self.Rein(self.mdot2,self.n2,rin2,self.rho2,self.nu2)

            if Rein1 and Rein2 > 10^4: # convergence condition
                break

        print('\nRein1 = %.2f\nRein2 = %.2f' % (Rein1, Rein2))
        print('\nRn1 = %.5f cm\nRn2 = %.5f cm' % (Rn1*100, Rn2*100))
        print('\nRin1 = %.5f cm\nRin2 = %.5f cm' %(Rin1*100, Rin2*100))
        print('\nrin1 = %.5f cm\nrin2 = %.5f cm' %(rin1*100, rin2*100))

    def Rn(self, mdot, rho, deltaP, mu):
        """ Calculate the nozzle radius

            Bazarov Eq. 103
        """
        return 0.475*np.sqrt(mdot/(mu*np.sqrt(rho*deltaP)))

    def rin(self, Rn, Rin, n, A):
        """ Calculate the radius of the inlet passages
        
            Bazarov Eq. 104
        """
        return np.sqrt(Rin*Rn/(n*A))

    def Rein(self, mdot, n, rin, rho, nu):
        """ Calculate the reynolds number in the inlet passages
        
            Bazarov Eq. 101
        """
        return 2*mdot/(np.pi*np.sqrt(n)*rin*rho*nu)
