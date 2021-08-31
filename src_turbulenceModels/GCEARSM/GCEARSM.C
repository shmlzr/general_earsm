 /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "GCEARSM.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(GCEARSM, 0);
addToRunTimeSelectionTable(RASModel, GCEARSM, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> GCEARSM::F1(const volScalarField& CDkOmega) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


tmp<volScalarField> GCEARSM::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


tmp<volScalarField> GCEARSM::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*nu()/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


tmp<volScalarField> GCEARSM::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GCEARSM::GCEARSM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            coeffDict_,
            false
        )
    ),
    Theta_
    (
        coeffDict_.lookup("Theta")
    ),

    y_(mesh_),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),
    bDelta_
    (
        IOobject
        (
            "bDelta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*symm(fvc::grad(U_))/omega_
    ),
    aDelta_
    (
        IOobject
        (
            "aDelta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*bDelta_*2*k_
    ),
    Rall_
    (
        IOobject
        (
            "Rall",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*(((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)))
    ),
    aij_
    (
        IOobject
        (
            "aij",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*(((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)))
    ),
    bij_
    (
        IOobject
        (
            "bij",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*(((1.0/3.0)*I) - nut_*twoSymm(fvc::grad(U_))/(2*k_))
    ),
    Pk_
    (
        IOobject
        (
            "Pk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
        "Pk",
        dimensionSet(0,2,-3,0,0,0,0),
        0.0
        ) 
        //0.0*(nut_*2*magSqr(symm(fvc::grad(U_))) -  2.0*k_*(bDelta_&&symm(fvc::grad(U_))))
    ),
    PkDelta_
    (
        IOobject
        (
            "PkDelta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
        "PkDelta",
        dimensionSet(0,2,-3,0,0,0,0),
        0.0
        ) 
        //0.0*(nut_*2*magSqr(symm(fvc::grad(U_))) -  2.0*k_*(bDelta_&&symm(fvc::grad(U_))))
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    nut_ =
    (
        a1_*k_
      / max
        (
            a1_*omega_,
            b1_*F23()*sqrt(2.0)*mag(symm(fvc::grad(U_)))
        )
    );
    nut_.correctBoundaryConditions();
    Rall_ = ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + 2.0*k_*bDelta_;
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> GCEARSM::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + 2*k_*bDelta_,
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> GCEARSM::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_))) //linear part
           + dev(2.0*k_*bDelta_) //non-linear part
        )
    );
}


tmp<fvVectorMatrix> GCEARSM::divDevReff(volVectorField& U) const
{
    return
    
    
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U)))) //linear part
      + fvc::div(dev(2.0*k_*bDelta_)) //non-linear part
    );
    
}


tmp<fvVectorMatrix> GCEARSM::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U)))) //linear part
      + fvc::div(dev(rho*2.0*k_*bDelta_)) //non-linear part
    );
}


bool GCEARSM::read()
{
    if (RASModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        b1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        F3_.readIfPresent("F3", coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void GCEARSM::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    volTensorField gradU = fvc::grad(U_);
    volSymmTensorField Sij = dev(symm(gradU));
    volTensorField Oij = -0.5*(gradU - gradU.T());
    volScalarField S = sqrt(2*magSqr(symm(fvc::grad(U_))));
     
    volScalarField tau = 1./max( S/a1_ + omegaMin_,omega_ + omegaMin_);
    volScalarField tau2 = sqr(tau);
    volScalarField tau3 = tau*tau2;
    volScalarField tau4 = tau*tau3;
    volScalarField tau5 = tau*tau4;

    volScalarField l1 = tau2 * tr(Sij & Sij);
    volScalarField l2 = tau2 * tr(Oij & Oij);
    volScalarField l3 = tau3 * tr((Sij & Sij) & Sij);
    volScalarField l4 = tau4 * tr((Oij & Oij) & Sij);
    volScalarField l5 = tau4 * tr((Oij & Oij) & (Sij & Sij));

    volSymmTensorField T1 = tau* Sij;
    volSymmTensorField T2 = tau2 * symm((Sij & Oij) - (Oij & Sij));
    volSymmTensorField T3 = tau2 * symm(Sij & Sij) - scalar(1.0/3.0)*I*l1;
    volSymmTensorField T4 = tau2 * symm(Oij & Oij) - scalar(1.0/3.0)*I*l2;
   
    tmp<volTensorField> tgradU = fvc::grad(U_);
    
    // Include nonlinear models. 
    if (Theta_.size() == 84){
        #include "nonLinearModel.H"    
        Info << "\nSize of Theta vector " << Theta_.size() << nl << endl; 
    }
    else{
        FatalError << "The coefficient vector Theta does not contain 84 elements, but " << Theta_.size() << nl << "Change entries in constant/RASProperties" << nl << exit(FatalError);
    }
    
    bDelta_.correctBoundaryConditions();
    aDelta_ = 2*k_*bDelta_;
    aij_  = -nut_*twoSymm(fvc::grad(U_)) + 2*k_*bDelta_;
    bij_  = aij_/(2*k_);
    Rall_ = ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + 2*k_*bDelta_;


    dimensionedScalar nutSmall
    (
        "nutSmall",
        dimensionSet(0, 2, -1, 0, 0, 0 ,0),
        1e-10
    );

    
    volScalarField S2(2*magSqr(symm(tgradU())));
    
    // Production
    PkDelta_ = 2.0*k_*(bDelta_) && symm(tgradU());
    Pk_ = nut_*S2 - PkDelta_;
    volScalarField G2(GName(), Pk_);
    Info << "Production max" << max(G2) << endl;

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)
        //*min(S2, (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F23()*sqrt(S2)))
      *min(G2/(nut_+nutSmall), (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F23()*sqrt(S2)))
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G2, c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);
    

    // Re-calculate viscosity
    nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    nut_.correctBoundaryConditions();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
