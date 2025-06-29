/*---------------------------------------------------------------------------*\

Class
    Foam::PCGBandit

Group
    grpLduMatrixSolvers

Description
    Preconditioner conjugate gradient solver for symmetric lduMatrices 
    that uses Tsallis-INF (as described in Zimmert & Seldin (2021) to 
    select which preconditioner to use.

SourceFiles
    PCGBandit.C

\*---------------------------------------------------------------------------*/

#ifndef PCGBandit_H
#define PCGBandit_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class PCGBandit Declaration
\*---------------------------------------------------------------------------*/

class PCGBandit
:
    public lduMatrix::solver
{
    // Private Member Functions

        //- No copy construct
        PCGBandit(const PCGBandit&) = delete;

        //- No copy assignment
        void operator=(const PCGBandit&) = delete;

    word banditName_;
    word lossEstimator_;
    scalar learningRate_;
    bool deterministic_;
    bool randomUniform_;
    label seed_;
    label backstop_;
    label static_;

    scalar maxLogDroptol_;
    scalar minLogDroptol_;
    label numDroptols_;
    label dGAMG_;
    bool cacheAgglomeration_;

    scalar Newton(scalarField& probs, scalar x, const scalarField& k, const scalar eta) const;
    void queryTsallisINF(const scalar initialResidual) const;

    scalar lossEstimate(const scalar loss, const scalar p, const scalar t, const scalar scale) const;
    void updateTsallisINF(const scalar loss) const;

    scalar perIterationCostEstimate(const word preconditioner) const;
    scalar totalCostEstimate(const label nIterations) const;


public:

    //- Runtime type information
    TypeName("PCGBandit");

    // Constructors

        //- Construct from matrix components and solver controls
        PCGBandit
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~PCGBandit() = default;


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance scalarSolve
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt=0
        ) const;

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
