IOdictionary writeInfo
(
    IOobject
    (
        subdir.str() + "info",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    )
);
writeInfo.set<scalar>("relativeTolerance", relTol_);
writeInfo.set<scalar>("tolerance", tolerance_);
writeInfo.set<scalar>("initialResidual", solverPerf.initialResidual());
writeInfo.set<scalar>("finalResidual", solverPerf.finalResidual());
writeInfo.set<scalar>("nIterations", solverPerf.nIterations());
writeInfo.set<scalar>("solverTime", -solverTime);
writeInfo.Foam::regIOobject::write();

dump<solveScalarField, IOField<solveScalar>>
(
    psi,
    subdir.str() + "sol",
    mesh
);

lowerAddrDict = check_and_dump<List<label>, IOList<label>>
(
    banditName_,
    lowerAddrDict,
    List<label>(matrix_.lduAddr().lowerAddr()),
    subdir.str() + "lowerAddr",
    mesh
);

upperAddrDict = check_and_dump<List<label>, IOList<label>>
(
    banditName_,
    upperAddrDict,
    List<label>(matrix_.lduAddr().upperAddr()),
    subdir.str() + "upperAddr",
    mesh
);

lowerDict = check_and_dump<solveScalarField, IOField<solveScalar>>
(
    banditName_,
    lowerDict,
    matrix_.lower(),
    subdir.str() + "lower",
    mesh
);

autoPtr<solveScalarField> AxPtr(new solveScalarField(nCells, Foam::Zero));
solveScalarField& Ax = AxPtr();
matrix_.Amul(Ax, solveScalarField(nCells, Foam::One), interfaceBouCoeffs_, interfaces_, cmpt);
A1Dict = check_and_dump<solveScalarField, IOField<solveScalar>>
(
    banditName_,
    A1Dict,
    Ax,
    subdir.str() + "A1",
    mesh
); 

autoPtr<lduMatrix> identityPtr(new lduMatrix(mesh));
lduMatrix& identity = identityPtr();
identity.diag() = 1.0;
identity.lower() = 0.0;
identity.upper() = 0.0;
identity.Amul(Ax, source, interfaceBouCoeffs_, interfaces_, cmpt);
bDict = check_and_dump<solveScalarField, IOField<solveScalar>>
(
    banditName_,
    bDict,
    Ax,
    subdir.str() + "b",
    mesh
);
