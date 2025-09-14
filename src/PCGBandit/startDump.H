const auto* meshPtr = isA<fvMesh>(matrix().mesh());
const auto& mesh = *meshPtr;

label numDumps = numDumpsDict.getOrAdd<label>(banditName_, 0) + 1;
numDumpsDict.set<label>(banditName_, numDumps);

std::stringstream subdir;
subdir << fieldName_ << "/" << numDumps << "/";

dump<solveScalarField, IOField<solveScalar>>
(
    psi,
    subdir.str() + "init",
    mesh
);

diagDict = check_and_dump<solveScalarField, IOField<solveScalar>>
(
    banditName_,
    diagDict,
    matrix_.diag(),
    subdir.str() + "diag",
    mesh
);
