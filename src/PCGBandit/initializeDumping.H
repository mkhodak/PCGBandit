#include "IOdictionary.H"
#include "objectRegistry.H"


dictionary numDumpsDict;
dictionary diagDict;
dictionary lowerAddrDict;
dictionary upperAddrDict; 
dictionary lowerDict;
dictionary A1Dict;
dictionary bDict;


template <typename DataType, typename IOType>
void dump(
    const DataType& current,
    const word fname,
    const fvMesh& mesh) {

    IOType writeCurrent(
        IOobject(
            fname,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
            ),
        current
        );
    writeCurrent.write();

}


template <typename DataType, typename IOType>
dictionary check_and_dump(
    const word key,
    dictionary dataDict,
    const DataType& current,
    const word fname,
    const fvMesh& mesh
    ) {

    bool different = false;
    if (numDumpsDict.get<label>(key) > 1) {
        DataType previous = dataDict.get<DataType>(key);
        if (previous.size() == current.size()) {
            label i;
            for (i = 0; i < current.size(); i++) {
                if (previous[i] != current[i]) {
                    different = true;
                    break;
                }
            }
            for (i = i; i < current.size(); i++) {
                previous[i] = current[i];
            }
        }
    } else {
        different = true;
        dataDict.set<DataType>(key, current);
    }

    if (different) {
        dump<DataType, IOType>(current, fname, mesh);
    }

    return dataDict;

}
