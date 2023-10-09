/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "cylinderGeomDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "globalIndex.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderGeomDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        cylinderGeomDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// A list compare binary predicate for normal sort by vector component
struct vectorLessOpCylinder
{
    typedef Foam::cylinderGeomDecomp::cylinderDir cycDir;
    const UList<vector>& values;
    const vector rotAxis_;  // the rotational axis
    vector radVec_;         // the radial vector
    cycDir sortCmpt;

    vectorLessOpCylinder(const UList<vector>& list, const vector& rotAxis, cycDir cmpt = cycDir::axial)
    :
        values(list),
        rotAxis_(rotAxis),
        sortCmpt(cmpt)
    {
        radVec_.x() = copysign(rotAxis_.z(),rotAxis_.x());
        radVec_.y() = copysign(rotAxis_.z(),rotAxis_.y());
        radVec_.z() = 
            -copysign(rotAxis_.x(),rotAxis_.z())
            -copysign(rotAxis_.y(),rotAxis_.z());
    }

    scalar copysign(const scalar& a,const scalar& b)
    {
        if (b >= 0)
            return a;
        return -1.0*a;
    }

    void setComponent(cycDir cmpt)
    {
        sortCmpt = cmpt;
    }

    bool operator()(const label a, const label b) const
    {
        switch (sortCmpt)
        {
            case cycDir::axial:
            {
                scalar axA = values[a] & rotAxis_;
                scalar axB = values[b] & rotAxis_;
                return axA < axB;
                break;
            }
            case cycDir::radial:
            {
                // subtract axial component
                vector pA = values[a];
                pA = pA - rotAxis_*(pA&rotAxis_);

                vector pB = values[b];
                pB = pB - rotAxis_*(pB&rotAxis_);

                scalar magpA = mag(pA);
                scalar magpB = mag(pB);

                return magpA < magpB;
                break;
            }
            case cycDir::circumferential:
            {
                // Get a radial vector

                // === TESTING  ===
                // Hard code the angle
                const vector& pA = values[a];
                const vector& pB = values[b];
                scalar radA = atan2(pA.x(),pA.y());
                scalar radB = atan2(pB.x(),pB.y());
                return radA < radB;
                break;
            }
        }
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// assignToProcessorGroup : given nCells cells and nProcGroup processor
// groups to share them, how do we share them out? Answer : each group
// gets nCells/nProcGroup cells, and the first few get one
// extra to make up the numbers. This should produce almost
// perfect load balancing

void Foam::cylinderGeomDecomp::assignToProcessorGroup
(
    labelList& processorGroup,
    const label nProcGroup
)
{
    label jump = processorGroup.size()/nProcGroup;
    label jumpb = jump + 1;
    label fstProcessorGroup = processorGroup.size() - jump*nProcGroup;

    label ind = 0;
    label j = 0;

    // assign cells to the first few processor groups (those with
    // one extra cell each
    for (j=0; j<fstProcessorGroup; j++)
    {
        for (label k=0; k<jumpb; k++)
        {
            processorGroup[ind++] = j;
        }
    }

    // and now to the `normal' processor groups
    for (; j<nProcGroup; j++)
    {
        for (label k=0; k<jump; k++)
        {
            processorGroup[ind++] = j;
        }
    }
}


void Foam::cylinderGeomDecomp::assignToProcessorGroup
(
    labelList& processorGroup,
    const label nProcGroup,
    const labelList& indices,
    const scalarField& weights,
    const scalar summedWeights
)
{
    // This routine gets the sorted points.
    // Easiest to explain with an example.
    // E.g. 400 points, sum of weights : 513.
    // Now with number of divisions in this direction (nProcGroup) : 4
    // gives the split at 513/4 = 128
    // So summed weight from 0..128 goes into bin 0,
    //     ,,              128..256 goes into bin 1
    //   etc.
    // Finally any remaining ones go into the last bin (3).

    const scalar jump = summedWeights/nProcGroup;
    const label nProcGroupM1 = nProcGroup - 1;
    scalar sumWeights = 0;
    label ind = 0;
    label j = 0;

    // assign cells to all except last group.
    for (j=0; j<nProcGroupM1; j++)
    {
        const scalar limit = jump*scalar(j + 1);
        while (sumWeights < limit)
        {
            sumWeights += weights[indices[ind]];
            processorGroup[ind++] = j;
        }
    }
    // Ensure last included.
    while (ind < processorGroup.size())
    {
       processorGroup[ind++] = nProcGroupM1;
    }
}


Foam::labelList Foam::cylinderGeomDecomp::decomposeOneProc
(
    const pointField& points
) const
{
    // construct a list for the final result
    labelList finalDecomp(points.size());

    labelList processorGroups(points.size());

    labelList pointIndices(identity(points.size()));

    const pointField rotatedPoints(adjustPoints(points));

    vectorLessOpCylinder sorter(rotatedPoints,rotAxis_);

    // and one to take the processor group id's. For each direction.
    // we assign the processors to groups of processors labelled
    // 0..nX to give a banded structure on the mesh. Then we
    // construct the actual processor number by treating this as
    // the units part of the processor number.


    // n_.x() is the axial component
    // n_.y() is the radial component
    // n_.z() is the circumferential component

    sorter.setComponent(cylinderDir::axial);
    Foam::sort(pointIndices, sorter);

    assignToProcessorGroup(processorGroups, n_.x());

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] = processorGroups[i];
    }


    // now do the same thing in the Y direction. These processor group
    // numbers add multiples of nX to the proc. number (columns)

    sorter.setComponent(cylinderDir::radial);

    Foam::sort(pointIndices, sorter);

    assignToProcessorGroup(processorGroups, n_.y());

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*processorGroups[i];
    }


    // finally in the Z direction. Now we add multiples of nX*nY to give
    // layers

    sorter.setComponent(cylinderDir::circumferential);
    Foam::sort(pointIndices, sorter);

    assignToProcessorGroup(processorGroups, n_.z());

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*n_.y()*processorGroups[i];
    }

    return finalDecomp;
}


Foam::labelList Foam::cylinderGeomDecomp::decomposeOneProc
(
    const pointField& points,
    const scalarField& weights
) const
{
    // construct a list for the final result
    labelList finalDecomp(points.size());

    NotImplemented;

    // labelList processorGroups(points.size());

    // labelList pointIndices(identity(points.size()));

    // const pointField rotatedPoints(adjustPoints(points));

    // vectorLessOpCylinder sorter(rotatedPoints);

    // // and one to take the processor group id's. For each direction.
    // // we assign the processors to groups of processors labelled
    // // 0..nX to give a banded structure on the mesh. Then we
    // // construct the actual processor number by treating this as
    // // the units part of the processor number.

    // sorter.setComponent(vector::X);
    // Foam::sort(pointIndices, sorter);

    // const scalar summedWeights = sum(weights);
    // assignToProcessorGroup
    // (
    //     processorGroups,
    //     n_.x(),
    //     pointIndices,
    //     weights,
    //     summedWeights
    // );

    // forAll(points, i)
    // {
    //     finalDecomp[pointIndices[i]] = processorGroups[i];
    // }


    // // now do the same thing in the Y direction. These processor group
    // // numbers add multiples of nX to the proc. number (columns)

    // sorter.setComponent(vector::Y);
    // Foam::sort(pointIndices, sorter);

    // assignToProcessorGroup
    // (
    //     processorGroups,
    //     n_.y(),
    //     pointIndices,
    //     weights,
    //     summedWeights
    // );

    // forAll(points, i)
    // {
    //     finalDecomp[pointIndices[i]] += n_.x()*processorGroups[i];
    // }


    // // finally in the Z direction. Now we add multiples of nX*nY to give
    // // layers

    // sorter.setComponent(vector::Z);
    // Foam::sort(pointIndices, sorter);

    // assignToProcessorGroup
    // (
    //     processorGroups,
    //     n_.z(),
    //     pointIndices,
    //     weights,
    //     summedWeights
    // );

    // forAll(points, i)
    // {
    //     finalDecomp[pointIndices[i]] += n_.x()*n_.y()*processorGroups[i];
    // }

    return finalDecomp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderGeomDecomp::cylinderGeomDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    geomDecomp(typeName, decompDict, regionName),
    rotAxis_(decompDict.lookupOrDefault<vector>("rotAxis",vector(0,0,1)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::cylinderGeomDecomp::decompose
(
    const pointField& points
) const
{
    if (!Pstream::parRun())
    {
        return decomposeOneProc(points);
    }
    else
    {
        const globalIndex globalNumbers(points.size());

        pointField allPoints(globalNumbers.gather(points));

        labelList allDecomp;
        if (Pstream::master())
        {
            allDecomp = decomposeOneProc(allPoints);
            allPoints.clear();  // Not needed anymore
        }

        return globalNumbers.scatter(allDecomp);
    }
}


Foam::labelList Foam::cylinderGeomDecomp::decompose
(
    const pointField& points,
    const scalarField& weights
) const
{
    if (returnReduceOr(points.size() != weights.size()))
    {
        // Ignore zero-sized weights ... and poorly sized ones too
        return decompose(points);
    }
    else if (!Pstream::parRun())
    {
        return decomposeOneProc(points, weights);
    }
    else
    {
        const globalIndex globalNumbers(points.size());

        pointField allPoints(globalNumbers.gather(points));
        scalarField allWeights(globalNumbers.gather(weights));

        labelList allDecomp;
        if (Pstream::master())
        {
            allDecomp = decomposeOneProc(allPoints, allWeights);
            allPoints.clear();  // Not needed anymore
            allWeights.clear(); // ...
        }

        return globalNumbers.scatter(allDecomp);
    }
}


// ************************************************************************* //
