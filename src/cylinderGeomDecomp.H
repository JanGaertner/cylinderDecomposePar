/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

Class
    Foam::cylinderGeomDecomp

Description
    Simple cylinder decomposition, selectable as \c cylinder

    Usage

    In the system/decomposeParDict set:

        numberOfSubdomains 64;

        method         cylinder;

        cylinderCoeffs
        {
            n       (4 4 4);    // axial, radial, circumferential
            rotAxis z;          // rotation axis
        }




    Method coefficients:
    \table
        Property  | Description                             | Required | Default
        n         | (nAxial nRadial nCircumferential)       | yes |
        rotAxis   | Select axis of rotation                 | yes |
        delta     | delta (jitter) for rotation matrix      | no  | 0.001
        transform | cartesian coordinate transformation     | no  |
    \endtable

SourceFiles
    cylinderGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderGeomDecomp_H
#define cylinderGeomDecomp_H

#include "geomDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cylinderGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class cylinderGeomDecomp
:
    public geomDecomp
{
    public:

        enum cylinderDir
        {
            radial,
            circumferential,
            axial
        };
    
    private:

        //- Rotation axis of cylinder - default z axis
        const vector rotAxis_;

    // Private Member Functions

        static void assignToProcessorGroup
        (
            labelList& processorGroup,
            const label nProcGroup
        );

        static void assignToProcessorGroup
        (
            labelList& processorGroup,
            const label nProcGroup,
            const labelList& indices,
            const scalarField& weights,
            const scalar summedWeights
        );

        labelList decomposeOneProc(const pointField& points) const;

        labelList decomposeOneProc
        (
            const pointField& points,
            const scalarField& weights
        ) const;

        //- Set the rotation axis from word or vector entry
        vector readRotAxis(const dictionary& dict);


public:

    //- No copy construct
    void operator=(const cylinderGeomDecomp&) = delete;

    //- No copy assignment
    cylinderGeomDecomp(const cylinderGeomDecomp&) = delete;

    //- Runtime type information
    TypeName("cylinder");


    // Constructors

        //- Construct given decomposition dictionary and optional region name
        explicit cylinderGeomDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~cylinderGeomDecomp() = default;


    // Member Functions

        //- Simple sends all points to the master for decomposition.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Decompose with uniform weights.
        virtual labelList decompose(const pointField& points) const;

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& weights
        ) const;

        //- Decompose with uniform weights.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points
        ) const
        {
            checkDecompositionDirections(mesh.geometricD());
            return decompose(points);
        }

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& weights
        ) const
        {
            checkDecompositionDirections(mesh.geometricD());
            return decompose(points, weights);
        }

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,  // unused
            const pointField& cc,
            const scalarField& cWeights
        ) const
        {
            return decompose(cc, cWeights);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
